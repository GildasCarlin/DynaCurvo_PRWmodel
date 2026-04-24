# ================================================================================================================
# ================================================================================================================
#                                                    packages importation
# ================================================================================================================
# ================================================================================================================

import os
import regex as re
import json
import pandas as pd

import numpy as np



# ================================================================================================================
# ================================================================================================================
#                                                    intermediate functions
# ================================================================================================================
# ================================================================================================================


def searchFilesToAnalyse(root):
    paths_to_files = []

    for root_folder, sub_folder, files in os.walk(root):
        for file in files:
            if re.search('data.csv|data.txt', file):
                paths_to_files.append(os.path.join(root_folder, file))
    return paths_to_files
#

def loadTrajectoryData(filename, list):
    # Read metadata

    metadata = {}
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith("#"):
                if "=" in line:
                    # Parse metadata in the format "# key=value"
                    key, value = line[1:].strip().split("=", 1)
                    key=key.strip()

                    try:
                        value=json.loads(value)
                        if key in list:                    # metadata is an integer
                            metadata[key] = value
                        else:                              # metadata is a float
                            metadata[key] = float(value)
                    except json.JSONDecodeError:           # metadata is a string
                        metadata[key] = value
            else:
                break
    #


    # Load the numerical data from the CSV file
    df = pd.read_csv(filename, comment='#')

    return metadata, df
#

def formatTrajectories(df, center=False):
    # trajectories are centered at origin.
    if center:
        for i in df['cell'].unique():
            cell_i = df.loc[df['cell']==i,['x','y']].copy()
            df.loc[df['cell']==i,['x', 'y']] -= cell_i.iloc[0]
    #

    # js frame vs python frame: the y-axis is reversed
    df['y'] *= -1.0

    return df
#

def get_Dx(metadata, data):
    Dx_by_cell = []
    for i in data['cell'].unique(): # calculates Dx for each cell
        Sx, Sy = data.loc[data['cell']==i,['x','y']].diff().abs().sum(axis=0)
        Dx_by_cell.append(Sx/(Sx+Sy))
    #

    df = pd.DataFrame(Dx_by_cell, columns=['dx'])
    df.insert(0, 'celerity', len(df)*[metadata['Wave Celerity']])
    df.insert(0, 'wavelength', len(df)*[int(metadata['Wavelength'])])
    return df
# 

def get_vd(metadata, data):

    df = data[["speed", "speedCorr", "VzSub"]].copy()
    df.insert(0, 'celerity', metadata['Wave Celerity'])
    df.insert(0, 'wavelength', int(metadata['Wavelength']))
    return df
#

def get_last_point(metadf, df):
    # wave celerity
    c = metadf['Wave Celerity']

    # center trajectories
    first_points = df.loc[df['frame']==1, 'y']
    last_points = df.loc[df['frame']==df['frame'].max(), 'y'].reset_index(drop=True)

    # convert to dataframe
    Last = pd.DataFrame(last_points.sub(first_points))
    Last.columns=['last_y']
    Last['r_d']=Last['last_y'].div(c*metadata['Simulation Duration'])

    Last.insert(0, 'celerity', c)
    Last.insert(0, 'wavelength', int(metadata['Wavelength']))
    return Last
#

def computeCellVelocity(metadf, df, check_unity=False, report=False):
    dt=metadata["Time Interval"]
    
    pos0 = df.columns.get_loc("z")

    ## FIRST: the cell velocity is computed in the external fixed frame
    ##  v=(r(t_{i+1})-r(t_i) )/Delta t
    
    ## displacement components in the external fixed frame
    diff_result = df.groupby("cell")[["x", "y", "z"]].diff()
    df.insert(pos0+1, "Dx", diff_result["x"])
    df.insert(pos0+2, "Dy", diff_result["y"])
    df.insert(pos0+3, "Dz", diff_result["z"])

    # velocity components in the external fixed frame
    df.insert(pos0+4, "Vx", diff_result["x"]/dt)
    df.insert(pos0+5, "Vy", diff_result["y"]/dt)
    df.insert(pos0+6, "Vz", diff_result["z"]/dt)
    df = df.dropna().reset_index(drop=True)

    # speed=sqrt(vx^2+vy^2++vz^2) of the cell in the external fixed frame
    df.insert(pos0+7, "speed", (df["Vx"].pow(2) + df["Vy"].pow(2) + df["Vz"].pow(2)).pow(.5))

    ## SECOND: the cell velocity is computed in the frame attached to surface but parallel to the external fixed frame
    ## v=(r(t_{i+1})-r(t_i) )/Delta t - (0, 0, h(x(t_i), t_{i+1}) - h(x(t_i), t_{i}))/Delta t

    k=2.*np.pi/metadf["Wavelength"]
    df.insert(pos0+8, "zSub", metadf["Amplitude"]*np.sin(k*df["y"] - k*metadf["Wave Celerity"]*(df["time"]+dt)))
    df.insert(pos0+9, "VzSub", (df["zSub"] - df["z"])/dt)
    df.insert(pos0+10, "VzCorr", df["Vz"] - df["VzSub"])
    df.insert(pos0+11, "speedCorr", (df["Vx"].pow(2) + df["Vy"].pow(2) + df["VzCorr"].pow(2)).pow(.5))
    df.insert(pos0+12, "speedSub", df["VzSub"].pow(2).pow(.5))

    if check_unity:
        df.insert(pos0+12, "check_unity", (df["Vx_hat"].pow(2) + df["Vy_hat"].pow(2) + df["Vz_hat"].pow(2)).pow(.5))
        # normalized velocity vx_hat=vx/speed
        df.insert(pos0+13, "Vx_hat", df["Vx"]/df["speed"])
        df.insert(pos0+14, "Vy_hat", df["Vy"]/df["speed"])
        df.insert(pos0+15, "Vz_hat", df["Vz"]/df["speed"])

    if report:
        print("\n-- report on velocity")
        print(df.head(), "\n")

    return df
#




# ================================================================================================================
# ================================================================================================================
#                                                    User interface
# ================================================================================================================
# ================================================================================================================

## ------------------------------------ folder and paths

# study
substrate_regime="dyna"  # dyna or static
study="NDSC_SpeedDistribution" # os.path.join("SpeedDistribution", "c=0.0")

# assets
path_to_assets=os.path.join(os.getcwd(), substrate_regime, "assets", study)

# list of metadata which must be integers
must_be_int = ['Surface Index', 'Number of Cells']

## list of wave celerity to plot distributions
indiv=False
c_list = [.1, .4, .6, 2.]

## velocity mean or med
vel_cent = 'mean'

# ================================================================================================================
# ================================================================================================================
#                                                    MAIN function
# ================================================================================================================
# ================================================================================================================

## look for folders containing a data.csv to analyse
files_to_analyse = searchFilesToAnalyse(path_to_assets)

## data analysis
dx_list, vd_list, last_list = [], [], []
for path in files_to_analyse:
    # get data and metadata
    print('\n------------------------------------ open data')
    metadata, data = loadTrajectoryData(path, must_be_int)

    wl, c = metadata['Wavelength'], metadata['Wave Celerity']
    print('------------------------------------ data analysis: wavelength={}µm, celerity={}µm/s'.format(wl, c))
    # trajectories
    print('... trajectories')
    data=formatTrajectories(data, center=False)
    data=computeCellVelocity(metadata, data, check_unity=False, report=False)

    if indiv and wl==20. and c in c_list: 
        print('... export')
        data.to_csv(path_to_assets + '\\' + 'c=' + str(c) + '.csv', index=False)
    #

    # computes dx trajectory bias and vd velocity distribution
    print('... dx, vd, last point')
    
    dx_list.append(get_Dx(metadata, data))
    vd_list.append(get_vd(metadata, data))
    last_list.append(get_last_point(metadata, data))
#

## concat
dx_df = pd.concat(dx_list, ignore_index=True)
vd_df = pd.concat(vd_list, ignore_index=True)
last_df = pd.concat(last_list, ignore_index=True)

dx_df = pd.DataFrame(dx_df)
vd_df = pd.DataFrame(vd_df)
last_df = pd.DataFrame(last_df)

### get the stats
central_dx = dx_df.groupby(['wavelength', 'celerity']).agg(["mean", "median", "std"]).reset_index()
central_vd = vd_df.groupby(['wavelength', 'celerity']).agg(["mean", "median", "std"]).reset_index()
central_last = last_df.groupby(['wavelength', 'celerity']).agg(["mean", "median", "std"]).reset_index()


print('\n------------------------------------ global data export: ')
dx_df.to_csv(path_to_assets + '\\' + 'dx.csv', index=False)
vd_df.to_csv(path_to_assets + '\\' + 'vd.csv', index=False)
last_df.to_csv(path_to_assets + '\\' + 'last.csv', index=False)

central_dx.to_csv(path_to_assets + '\\' + 'central_dx.csv', index=False)
central_vd.to_csv(path_to_assets + '\\' + 'central_vd.csv', index=False)
central_last.to_csv(path_to_assets + '\\' + 'central_last.csv', index=False)









