# ================================================================================================================
# ================================================================================================================
#                                                    packages importation
# ================================================================================================================
# ================================================================================================================

import os
import regex as re
import json
import pandas as pd



# ================================================================================================================
# ================================================================================================================
#                                                    intermediate functions
# ================================================================================================================
# ================================================================================================================


def search_files_to_analyse(root):
    paths_to_files = []

    for root_folder, sub_folder, files in os.walk(root):
        for file in files:
            if re.search('data.csv|data.txt', file):
                paths_to_files.append(os.path.join(root_folder, file))
    return paths_to_files
#

def load_trajectory_data(filename, list):
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

def format_traj(data, center=False):
    # trajectory are centered at origin.
    if center:
        for i in data['cell'].unique():
            cell_i = data.loc[data['cell']==i,['x','y']].copy()
            data.loc[data['cell']==i,['x', 'y']] -= cell_i.iloc[0]
    #

    # js frame vs python frame: the y-axis is reversed
    data['y'] *= -1.0
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
    dt = metadata['Time Interval']

    v_list = []
    for i in data['cell'].unique():
        v_list.append(data.loc[data['cell']==i, ['x', 'y']].diff().dropna().pow(2).sum(axis=1).pow(.5).div(dt).reset_index(drop=True))
    #

    df = pd.concat(v_list, ignore_index=True)
    df = pd.DataFrame(df, columns=['vd'])
    df.insert(0, 'celerity', len(df)*[metadata['Wave Celerity']])
    df.insert(0, 'wavelength', len(df)*[int(metadata['Wavelength'])])
    porcentage = len(df[df['vd']>4.5])/len(df)*100.
    print('...porcentage of speed above 4.5µm/s: %.2f' %porcentage)
    return df
#

def get_last_point(metadata, data):
    # wave celerity
    c = metadata['Wave Celerity']

    # center trajectories
    first_points = data.loc[data['frame']==0, 'y']
    last_points = data.loc[data['frame']==data['frame'].max(), 'y'].reset_index(drop=True)

    # convert to dataframe
    df = pd.DataFrame(last_points.sub(first_points))
    df.columns=['last_y']
    df['r_d']=df['last_y'].div(c*metadata['Simulation Duration'])

    df.insert(0, 'celerity', len(df)*[c])
    df.insert(0, 'wavelength', len(df)*[int(metadata['Wavelength'])])
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
study="DSC_RightSinus_20min_VariousCelerity"

# assets
path_to_assets=os.path.join(os.getcwd(), substrate_regime, "assets", study)



# list of metadata which must be integers
must_be_int = ['Surface Index', 'Number of Cells']

## list of wave celerity to plot distributions
indiv=False
c_list = [.1, .6, 3., 10.]

## velocity mean or med
vel_cent = 'mean'

'''
## source    cell_is_surfer  /  cell_is_not_surfer  /  role_of_Vmax\\Vmax=...\\cell_is_surfer
para1='tau_dyna'
para1_value=24.0
para2='Vmax'
para2_value=6.0
assumption='cell_is_surfer'








## list of wave celerity to plot distributions
#if source=='cell_is_not_surfer'
indiv=True
c_list = [.1, .6, 3., 10.]

if assumption=='cell_is_surfer':
    source='role_of_' + para1 + '\\' + para1 + '=' + str(para1_value) + '\\' + assumption + '\\' + para2 + '=' + str(para2_value)
elif assumption=='cell_is_not_surfer':
    source = 'role_of_' + para1 + '\\' + para1 + '=' + str(para1_value) + '\\' + assumption
#
'''

# ================================================================================================================
# ================================================================================================================
#                                                    MAIN function
# ================================================================================================================
# ================================================================================================================

## look for folders containing a data.csv to analyse
files_to_analyse = search_files_to_analyse(path_to_assets)

## data analysis
dx_list, vd_list, last_list = [], [], []
for path in files_to_analyse:
    # get data and metadata
    print('\n------------------------------------ open data')
    metadata, data = load_trajectory_data(path, must_be_int)

    wl, c = metadata['Wavelength'], metadata['Wave Celerity']
    print('------------------------------------ data analysis: wavelength={}µm, celerity={}µm/s'.format(wl, c))
    # trajectories
    print('... trajectories')
    format_traj(data, center=True)


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









