# ================================================================================================================
# ================================================================================================================
#                                                    packages importation
# ================================================================================================================
# ================================================================================================================
import json
import numpy as np

import os
import regex as re
import pandas as pd

import random

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import matplotlib.cbook as cbook  # image insets

import seaborn as sns
from matplotlib.gridspec import GridSpec
from matplotlib.collections import LineCollection
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

# ================================================================================================================
# ================================================================================================================
#                                                    Open figure
# ================================================================================================================
# ================================================================================================================
NbLignes, NbColonnes=4, 4
Largeur, Hauteur = 23, 20
Largeur_ins, Hauteur_ins = 1.1, 1.1
fig1 = plt.figure(num="test_on_speed", figsize=(Largeur, Hauteur), constrained_layout=True)
gs1 = GridSpec(NbLignes, NbColonnes, figure=fig1)

## font size
plt.rcParams.update({"font.size": 15})


# other stats
ax1 = fig1.add_subplot(gs1[0, 0])
ax2 = fig1.add_subplot(gs1[0, 1])
ax3 = fig1.add_subplot(gs1[0, 2])
ax4 = fig1.add_subplot(gs1[0, 3])
ax5 = fig1.add_subplot(gs1[1, 0])
ax6 = fig1.add_subplot(gs1[1, 1])
ax7 = fig1.add_subplot(gs1[1, 2])

## add labels
def add_label(ax, txt, x=-0.05, y=1.15):
    ax.text(x, y, txt, transform=ax.transAxes, fontsize=20, fontweight="bold", va="top", ha="right")
#

#add_label(ax1_NDSC, "A", y=1.1)
#add_label(ax2_NDSC, "B", y=1.1)
#add_label(ax3_NDSC, "C")
#add_label(ax4_NDSC, "D")
#add_label(ax1, "E")
#add_label(ax2, "F")
#add_label(ax3, "G")
#add_label(ax4, "H")


# ================================================================================================================
# ================================================================================================================
#                                                    User interface
# ================================================================================================================
# ================================================================================================================

# study
substrate_regime="dyna"  # dyna or static
study="checkZspeed"
assumption="NDSC"

# assets
root=os.getcwd()
path_to_assets=os.path.join(root, substrate_regime, "assets", study)

## list of wave celerity
celerities = [0.0, 0.5]
regimes = ["Static", "NDSC", "DSC"]

# list of surfaces
surfaces=["converging_sinus", "diverging_sinus", "Snake_sinus", "wells", "rotation"]

# list of surfaces to hist
surfaces_hist=["converging_sinus", "diverging_sinus", "Snake_sinus", "wells"]
surfaces_labels=["Diverging", "Converging","Wells","Snake-like"]

## mapping of axes
#traj_mapping={("convr_sinus", 0, 1.0):ax1_NDSC,
#              ("convr_sinus", 0, -1.0):ax2_NDSC,
#              ("wells_sinus", 0, None):ax4_NDSC,
#              ("snake_sinus", 0, None):ax3_NDSC,
#              ("convr_sinus", 1, 1.0):ax1_DSC,
#              ("convr_sinus", 1, -1.0):ax2_DSC,
#              ("wells_sinus", 1, None):ax4_DSC,
#              ("snake_sinus", 1, None):ax3_DSC,
#              ("right_sinus", 0, None):ax5_NDSC
#             }

#patt_mapping={("convr_sinus", "dyna"):ax1_pattern,
#              ("wells_sinus", "dyna"):ax4_pattern,
#              ("snake_sinus", "dyna"):ax3_pattern
#              }


## color palette
rege_palette=["#8F2913", '#E1BE6A','#40B0A6']
bary_palette=["#3F2773", "#1F9E89", "#F2D027"]

# figure
show=False
save=True
path_to_save=os.path.join(os.getcwd(), substrate_regime, "test.svg")

# ================================================================================================================
# ================================================================================================================
#                                                    Functions
# ================================================================================================================
# ================================================================================================================

def searchFilesToAnalyse(root):
    paths_to_files = []

    for root_folder, sub_folder, files in os.walk(root):
        for file in files:
            if re.search('data', file):
                paths_to_files.append(os.path.join(root_folder, file))
    return paths_to_files
#

def getWaveAssumption(metadf):
    mapping={0:"NDSC", 1:"DSC"}
    return mapping[metadf["surf_assumption"]]
#

def getCurvAssumption(metadf):
    mapping={0:"0", 1:"1"}
    return mapping[metadf["curvature_assumption"]]
#

def getSign(metadf):
    surface=metadf["Surface"]
    if surface=="convr_sinus":
        return np.sign(metadf["Angular Velocity"])
    else:
        return None
#

def loadTrajectoryData(filename, list=["Surface Index", "Number of Cells", "surf_assumption", "surf_assumption"], comment="#", sep=",", eg="="):
    # Read metadata

    metadata = {}
    with open(filename, "r") as f:
        for line in f:
            if line.startswith(comment):
                if eg in line:
                    # Parse metadata in the format "comment key eg value"
                    key, value = line[1:].strip().split(eg, 1)
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
    df = pd.read_csv(filename, comment=comment, sep=sep)

    return metadata, df
#

def formatTrajectories(df, center=False, convert=True):
    
    # trajectories are centered at origin.
    if center:
        df[['x', 'y']] -= df.groupby('cell')[['x', 'y']].transform('first')

    # P5.js frame vs python frame: the y-axis is reversed
    data['y'] *= -1.0

    # convert µm to mm
    if convert:
        pos0 = data.columns.get_loc("x")
        data.insert(pos0, "X", data["x"]/1000.)
        data.insert(pos0+1, "Y", data["y"]/1000.)       
#

def plotTrajectories(metadf, df, mapping, xticks=[-.8, 0., .8], yticks=[-.8, 0., .8], N=30, centered=True):

    # spec
    if metadf["Surface"]=="right_sinus":
        xticks=[-1.7, 0., 1.3]
        yticks=[-.5, 0., 2.5]
    #
    
    ## random cells id to plot
    IDlist = np.random.choice(df['cell'].unique(), size=N, replace=False)

    # select axe to plot
    submetadf=(metadf["Surface"], metadf["surf_assumption"], getSign(metadf))  
    ax = mapping[submetadf]
    
    # extract trajectories to plot
    df_to_plot = df[df['cell'].isin(IDlist)].copy()
 
    sns.scatterplot(data=df_to_plot, 
                    x='X', 
                    y='Y', 
                    hue='cell', 
                    ax=ax, 
                    palette='deep', 
                    marker='o', 
                    edgecolor='none'
                    )
    
    ax.get_legend().set_visible(False)
    ax.set_xlabel('x (mm)', fontsize=22)
    ax.set_ylabel('y (mm)', fontsize=22)
    ax.tick_params(axis='both', which='major', labelsize=15)

    ax.set_yticks(yticks) 
    ax.set_xticks(xticks)
    ax.set_xlim(xticks[0], xticks[-1])
    ax.set_ylim(yticks[0], yticks[-1])

    # draw the initial condition either a point (if trajectories are centered) or a square
    if centered:
        ax.scatter([0.0], [0.0], s=200, color='deeppink')
    else:
        rect=mpatches.Rectangle([-0.2, -0.2], width=0.4, height=0.4, ec='deeppink', ls="--", alpha=1., fill=0, lw=4, zorder=2)
        ax.add_patch(rect)
#

def plotBarycenters(metadf, df, mapping, Timer=[10, 40, 80], plot_legend=False, xticks=[-.8, 0., .8], yticks=[-.8, 0., .8], palette=bary_palette):
    
    # select axe to plot
    submetadf=(metadf["Surface"], "dyna")  # , np.sign(metadf["Angular Velocity"])
    ax = mapping[submetadf]
    if submetadf[0]=="convr_sinus" and np.sign(metadf["Angular Velocity"])<0.0:
        ax=ax2_DSC

    pop = df[df['frame'].isin(Timer)]
    sns.kdeplot(data=pop,
                x='X', 
                y='Y',
                hue='frame',
                ax=ax, 
                palette=palette,   
                alpha=.8,
                levels=10,
                thresh=.01,
                zorder=1,
                legend=False
                )
    #

    # plot barycenter
    baryc = pop.groupby('frame').mean()
    ax.scatter(baryc['X'], baryc['Y'], marker='o', s=150, c=palette, alpha = 1., edgecolors='black', zorder=2)

    if plot_legend:
        # legend and labels
        handles, labels = plt.gca().get_legend_handles_labels()
        items = []
        for i in range(len(Timer)):
            items.append(mpatches.Patch(facecolor=palette[i], label=float(Timer[i])*dt/60., edgecolor='black'))
        #
        handles.extend(items)
        ax.legend(handles = handles, fontsize=15, title='t (min)', title_fontsize=15)
    #

    ax.set_xlabel('x (mm)', fontsize=20)
    ax.set_ylabel('y (mm)', fontsize=20) 

    ax.set_yticks(yticks) 
    ax.set_xticks(xticks)
    ax.set_xlim(xticks[0], xticks[-1])
    ax.set_ylim(yticks[0], yticks[-1])
    ax.tick_params(axis='both', which='major', labelsize=13)


def plotPattern(metadf, path, mapping):

    surface=metadf["Surface"]
    if surface=="right_sinus":
        return
    #

    
    # select axe to plot
    submetadf=(surface, "dyna") 
    ax = mapping[submetadf]
    if submetadf[0]=="convr_sinus" and np.sign(metadf["Angular Velocity"])<0.0:
        ax=ax2_pattern

    # image import
    with cbook.get_sample_data(path) as image_file:
        image = plt.imread(image_file)

    ax.imshow(image)
    ax.axis("off")  # clear x-axis and y-axis
#

def getEndingPoint(metadf, df, centered=True):

    # if trajectories has not been previously centered, it is done here
    if not centered:
        first_points = df.loc[df["frame"]==1, "Y"]
        last_points = df.loc[df["frame"]==df["frame"].max(), "Y"].reset_index(drop=True)
        Last = pd.DataFrame(last_points.sub(first_points))
    else:
        Last = df[df["time"] == df["time"].max()][["Y"]].reset_index(drop=True)
    #
    
    Last.columns=["last_y"]
    Last.insert(0, "celerity", metadf["Wave Celerity"])
    Last.insert(0, "pattern", metadf["Surface"])
    Last.insert(0, "sign", np.sign(metadf["Angular Velocity"]))
    Last.insert(0, "assumption", getAssumption(metadf))
    return Last
#

def getDx(df, globaldf):
 
    # compute Dx for each cell
    Dx_by_cell = df.groupby("cell")[["x", "y"]].apply(lambda g: (s := g.diff().abs().sum(axis=0))["x"] / s.sum())

    globaldf.insert(len(globaldf.columns), "dx", Dx_by_cell.reset_index(drop=True))
    return globaldf
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
    df.insert(pos0+13, "VzAna", -metadf["Amplitude"]*k*metadf["Wave Celerity"]*np.cos(k*df["y"] - k*metadf["Wave Celerity"]*(df["time"])))

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

def getSpeed(metadf, df, c="speed"):

    Speed = df[["time","z", "Vx", "Vy", "Vz", "VzCorr", "VzSub", "zSub", "VzAna"]].copy()
    Speed.insert(0, 'celerity', metadf['Wave Celerity'])
    Speed.insert(0, 'wavelength', int(metadf['Wavelength']))

    Speed.insert(0, "pattern", metadf["Surface"])
    Speed.insert(0, "sign", np.sign(metadf["Angular Velocity"]))
    Speed.insert(0, "Wave_assumption", getWaveAssumption(metadf))
    Speed.insert(0, "Curv_assumption", getCurvAssumption(metadf))
    return Speed
#

def prepareToPlot(df):
    # this functions rename more elegant labels for surface patterns
    labels={"convr_sinus":"Converging",
            "wells_sinus":"Wells",
            "snake_sinus": "Snake-like",
            "right_sinus": "Right Sinus"}

    df["pattern"] = df["pattern"].replace(labels)
    df.loc[(df["pattern"] == "Converging") & (df["sign"] < 0.0), "pattern"] = "Diverging"
    df.drop("sign", axis=1, inplace=True)

    # defines regime 
    labels={"NDSC":"Static"}
    pos0 = df.columns.get_loc("celerity")

    df.insert(pos0+1, "regime", df["Curv_assumption"])
    df.loc[(df["celerity"] == 0.0) & (df["Wave_assumption"]=="NDSC"), "regime"] = "Static"

    return df
#

def plotDist(metadf, df, stat, ax, plot_legend=True, marker=None):

    dt=metadf["Time Interval"]

    # PLOT
    # -------------------
    spec={"Vx":(r'$v^{obs}_x$ (µm/s)', "o"),
            "Vy":(r'$v^{obs}_y$ (µm/s)', "o"),
            "Vz":("", "o"),
            "VzSub":("", "s"),
            "VzAna":(r'$v^{obs}_z$ vs. $v^{sub,obs}_{z}$ vs. $v^{sub,th}_{z}$ (µm/s)', "x"),
            "VzCorr":(r'$v^{cell}_z$ (µm/s)', "o")
    }

    C={0.1:'steelblue',
       0.01:'orangered',
       1.0:"forestgreen"
    }

    if stat in ["VzAna"]:
        sns.lineplot(data=df,
                     x="time",
                     y=stat,
                     legend=False,
                     color=C[dt],
                     ax=ax,
                     )
    else:
        sns.scatterplot(data=df,
                    x="time",
                    y=stat,
                    legend=plot_legend,
                    color=C[dt],
                    marker=spec[stat][1],
                    ax=ax,
                    label=dt
                    )
    


    # labels fig 
    ax.set_ylabel(spec[stat][0], fontsize=17)
    ax.set_xlabel("t (s)")
    if plot_legend:
        handles, labels = plt.gca().get_legend_handles_labels()
        point1 = Line2D([0], [0], label=r'$v^{obs}_z$', marker='o', markersize=10, markerfacecolor='black', linestyle='', markeredgecolor='white')
        point2 = Line2D([0], [0], label=r'$v^{sub,obs}_z$', marker='s', markersize=10, markerfacecolor='black', linestyle='', markeredgecolor='white')
        point3 = Line2D([0], [0], label=r'$v^{sub,th}_z$', marker='', markersize=10, markerfacecolor='black', linestyle='-', color="black", markeredgecolor='white')
        mpatch1 = mpatches.Patch(facecolor='forestgreen', label=r"$\Delta t=1$s", edgecolor='black')
        mpatch2 = mpatches.Patch(facecolor='steelblue', label=r"$\Delta t=0.1$s", edgecolor='black')
        mpatch3 = mpatches.Patch(facecolor='orangered', label=r"$\Delta t=0.01$s", edgecolor='black')
        handles.extend([point1, point2, point3, mpatch1, mpatch2, mpatch3])

        ax.legend(handles=handles, fontsize=15)
        ax.tick_params(axis='both', which='major', labelsize=20)
    #

    ax.tick_params(axis='both', which='major', labelsize=20)   
#

def plotError(metadf, df, ax):
    dt=metadf["Time Interval"]


    ax.scatter([dt], [df["VzCorr"].abs().max()], marker="x", color="black", s=150)

    ax.set_ylabel(r"$max(|v^{cell}_z|)$ (µm/s)", fontsize=20)
    ax.set_xlabel(r"$\Delta t$ (s)")
    ax.set_yscale('log')
    ax.set_xscale('log')
#







    




# ================================================================================================================
# ================================================================================================================
#                                                    MAIN function
# ================================================================================================================
# ================================================================================================================

print("Plot TEST WHY CELLS GO FASTER WITH CELERITIES\n---------------------------------------\n")

list1, list2=[], []

paths=searchFilesToAnalyse(path_to_assets)

for path in paths:

    # each simulation is treated independantly
    metadata, data=loadTrajectoryData(path)
    dt=metadata["Time Interval"]
    print("\nCurrent treated data:")
    print("Sampling rate:", dt, "   Wave celerity (µm/s)", metadata["Wave Celerity"], "   Assumption", getWaveAssumption(metadata), "\n")

    # format trajectories
    print("... format trajectories")
    formatTrajectories(data, convert=True)
    
    # compute statistic
    print("... compute velocity")
    data=computeCellVelocity(metadata, data, report=False)
    summary_by_time=getSpeed(metadata, data)
    summary_by_time=prepareToPlot(summary_by_time)
    summary_by_time=summary_by_time[summary_by_time['time']<=10.]

    # plot
    plotError(metadata, summary_by_time, ax=ax5)
    if dt==0.05:
        continue

    if dt==1.:
        plotDist(metadata, summary_by_time, stat="Vx", ax=ax1)
        plotDist(metadata, summary_by_time, stat="Vy", ax=ax2)
    
    plotDist(metadata, summary_by_time, stat="Vz", ax=ax3)
    plotDist(metadata, summary_by_time, stat="VzSub", ax=ax3, plot_legend=False)
    plotDist(metadata, summary_by_time, stat="VzAna", ax=ax3)
    plotDist(metadata, summary_by_time, stat="VzCorr", ax=ax4)

    

#
#




print("\n\n----------------------------------")
if show:
    print("... figure is showed")
    plt.show()

if save:
    print("... exporting figure")
    plt.savefig(path_to_save)
    print("... figure is exported")


print("\nPlot TEST END OF PROCESS5\n---------------------------------------\n")
