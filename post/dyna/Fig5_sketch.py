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
Largeur_ins, Hauteur_ins = 1., 1.
fig1 = plt.figure(num="exotic_surfaces", figsize=(Largeur, Hauteur), constrained_layout=True)
gs1 = GridSpec(NbLignes, NbColonnes, figure=fig1)

## font size
plt.rcParams.update({"font.size": 15})

# axes for trajectories
ax1_traj = fig1.add_subplot(gs1[0, 0])
ax2_traj = fig1.add_subplot(gs1[0, 1])
ax3_traj = fig1.add_subplot(gs1[0, 2])
ax4_traj = fig1.add_subplot(gs1[0, 3])

# inset for surface pattern
ax1_pattern = inset_axes(ax1_traj, width=Largeur_ins, height=Hauteur_ins)
ax2_pattern = inset_axes(ax2_traj, width=Largeur_ins, height=Hauteur_ins)
ax3_pattern = inset_axes(ax3_traj, width=Largeur_ins, height=Hauteur_ins)
ax4_pattern = inset_axes(ax4_traj, width=Largeur_ins, height=Hauteur_ins)


ax1 = fig1.add_subplot(gs1[1, :2])
ax2 = fig1.add_subplot(gs1[1, 2:])
ax3 = fig1.add_subplot(gs1[2, :2])
#ax8 = fig1.add_subplot(gs1[1, 4])

## add labels
def add_label(ax, txt, x=-0.05, y=1.15):
    ax.text(x, y, txt, transform=ax.transAxes, fontsize=20, fontweight="bold", va="top", ha="right")
#

add_label(ax1_traj, "A", y=1.1)
add_label(ax2_traj, "B", y=1.1)
add_label(ax3_traj, "C")
add_label(ax4_traj, "D")
add_label(ax1, "E")
add_label(ax2, "F")
add_label(ax3, "G")
#add_label(ax8, "H")


# ================================================================================================================
# ================================================================================================================
#                                                    User interface
# ================================================================================================================
# ================================================================================================================

# study
substrate_regime="dyna"  # dyna or static
study="Other_surfaces"

# assets
root=os.getcwd()
path_to_assets=os.path.join(root, substrate_regime, "assets", study)

## list of wave celerity
celerities = [0.0, 0.4]

# list of surfaces
surfaces=["converging_sinus", "diverging_sinus", "Snake_sinus", "wells"]





## mapping of axes
traj_mapping={("convr_sinus", 0.4):ax1_traj,
              ("wells_sinus", 0.4):ax4_traj,
              ("snake_sinus", 0.4):ax3_traj
              }

patt_mapping={("convr_sinus", 0.4):ax1_pattern,
              ("wells_sinus", 0.4):ax4_pattern,
              ("snake_sinus", 0.4):ax3_pattern
              }


## color palette
cele_palette=["#007172", "#F29325"]


# figure
show=False
save=True
path_to_save=os.path.join(os.getcwd(), substrate_regime, "exotic_surfaces.svg")

# ================================================================================================================
# ================================================================================================================
#                                                    Functions
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

def loadTrajectoryData(filename, list=["Surface Index", "Number of Cells"], comment="#", sep=",", eg="="):
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

def formatTrajectories(df, center=True, convert=True):
    
    # trajectories are centered at origin.
    if center:
        df[['x', 'y']] -= df.groupby('cell')[['x', 'y']].transform('first')

    # P5.js frame vs python frame: the y-axis is reversed
    data['y'] *= -1.0

    # convert µm to mm
    if convert:
        data[['x', 'y']] /= 1000.0
#

def plotTrajectories(metadf, df, mapping, xticks=[-.6, 0., .6], yticks=[-.6, 0., .6], N=30, centered=True):

    ## random cells id to plot
    IDlist = np.random.choice(df['cell'].unique(), size=N, replace=False)

    # select axe to plot
    submetadf=(metadf["Surface"], metadf["Wave Celerity"])  # , np.sign(metadf["Angular Velocity"])
    ax = mapping[submetadf]
    if submetadf[0]=="convr_sinus" and np.sign(metadf["Angular Velocity"])<0.0:
        ax=ax2_traj
    
    # extract trajectories to plot
    df_to_plot = df[df['cell'].isin(IDlist)].copy()
 
    sns.scatterplot(data=df_to_plot, 
                    x='x', 
                    y='y', 
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

def plotPattern(metadf, path, mapping):

    # select axe to plot
    submetadf=(metadf["Surface"], metadf["Wave Celerity"]) 
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
        first_points = df.loc[df["frame"]==1, "y"]
        last_points = df.loc[df["frame"]==df["frame"].max(), "y"].reset_index(drop=True)
        Last = pd.DataFrame(last_points.sub(first_points))
    else:
        Last = df[df["time"] == df["time"].max()][["y"]].reset_index(drop=True)
    #

    Last.columns=["last_y"]
    Last.insert(0, "celerity", metadf["Wave Celerity"])
    Last.insert(0, "pattern", metadf["Surface"])
    Last.insert(0, "sign", np.sign(metadf["Angular Velocity"]))
    return Last
#

def getDx(df, globaldf):
 
    # compute Dx for each cell
    Dx_by_cell = df.groupby("cell")[["x", "y"]].apply(lambda g: (s := g.diff().abs().sum(axis=0))["x"] / s.sum())

    globaldf.insert(len(globaldf.columns), "dx", Dx_by_cell.reset_index(drop=True))
    return globaldf
#

def computeCellVelocity(df, dt=15., check_unity=False, report=False):
    pos0 = df.columns.get_loc("z")

    ## dx and dy
    diff_result = df.groupby("cell")[["x", "y", "z"]].diff()
    df.insert(pos0+1, "Dx", diff_result["x"])
    df.insert(pos0+2, "Dy", diff_result["y"])
    df.insert(pos0+3, "Dz", diff_result["z"])

    # vx=dx/dt and vy
    df.insert(pos0+4, "Vx", diff_result["x"]/dt)
    df.insert(pos0+5, "Vy", diff_result["y"]/dt)
    df.insert(pos0+6, "Vz", diff_result["z"]/dt)
    df = df.dropna().reset_index(drop=True)

    # speed=sqrt(vx^2+vy^2++vz^2)
    df.insert(pos0+7, "speed", (df["Vx"].pow(2) + df["Vy"].pow(2) + df["Vz"].pow(2)).pow(.5))

    # normalized velocity vx_hat=vx/speed
    df.insert(pos0+8, "Vx_hat", df["Vx"]/df["speed"])
    df.insert(pos0+9, "Vy_hat", df["Vy"]/df["speed"])
    df.insert(pos0+10, "Vz_hat", df["Vz"]/df["speed"])
    if check_unity:
        df.insert(pos0+7, "check_unity", (df["Vx_hat"].pow(2) + df["Vy_hat"].pow(2) + df["Vz_hat"].pow(2)).pow(.5))

    if report:
        print("\n-- report on velocity")
        print(df[["Vx", "Vy", "Vz", "speed"]].describe(), "\n")

    return df
#

def getSpeed(metadf, df):

    Speed=df[["speed"]].copy()

    Speed.insert(0, "celerity", metadf["Wave Celerity"])
    Speed.insert(0, "pattern", metadf["Surface"])
    Speed.insert(0, "sign", np.sign(metadf["Angular Velocity"]))
    return Speed
#

def prepareToPlot(df):
    # this functions rename more elegant labels for surface patterns
    labels={"convr_sinus":"Converging",
            "wells_sinus":"Wells",
            "snake_sinus": "Snake-like"}

    df["pattern"] = df["pattern"].replace(labels)
    df.loc[(df["pattern"] == "Converging") & (df["sign"] < 0.0), "pattern"] = "Diverging"
    df.drop("sign", axis=1, inplace=True)
    return df
#

def plotDist(df, stat, ax, palette, whis=1.5, plot_legend=True):

    # specify options
    # --------------------------
    if stat=="dx":
        infos={'y_label1':r'$D_{\perp}$'}

        ## dx=0.5 ==> no bias
        ax.axhline(y=.5, linewidth = 2, color='#8C0335', linestyle='--')
        ## dx belongs to [0, 1]
        ax.set_ylim(bottom=0., top=1.)
    #
        
    if stat=='speed':
        infos={'y_label1':'Cell speed (µm/s)'}
        ax.set_ylim(bottom=0., top=3.)
    
    if stat=="last_y":
        infos={'y_label1':r'Cell position at $T$=20min (mm)'}
        ax.axhline(y=0., linewidth = 2, color='#8C0335', linestyle='--')
        ax.set_ylim(bottom=-1., top=1.)
    #
        
    # plot distributions
    # ------------------------------
    sns.boxplot(data=df, 
                   x="pattern", 
                   y=stat, 
                   hue='celerity', 
                   linewidth=1, 
                   palette=palette, 
                   ax=ax,
                   whis = whis,
                   legend=plot_legend
                   )
    #

    # labels fig 
    ax.set_ylabel(infos['y_label1'], fontsize=25)
    ax.set_xlabel("")
    if plot_legend:
        ax.legend(fontsize=20, title="Celerity (µm/s)", title_fontsize=21, loc = "best")
    #

    ax.tick_params(axis='both', which='major', labelsize=20)    
#


    




# ================================================================================================================
# ================================================================================================================
#                                                    MAIN function
# ================================================================================================================
# ================================================================================================================

print("Plot Figure 5\n---------------------------------------\n")

list1, list2=[], []
for surface in surfaces:
    # get paths for all celerities for each surface
    working_folder=os.path.join(path_to_assets, surface)
    paths=searchFilesToAnalyse(working_folder)

    # filters the selected celerities
    paths=[path for path in paths if float(re.findall(r"c=(\d*\.?\d*)", path)[0]) in celerities]

    # each simulation is treated independantly

    for path in paths:
        metadata, data=loadTrajectoryData(path)

        print("\nCurrent treated data:")
        print("surface:", metadata["Surface"], "   Wave celerity (µm/s)", metadata["Wave Celerity"], "\n")

        # format trajectories
        print("... format trajectories")
        formatTrajectories(data)

        # plot trajectories and surface pattern
        if metadata["Wave Celerity"]==0.4:
            print("... trajectories plotting")
            plotTrajectories(metadf=metadata, df=data, mapping=traj_mapping, N=30)

            plotPattern(metadf=metadata, path=os.path.join(working_folder, "pattern.jpg"), mapping=patt_mapping)
        #

        # compute statistic
        print("... compute cell-statistics")
        summary_by_cell_local=getEndingPoint(metadata, data, centered=True)
        summary_by_cell_local=getDx(data, summary_by_cell_local)
        list1.append(summary_by_cell_local)

        print("... compute time-statistics")
        data=computeCellVelocity(data, report=True)
        list2.append(getSpeed(metadata, data))
    #
#

summary_by_cell_global=pd.concat(list1, ignore_index=True)
summary_by_time_global=pd.concat(list2, ignore_index=True)


print("\n--------------------------------------------\n")
print("\n ... plot statistics")

summary_by_cell_global=prepareToPlot(summary_by_cell_global)
summary_by_time_global=prepareToPlot(summary_by_time_global)

plotDist(summary_by_cell_global, stat="last_y", ax=ax1, palette=cele_palette, whis=1.5, plot_legend=True)
plotDist(summary_by_cell_global, stat="dx", ax=ax2, palette=cele_palette, whis=1.5, plot_legend=True)
plotDist(summary_by_time_global, stat="speed", ax=ax3, palette=cele_palette, whis=3., plot_legend=True)






        
        










        
    










print("\n\n----------------------------------")
if show:
    print("... figure is showed")
    plt.show()

if save:
    print("... exporting figure")
    plt.savefig(path_to_save)
    print("... figure is exported")


print("\nPlot Figure 5 END OF PROCESS5\n---------------------------------------\n")