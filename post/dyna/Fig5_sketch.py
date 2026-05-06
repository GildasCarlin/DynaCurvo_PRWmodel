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
Largeur_ins, Hauteur_ins = 1.1, 1.1
fig1 = plt.figure(num="exotic_surfaces", figsize=(Largeur, Hauteur), constrained_layout=True)
gs1 = GridSpec(NbLignes, NbColonnes, figure=fig1)

## font size
plt.rcParams.update({"font.size": 15})

# axes for trajectories
ax1_NDSC = fig1.add_subplot(gs1[0, 0])
ax2_NDSC = fig1.add_subplot(gs1[0, 1])
ax3_NDSC = fig1.add_subplot(gs1[0, 2])
ax4_NDSC = fig1.add_subplot(gs1[0, 3])
ax5_NDSC = fig1.add_subplot(gs1[2, 0])

# inset for surface pattern
ax1_pattern = inset_axes(ax1_NDSC, width=Largeur_ins, height=Hauteur_ins)
ax2_pattern = inset_axes(ax2_NDSC, width=Largeur_ins, height=Hauteur_ins)
ax3_pattern = inset_axes(ax3_NDSC, width=Largeur_ins, height=Hauteur_ins)
ax4_pattern = inset_axes(ax4_NDSC, width=Largeur_ins, height=Hauteur_ins)

# axes for barycenters
ax1_DSC = fig1.add_subplot(gs1[1, 0])
ax2_DSC = fig1.add_subplot(gs1[1, 1])
ax3_DSC = fig1.add_subplot(gs1[1, 2])
ax4_DSC = fig1.add_subplot(gs1[1, 3])

# other stats
ax1 = fig1.add_subplot(gs1[2, 1])
ax2 = fig1.add_subplot(gs1[2, 2])
ax3 = fig1.add_subplot(gs1[2, 3])
#ax4 = fig1.add_subplot(gs1[3, 2:])

## add labels
def add_label(ax, txt, x=-0.05, y=1.15):
    ax.text(x, y, txt, transform=ax.transAxes, fontsize=20, fontweight="bold", va="top", ha="right")
#

add_label(ax1_NDSC, "A.1")
add_label(ax2_NDSC, "A.2")
add_label(ax3_NDSC, "A.3")
add_label(ax4_NDSC, "A.4")
add_label(ax1_DSC, "B.1")
add_label(ax2_DSC, "B.2")
add_label(ax3_DSC, "B.3")
add_label(ax4_DSC, "B.4")
add_label(ax1, "D")
add_label(ax2, "E")
add_label(ax3, "F")
add_label(ax5_NDSC, "C")


# ================================================================================================================
# ================================================================================================================
#                                                    User interface
# ================================================================================================================
# ================================================================================================================

# study
substrate_regime="dyna"  # dyna or static
study="test"

# assets
root=os.getcwd()
path_to_assets=os.path.join(root, substrate_regime, "assets", study)

## list of wave celerity
celerities = [0.0, 0.5]
regimes = ["Static", "DSC(-)", "DSC(+)"]

# list of surfaces
surfaces=["converging_sinus", "diverging_sinus", "Snake_sinus", "wells", "rotation"]

# list of surfaces to hist
surfaces_hist=["converging_sinus", "diverging_sinus", "Snake_sinus", "wells"]
surfaces_labels=["Diverging", "Converging","Wells","Snake-like"]

## mapping of axes
traj_mapping={("convr_sinus", 0, 1.0):ax1_NDSC,
              ("convr_sinus", 0, -1.0):ax2_NDSC,
              ("wells_sinus", 0, None):ax4_NDSC,
              ("snake_sinus", 0, None):ax3_NDSC,
              ("convr_sinus", 1, 1.0):ax1_DSC,
              ("convr_sinus", 1, -1.0):ax2_DSC,
              ("wells_sinus", 1, None):ax4_DSC,
              ("snake_sinus", 1, None):ax3_DSC,
              ("right_sinus", 0, None):ax5_NDSC
             }

patt_mapping={("convr_sinus", "dyna"):ax1_pattern,
              ("wells_sinus", "dyna"):ax4_pattern,
              ("snake_sinus", "dyna"):ax3_pattern
              }


## color palette
rege_palette=["#8F2913", '#E1BE6A','#40B0A6']
bary_palette=["#3F2773", "#1F9E89", "#F2D027"]

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

def getAssumption(metadf):
    mapping={0:"DSC(-)", 1:"DSC(+)"}
    return mapping[metadf["surf_assumption"]]
#

def getSign(metadf):
    surface=metadf["Surface"]
    if surface=="convr_sinus":
        return np.sign(metadf["Angular Velocity"])
    else:
        return None
#

def loadTrajectoryData(filename, list=["Surface Index", "Number of Cells", "surf_assumption"], comment="#", sep=",", eg="="):
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
    df['y'] *= -1.0

    # convert µm to mm
    if convert:
        pos0 = df.columns.get_loc("x")
        df.insert(pos0, "X", df["x"]/1000.)
        df.insert(pos0+1, "Y", df["y"]/1000.)    

    return df   
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
    
    if ax == ax1_NDSC:
        ax.text(-.25, .5, 
                r"DSC(-) cells", 
                rotation = 'vertical', fontsize=20, fontweight='bold', va='center', ha='center', transform=ax.transAxes
                )
        #
        
        ax.arrow(0.1, .3, 0., 0.4, width=.04, head_length=.1, color='red', transform=ax.transAxes, length_includes_head=True)
        ax.text(0.18, 0.5, 
                r"Wave propagation",
                color='red', rotation = 'vertical', fontsize=15, fontweight='bold', va='center', ha='center', transform=ax.transAxes
                )
        #


    if ax == ax1_DSC:
        ax.text(-.25, .5, 
                r"DSC(+) cells", 
                rotation = 'vertical', fontsize=20, fontweight='bold', va='center', ha='center', transform=ax.transAxes
                )
        #

        ax.arrow(0.1, .3, 0., 0.4, width=.04, head_length=.1, color='red', transform=ax.transAxes, length_includes_head=True)
        ax.text(0.18, 0.5, 
                r"Wave propagation",
                color='red', rotation = 'vertical', fontsize=15, fontweight='bold', va='center', ha='center', transform=ax.transAxes
                )
        #
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
    print(df.head())
    dt=metadf["Time Interval"]
    
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

    print(df.head())

    ## SECOND: the cell velocity is computed in the frame attached to surface but parallel to the external fixed frame
    ## v=(r(t_{i+1})-r(t_i) )/Delta t - (0, 0, h(x(t_i), t_{i+1}) - h(x(t_i), t_{i}))/Delta t

    df.insert(pos0+8, "zSub", computeZPosition(metadf, df))
    df.insert(pos0+9, "VzSub", (df["zSub"].copy() - df["z"]).copy()/dt)
    df.insert(pos0+10, "VzCorr", df["Vz"].copy() - df["VzSub"].copy())
    df.insert(pos0+11, "speedCorr", (df["Vx"].pow(2) + df["Vy"].pow(2) + df["VzCorr"].pow(2)).pow(.5)) ## cell speed corrected 
    df.insert(pos0+12, "speedSub", df["VzSub"].pow(2).pow(.5))

    if check_unity:
        df.insert(pos0+12, "check_unity", (df["Vx_hat"].pow(2) + df["Vy_hat"].pow(2) + df["Vz_hat"].pow(2)).pow(.5))
        # normalized velocity vx_hat=vx/speed
        df.insert(pos0+13, "Vx_hat", df["Vx"]/df["speed"])
        df.insert(pos0+14, "Vy_hat", df["Vy"]/df["speed"])
        df.insert(pos0+15, "Vz_hat", df["Vz"]/df["speed"])

    if report:
        print("\n-- report on velocity")
        print(df[["Dx", "Dy", "Dz", "Vx", "Vy", "Vz", "z", "zSub", "speedCorr"]].head(), "\n")

    return df
#

def computeZPosition(metadf, df):

    k=2.*np.pi/metadf["Wavelength"]

    if metadf["Surface"]=="right_sinus":
        return metadf["Amplitude"]*np.sin(k*df["y"].copy() - k*metadf["Wave Celerity"]*(df["time"]+metadf["Time Interval"]))
    
    if metadf["Surface"]=="wells_sinus":
    



def getSpeed(metadf, df, cl="speedCorr"):

    Speed=df[[cl]].copy()

    Speed.insert(0, "celerity", metadf["Wave Celerity"])
    Speed.insert(0, "pattern", metadf["Surface"])
    Speed.insert(0, "sign", np.sign(metadf["Angular Velocity"]))
    Speed.insert(0, "assumption", getAssumption(metadf))
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
    labels={"DSC(-)":"Static"}
    pos0 = df.columns.get_loc("celerity")

    df.insert(pos0+1, "regime", df["assumption"])
    df.loc[(df["celerity"] == 0.0) & (df["assumption"]=="DSC(-)"), "regime"] = "Static"

    return df
#

def plotDist(df, stat, ax, palette=rege_palette, whis=1.5, plot_legend=True, hue_order=regimes):

    # specify options
    # --------------------------
    if stat=="dx":
        infos={'y_label1':r'$D_{\mathbf{k}}$'}

        ## dx=0.5 ==> no bias
        ax.axhline(y=.5, linewidth = 2, color='#8C0335', linestyle='--')
        ## dx belongs to [0, 1]
        ax.set_ylim(bottom=0., top=1.)
    #
        
    if stat=='speedCorr':
        infos={'y_label1':r'Cell speed $|| \mathbf{v} ||_2$ (µm/s)'}
        #ax.set_ylim(bottom=0., top=3.)

    if stat=='speed2D':
        infos={'y_label1':'Apparent cell speed (µm/s)'}
        #ax.set_ylim(bottom=0., top=3.)
    
    if stat=="last_y":
        infos={'y_label1':r"$B_T$"+" at "+r'$T$=20min (mm)'}
        ax.axhline(y=0., linewidth = 2, color='#8C0335', linestyle='--')
        ax.set_ylim(bottom=-1., top=1.)
    #

    if stat=="distance":
        infos={'y_label1':r'Apperent distance travelled (µm)'}
        ax.axhline(y=0., linewidth = 2, color='#8C0335', linestyle='--')
        ax.set_ylim(bottom=-1., top=1.)
    #
        
    # plot distributions
    # ------------------------------
    sns.boxplot(data=df, 
                   x="pattern", 
                   y=stat, 
                   hue='regime', 
                   hue_order=hue_order,
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
        ax.legend(fontsize=15, loc = "best")
    #

    ax.tick_params(axis='both', which='major', labelsize=20)   
    ax.set_xticklabels(labels=df["pattern"].unique(), rotation=45)
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
    paths_Zspeed=[path for path in paths if "Zspeed" in path]
    paths=[path for path in paths if "Zspeed" not in path]


    # each simulation is treated independantly
    for path in paths:
        metadata, data=loadTrajectoryData(path)
        surf, c=metadata["Surface"], metadata["Wave Celerity"]

        print("\n----------------------------------------------------------")
        print("Current treated data:")
        print("surface:", surf, "   Wave celerity (µm/s)", c, "   Assumption", getAssumption(metadata), "\n")

        # format trajectories
        print("... format trajectories")
        data=formatTrajectories(data, convert=True)

        # plot trajectories and surface pattern
        if c in [.5]:
            print("... trajectories plotting")
            plotTrajectories(metadf=metadata, df=data, mapping=traj_mapping, N=30)
            plotPattern(metadf=metadata, path=os.path.join(working_folder, "pattern.jpg"), mapping=patt_mapping)
            #plotBarycenters(metadata, data, mapping=DSC_bary_mapping)
        #

        # compute statistic
        print("... compute cell-statistics")
        summary_by_cell_local=getEndingPoint(metadata, data, centered=True)
        summary_by_cell_local=getDx(data, summary_by_cell_local)
        list1.append(summary_by_cell_local)

        if surf not in ["right_sinus"]:
            # open file for the speed computation
            print("... format trajectories --- Zspeed")
            metadataZ, dataZ=loadTrajectoryData(re.sub(rf"c={c}", rf"c={c}_Zspeed", path))
            dataZ=formatTrajectories(dataZ, convert=False, center=False)
        
            print("... compute time-statistics")
            dataZ=computeCellVelocity(metadataZ, dataZ, report=False)
            summary_by_time_local=getSpeed(metadataZ, dataZ, cl="speedCorr")
            list2.append(summary_by_time_local)
    #
#

summary_by_cell_global=pd.concat(list1, ignore_index=True)
summary_by_time_global=pd.concat(list2, ignore_index=True)


print("\n\n--------------------------------------------")
print("... plot statistics")

# format to plot
summary_by_cell_global=prepareToPlot(summary_by_cell_global)
summary_by_time_global=prepareToPlot(summary_by_time_global)

# filter assumptions
summary_by_cell_to_plot=summary_by_cell_global[summary_by_cell_global["assumption"].isin(regimes)]
summary_by_time_to_plot=summary_by_time_global[summary_by_time_global["assumption"].isin(regimes)]

summary_by_cell_to_plot=summary_by_cell_to_plot[summary_by_cell_to_plot["pattern"].isin(surfaces_labels)]
summary_by_time_to_plot=summary_by_time_to_plot[summary_by_time_to_plot["pattern"].isin(surfaces_labels)]

# plots
plotDist(summary_by_cell_to_plot, stat="last_y", ax=ax1, whis=1.5, plot_legend=True)
plotDist(summary_by_cell_to_plot, stat="dx", ax=ax2, whis=1.5, plot_legend=True)
plotDist(summary_by_time_to_plot, stat="speedCorr", ax=ax3, whis=1.5, plot_legend=True) # .copy().sample(frac=1.)







        
        










        
    










print("\n\n----------------------------------")
if show:
    print("... figure is showed")
    plt.show()

if save:
    print("... exporting figure")
    plt.savefig(path_to_save)
    print("... figure is exported")


print("\nPlot Figure 5 END OF PROCESS5\n---------------------------------------\n")