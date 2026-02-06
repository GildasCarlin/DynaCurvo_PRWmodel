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

from scipy import stats
from scipy.stats import norm
import pingouin as pg
import scikit_posthocs as sp
import openturns as ot
import random

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from openturns.viewer import View
import seaborn as sns
from matplotlib.gridspec import GridSpec

# ================================================================================================================
# ================================================================================================================
#                                                    intermediate functions
# ================================================================================================================
# ================================================================================================================

def load_trajectory_data(filename, list, comment="#", sep=",", eg="="):
    # Read metadata

    metadata = {}
    with open(filename, 'r') as f:
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

def computeCellVelocity(df, dt=15., check_unity=False, report=False):
    pos0 = df.columns.get_loc("y")

    ## dx and dy
    diff_result = df.groupby("cell")[["x", "y"]].diff()
    df.insert(pos0+1, "Dx", diff_result["x"])
    df.insert(pos0+2, "Dy", diff_result["y"])

    # vx=dx/dt and vy
    df.insert(pos0+3, "Vx", diff_result["x"]/dt)
    df.insert(pos0+4, "Vy", diff_result["y"]/dt)
    df = df.dropna().reset_index(drop=True)

    # speed=sqrt(vx^2+vy^2)
    df.insert(pos0+5, "speed", (df["Vx"].pow(2) + df["Vy"].pow(2)).pow(.5))

    # normalized velocity vx_hat=vx/speed
    df.insert(pos0+6, "Vx_hat", df["Vx"]/df["speed"])
    df.insert(pos0+7, "Vy_hat", df["Vy"]/df["speed"])
    if check_unity:
        df.insert(pos0+7, "check_unity", (df["Vx_hat"].pow(2) + df["Vy_hat"].pow(2)).pow(.5))

    if report:
        print("\n-- report on velocity")
        print(df.head(), "\n")

    return df
#

def computeWaveVector(df, report=False):
    pos0 = df.columns.get_loc("theta")

    # normalized wave vector
    df.insert(pos0+1, "kx", np.cos(df["theta"]))
    df.insert(pos0+2, "ky", np.sin(df["theta"]))

    if report:
        print("\n-- report on wave vector")
        print(df.head(), "\n")
    return df
#

def filterConsecutiveSequences(group, min_length=5):
    ## This function filters for one trajectory (group) the sequence of frames that are too short
    group = group.sort_values("frame").reset_index(drop=True)
    
    # indexes the sequences and count length
    group["frame_diff"] = group["frame"].diff()
    group["new_sequence"] = (group["frame_diff"] != 1).cumsum()
    sequence_sizes = group.groupby("new_sequence").size()
    
    # isolate only long enough sequences
    valid_sequences = sequence_sizes[sequence_sizes >= min_length].index
    
    return group[group["new_sequence"].isin(valid_sequences)].drop(["frame_diff", "new_sequence"], axis=1)
#


def isolateFollowers(df, opening_angle=30., min_length=10, report=False):
    # the opening angle is expected in degrees
    threshold=np.cos(opening_angle*np.pi/180.)

    pos0 = df.columns.get_loc("ky")
    # compute scalar product between cell velocity (normalized) and wave vector (normalized)
    df.insert(pos0+1, "v_dot_k", df["Vx_hat"]*df["kx"]+df["Vy_hat"]*df["ky"])
    df.insert(pos0+2, "follower", df["v_dot_k"]>=threshold)

    # isolation of the follower only
    subdf=df[df["follower"]]
    subdf = subdf.groupby("cell", group_keys=False).apply(lambda x: filterConsecutiveSequences(x, min_length=min_length))

    if report:
        print("\n-- report on v.k")
        print("threshold (deg):", opening_angle, round(threshold, 3))
        print(df.head(), "\n")
        print(subdf.head(), "\n")

        sns.scatterplot(data=subdf, x="x", y="y", hue="cell", palette='deep', edgecolor="black", legend=False)
        plt.show()
    #

    return df, subdf
#

def normalizeCurvature(df, surface="right_sinus", a=10., L=20., report=False):
    idx={"right_sinus":.5*a*(2.*np.pi/L)**2.}
    maxH=idx[surface]

    df["meanCurvature_hat"]=df["meanCurvature"] / maxH



    if report:
        print("\n-- report on curvature")
        print("maximal mean curvature on the surface (1/µm):", round(maxH, 3))
        print(df[["meanCurvature", "meanCurvature_hat"]].describe().T, "\n")
    #
    
    return df
#

def plotCurvature(df, ax):

    # plot histogram of curvature
    sns.histplot(data=df,
                  x="meanCurvature_hat",
                  stat="density",
                  alpha=1.,
                  ax=ax,
                  zorder=10)
    #

    # background
    x_curve=np.linspace(-1., 1., 100)
    yb, yt=ax.get_ylim()
    ax.plot(x_curve,(yt-yb)*.5*(1.+np.sin(-np.pi/2.*x_curve)),
            color="black",
            linewidth=3.
            )
    #

    ax.arrow(.45, 0.9, -0.15, 0., width=.02, head_length=.04, color='dodgerblue', length_includes_head=True, transform=ax.transAxes)
    ax.arrow(.55, 0.9, 0.15, 0., width=.02, head_length=.04, color='dodgerblue', length_includes_head=True, transform=ax.transAxes)
    ax.text(0.37, 0.8,
            r"Convex",
            color='dodgerblue', rotation='horizontal', fontsize=20, fontweight='bold', va='center', ha='center', transform=ax.transAxes
            )
    ax.text(0.65, 0.8,
            r"Concave",
            color='dodgerblue', rotation='horizontal', fontsize=20, fontweight='bold', va='center', ha='center', transform=ax.transAxes
            )
    
    ax.set_xlim(left=-1., right=1.)
    ax.set_ylim(bottom=yb, top=yt)

    ax.set_xlabel(xlabel=r"Normalized $H$", fontsize=20)
    ax.set_ylabel(ylabel="Density", fontsize=20)
#

def plotCurvatureVsTime(df, ax, N=10):
    # number of traj to plot
    #N=min(N, len(df["cell"].unique()))

    # traj to plot
    #selected_id = random.sample(list(df["cell"].unique()), N)
    selected_id=[275, 26, 18, 264]
    subdf=df[df["cell"].isin(selected_id)]
    
    palette=sns.color_palette("deep", len(selected_id))

    # convert seconds to minutes
    pos0 = subdf.columns.get_loc("time")
    subdf.insert(pos0+1, "time_min", subdf["time"] / 60.)

    ## separation of sequences
    subdf = subdf.sort_values(["cell", "frame"])
    subdf["t_diff"] = subdf.groupby("cell")["frame"].diff()
    subdf["is_gap"] = (subdf["t_diff"] != 1) | subdf["t_diff"].isna()
    subdf["sequence_id"] = subdf.groupby("cell")["is_gap"].cumsum()

    #print(subdf[subdf["cell"]==275][["cell", "frame", "t_diff", "is_gap", "sequence_id"]])

    count=0
    for cell in subdf["cell"].unique():
        cell_i=subdf[subdf["cell"]==cell]
        for seq in cell_i["sequence_id"].unique():
            seq_i=cell_i[cell_i["sequence_id"]==seq]
            sns.lineplot(data=seq_i,
                         x="meanCurvature_hat",
                         y="time_min",
                         color=palette[count],
                         linewidth=5,
                         markers=True,
                         marker="o",
                         markersize=8,
                         legend=False,
                         orient="y",
                         zorder=2,
                         ax=ax
                     )
        count+=1
        #
    
    ax.set_yticks([0., 5., 10., 15., 20.])
    ax.set_xlim(left=-1., right=1.)
    ax.set_ylim(bottom=0., top=20.)
    ax.set_xlabel(xlabel=r"Normalized $H$", fontsize=20)
    ax.set_ylabel(ylabel="Time (min)", fontsize=20)

#

def plotTimeAsFollower(df, ax, report=False):
    group=df.groupby("cell")
    percentage_follower=100.*group["follower"].sum().div(group.size())
    stats=percentage_follower.describe()
    
    # plot histogram of follower
    sns.histplot(data=percentage_follower,
                  stat="density",
                  ax=ax,
                  zorder=1)
    #

    ax.axvline(x=stats["50%"], linewidth=4, color="indianred", zorder=2)

    ax.set_xlim(left=0., right=100.)
    ax.set_xticks([0., 25., stats["50%"], 50., 75., 100.])
    ax.set_xlabel("% of time spent following the wave", fontsize=15)
    ax.set_ylabel("Cell density", fontsize=20)
    
    if report:
        print("\n-- report on followers")
        print(stats.T, "\n")
    #
#

def plotVelocityVsGradient(df1, df2, axs, Nt=1):

    sns.scatterplot(data=df1[df1["follower"]==False],
                  x="meanCurvature_hat",
                  y="v_dot_e1",
                  color="lightcoral",
                  edgecolor="black",
                  alpha=.5,
                  ax=axs[0],
                  legend=False,
                  zorder=1)
    #

    sns.scatterplot(data=df2,
              x="meanCurvature_hat",
              y="v_dot_e1",
              color="cyan",
              edgecolor="black",
              alpha=.5,
              ax=axs[1],
              legend=False,
              zorder=1)
    #

    subdf=df2.groupby("cell").head(Nt)
    sns.scatterplot(data=subdf,
                    x="meanCurvature_hat",
                    y="v_dot_e1",
                    color="blueviolet",
                    edgecolor="black",
                    alpha=.5,
                    ax=axs[1],
                    legend=False,
                    zorder=2)


    for ax in axs:
        ax.set_xlim(left=-1., right=1.)
        ax.set_ylim(bottom=-1., top=1.)

        ax.set_xlabel(r"Normalized $H$", fontsize=20)
        ax.set_ylabel(r"$\hat{v}\cdot \vec{e}_1$", fontsize=20)
#

def plotBeginningOfFollowing(df, ax, Nt=1):
    subdf=df.groupby("cell").head(Nt)
    
    sns.scatterplot(data=subdf,
                    x="meanCurvature_hat",
                    y="v_dot_e1",
                    color="cyan",
                    edgecolor="black",
                    alpha=.5,
                    ax=ax,
                    legend=False,
                    zorder=1)
#

def plotVelocityGradientPhase(df, ax, N=10):
    # number of traj to plot
    N=min(N, len(df["cell"].unique()))

    # traj to plot
    selected_id = random.sample(list(df["cell"].unique()), N)
    selected_id=[275]#[275, 26, 18, 264]
    subdf=df[df["cell"].isin(selected_id)]
    
    #palette=sns.color_palette("deep", len(selected_id))

    # convert seconds to minutes
    pos0 = subdf.columns.get_loc("time")
    subdf.insert(pos0+1, "time_min", subdf["time"] / 60.)

    ### separation of sequences
    #subdf = subdf.sort_values(["cell", "frame"])
    #subdf["t_diff"] = subdf.groupby("cell")["frame"].diff()
    #subdf["is_gap"] = (subdf["t_diff"] != 1) | subdf["t_diff"].isna()
    #subdf["sequence_id"] = subdf.groupby("cell")["is_gap"].cumsum()

    #print(subdf[subdf["cell"]==275][["cell", "frame", "t_diff", "is_gap", "sequence_id"]])

    # sns.lineplot(data=subdf,
    #              x="meanCurvature_hat",
    #              y="v_dot_e1",
    #              hue="cell",
    #              palette="deep",
    #              legend=False,
    #              ax=ax,
    #              zorder=1
    #             )
    
    # sns.lineplot(data=subdf,
    #          x="meanCurvature_hat",
    #          y="v_dot_e1",
    #          hue="time",
    #          markers=True,
    #          marker="o",
    #          markersize=8,
    #          palette="viridis",
    #          legend=False,
    #          ax=ax,
    #          zorder=2
    #         )
    # Tracer les lignes en gris pour chaque trajectoire
    for idx_val in subdf["cell"].unique():
        subset = subdf[subdf["cell"] == idx_val].sort_values("frame")
        ax.plot(subset["meanCurvature_hat"], subset["v_dot_e1"], color='lightgray', linewidth=1, zorder=1)

# Superposer les points colorés selon le temps
    scatter = ax.scatter(subdf["meanCurvature_hat"], subdf["v_dot_e1"], c=subdf["frame"], 
                         cmap='viridis', s=64, zorder=2, edgecolors='white', linewidth=0.5)
    plt.colorbar(scatter, ax=ax, label="Time")
    


    
    #line = sns.lineplot(data=subdf, x="meanCurvature_hat", y="v_dot_e1", hue="frame", palette="viridis", 
    #                markers=True, marker="o", legend=False, ax=ax)


    #norm = plt.Normalize(subdf["time_min"].min(), subdf["time_min"].max())
    #sm = plt.cm.ScalarMappable(cmap="viridis", norm=norm)
    #sm.set_array([])
    #plt.colorbar(sm, ax=ax, label="Time")
    
    #sns.scatterplot(data=subdf.groupby("cell").head(1),
    #         x="meanCurvature_hat",
    #         y="v_dot_e1",
    #         color="blueviolet",
    #         s=100,
    #         legend=False,
    #         ax=ax,
    #         zorder=2
    #        )
    
    ax.set_xlim(left=-1.1, right=1.1)
    ax.set_ylim(bottom=-1.1, top=1.1)

    ax.set_xlabel(r"Normalized $H$", fontsize=20)
    ax.set_ylabel(r"$\hat{v}\cdot \vec{e}_1$", fontsize=20)
#  
    
    
#count=0
#for cell in subdf["cell"].unique():
#    cell_i=subdf[subdf["cell"]==cell]
#    for seq in cell_i["sequence_id"].unique():
#        seq_i=cell_i[cell_i["sequence_id"]==seq]
#        sns.lineplot(data=seq_i,
#                     x="meanCurvature_hat",
#                     y="time_min",
#                     color=palette[count],
#                     linewidth=5,
#                     markers=True,
#                     marker="o",
#                     markersize=8,
#                     legend=False,
#                     orient="y",
#                     zorder=2,
#                     ax=ax
#                 )
#    count+=1
#    #
    
def plotGradientAlignmentVsTime(df, ax, N=10):
    # number of traj to plot
    #N=min(N, len(df["cell"].unique()))

    # traj to plot
    #selected_id = random.sample(list(df["cell"].unique()), N)
    selected_id=[275, 26, 18, 264]
    subdf=df[df["cell"].isin(selected_id)]
    
    palette=sns.color_palette("deep", len(selected_id))

    # convert seconds to minutes
    pos0 = subdf.columns.get_loc("time")
    subdf.insert(pos0+1, "time_min", subdf["time"] / 60.)

    ## separation of sequences
    subdf = subdf.sort_values(["cell", "frame"])
    subdf["t_diff"] = subdf.groupby("cell")["frame"].diff()
    subdf["is_gap"] = (subdf["t_diff"] != 1) | subdf["t_diff"].isna()
    subdf["sequence_id"] = subdf.groupby("cell")["is_gap"].cumsum()

    #print(subdf[subdf["cell"]==275][["cell", "frame", "t_diff", "is_gap", "sequence_id"]])

    count=0
    for cell in subdf["cell"].unique():
        cell_i=subdf[subdf["cell"]==cell]
        for seq in cell_i["sequence_id"].unique():
            seq_i=cell_i[cell_i["sequence_id"]==seq]
            sns.lineplot(data=seq_i,
                         x="v_dot_e1",
                         y="time_min",
                         color=palette[count],
                         linewidth=5,
                         markers=True,
                         marker="o",
                         markersize=8,
                         legend=False,
                         orient="y",
                         zorder=2,
                         ax=ax
                     )
        count+=1
        #
    
    ax.set_yticks([0., 5., 10., 15., 20.])
    ax.set_xlim(left=-1., right=1.)
    ax.set_ylim(bottom=0., top=20.)
    ax.set_xlabel(xlabel=r"$\hat{v}\cdot \vec{e}_1$", fontsize=20)
    ax.set_ylabel(ylabel="Time (min)", fontsize=20)

#







# ================================================================================================================
# ================================================================================================================
#                                                    User interface
# ================================================================================================================
# ================================================================================================================

# study
substrate_regime="dyna"  # dyna or static
NDSC_study="test"
DSC_study="test"
# assets
root=os.getcwd()
NDSC_assets=os.path.join(root, substrate_regime, "assets", NDSC_study)
DSC_assets=os.path.join(root, substrate_regime, "assets", DSC_study)

# list of metadata which must be integers
must_be_int = ['Surface Index', 'Number of Cells']

assumptions=['NDSC', 'DSC']
mapping_path_assumption={"NDSC":NDSC_assets, "DSC":DSC_assets}

## conditions experimentales
wavelength = 20.

## list of wave celerity to plot distributions
celerities = [.7]




# ================================================================================================================
# ================================================================================================================
#                                                    Open figure
# ================================================================================================================
# ================================================================================================================
NbLignes, NbColonnes=2, 4
Largeur, Hauteur = 10, 5#30, 15
Largeur_ins, Hauteur_ins = 2.5, 2.5
fig1 = plt.figure(num='curvature_traj', figsize=(Largeur, Hauteur), constrained_layout=True)
gs1 = GridSpec(NbLignes, NbColonnes, figure=fig1)

ax1 = fig1.add_subplot(gs1[:, 0])
ax2 = fig1.add_subplot(gs1[0, 1])
ax3 = fig1.add_subplot(gs1[0, 2])
ax4 = fig1.add_subplot(gs1[:, 3])
ax5 = fig1.add_subplot(gs1[1, 1])
ax6 = fig1.add_subplot(gs1[1, 2])
#ax7 = fig1.add_subplot(gs1[1, 3])

## add labels
def add_label(ax, txt, x=-0.05, y=1.1):
    ax.text(x, y, txt, transform=ax.transAxes, fontsize=20, fontweight='bold', va='top', ha='right')
#

add_label(ax1, 'A')
add_label(ax2, 'B')
add_label(ax3, 'C')
add_label(ax4, 'D')
add_label(ax5, 'E')
add_label(ax6, 'F')
#add_label(ax7, 'G')

## font size
plt.rcParams.update({'font.size': 50})


show=True

# ================================================================================================================
# ================================================================================================================
#                                                    MAIN function
# ================================================================================================================
# ================================================================================================================

print('---------------------- open numerical data')

path=NDSC_assets + "\\" + "c=0.7" + "\\" + "data.txt"
metadata, data = load_trajectory_data(path, must_be_int)


## objectif : isoler les moments où la cellule suit la vague
##
## une cellule suit la vague si (les vecteurs sont normalisés) :
##    -- k.v>cos(a), avec a choisit préalablement
## 
## on regarde alors :
##    -- est-elle dans le concave (H>0) ? sur le convexe (H<0) ? 
##    -- sur quel pan de vague est-elle ? suit-elle la vague (a) ou alors est-elle portée par la vague (b) ?
##       (a) e1.v<0   (b) e1.v>0

## passage du repère js au repère python
data[["y", "theta"]] *= -1.

## calcul des vitesses
data=computeCellVelocity(data, dt=metadata["Time Interval"])

## calcul du vecteur d'onde
data=computeWaveVector(data)

## normalize mean curvature
data=normalizeCurvature(data,
                        surface=metadata["Surface"],
                        a=metadata["Amplitude"],
                        L=metadata["Wavelength"],
                        report=False
                        )
#

## extraction des frame où la cellule suit la vague (selon le critère vu plus haut)
data, data_followers=isolateFollowers(data, opening_angle=60., min_length=10, report=False)

##
plotTimeAsFollower(data, ax=ax2, report=False)

## 
plotCurvature(data_followers, ax=ax5)

##
plotCurvatureVsTime(data_followers, ax=ax1)

##
plotVelocityVsGradient(data, data_followers, Nt=1, axs=[ax3, ax6])

##
#plotVelocityGradientPhase(data_followers, ax=ax4, N=1)

##
plotGradientAlignmentVsTime(data_followers, ax=ax4)

# ##
# #plotBeginningOfFollowing(data_followers, Nt=1, ax=ax4)

if show:
    plt.show()








