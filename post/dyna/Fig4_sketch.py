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

#from scipy import stats
#from scipy.stats import norm
#import pingouin as pg
#import scikit_posthocs as sp
#import openturns as ot
import random
from scipy.spatial.distance import cdist

import matplotlib.pyplot as plt
#import matplotlib.patches as mpatches
#from openturns.viewer import View
import seaborn as sns
from matplotlib.gridspec import GridSpec
from matplotlib.collections import LineCollection
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

# ================================================================================================================
# ================================================================================================================
#                                                    intermediate functions
# ================================================================================================================
# ================================================================================================================

def loadTrajectoryData(filename, list, comment="#", sep=",", eg="="):
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

def computeCellVelocity(df, dt=15., check_unity=False, report=False):
    pos0 = df.columns.get_loc("y")

    ## dx and dy
    diff_result = df.groupby("cell")[["x", "y"]].diff()
    df.insert(pos0+1, "Dx", diff_result["x"])
    df.insert(pos0+2, "Dy", diff_result["y"])

    # vx=dx/dt and vy
    df.insert(pos0+3, "Vx", diff_result["x"]/dt)
    df.insert(pos0+4, "Vy", diff_result["y"]/dt)
    df.insert(pos0+5, "speed", np.sqrt(np.square(df[["Vx", "Vy"]]).sum(axis=1)))
    df = df.dropna().reset_index(drop=True)

    # normalized velocity vx_hat=vx/speed
    df.insert(pos0+6, "Vx_hat", df["Vx"]/df["speed"])
    df.insert(pos0+7, "Vy_hat", df["Vy"]/df["speed"])
    #if check_unity:
    #    df.insert(pos0+7, "check_unity", (df["Vx_hat"].pow(2) + df["Vy_hat"].pow(2)).pow(.5))


    if report:
        print("\n-- report on velocity")
        print(df[["Dx", "Dy", "Vx", "Vy"]].describe(), "\n")

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

        sns.scatterplot(data=subdf, x="x", y="y", hue="cell", palette="deep", edgecolor="black", legend=False)
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

def normalizeHeight(metadata, df, report=False):
    A=metadata["Amplitude"]
    L=metadata["Wavelength"]
    k=2.*np.pi / L

    pos0 = df.columns.get_loc("z")
    # normalize height
    df.insert(pos0+1, "z_hat", df["z"] / A)

    # get coresponding x
    df.insert(pos0+2, "xbis", np.arccos(df["z_hat"]) / (k*L))

    if report:
        print("\n-- report on height")
        print("Amplitude (µm):", round(A, 3))
        print(df[["z", "z_hat", "xbis"]].describe().T, "\n")
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

    ax.arrow(.45, 0.9, -0.15, 0., width=.02, head_length=.04, color="dodgerblue", length_includes_head=True, transform=ax.transAxes)
    ax.arrow(.55, 0.9, 0.15, 0., width=.02, head_length=.04, color="dodgerblue", length_includes_head=True, transform=ax.transAxes)
    ax.text(0.37, 0.8,
            r"Convex",
            color="dodgerblue", rotation="horizontal", fontsize=20, fontweight="bold", va="center", ha="center", transform=ax.transAxes
            )
    ax.text(0.65, 0.8,
            r"Concave",
            color="dodgerblue", rotation="horizontal", fontsize=20, fontweight="bold", va="center", ha="center", transform=ax.transAxes
            )
    
    ax.set_xlim(left=-1., right=1.)
    ax.set_ylim(bottom=yb, top=yt)

    ax.set_xlabel(xlabel=r"Normalized $H$", fontsize=20)
    ax.set_ylabel(ylabel="Density", fontsize=20)
#

def plotVsTime(metadata, df, axs, N=10, rand=False):

    spec={ax1:(r"$x/\lambda$", (0., .5)), ax2:(r"$\hat{v}\cdot \vec{e}_1$", (-1., 1.))}
    # traj to plot
    selected_id=[125, 233, 185]
    if rand:
        N=min(N, len(df["cell"].unique()))
        selected_id = random.sample(list(df["cell"].unique()), N)
        print("plotted cells:", selected_id)
    subdf=df[df["cell"].isin(selected_id)]
    

    # convert seconds to minutes
    T=metadata["Simulation Duration"]/60.
    pos0 = subdf.columns.get_loc("time")
    subdf.insert(pos0+1, "time_min", subdf["time"] / 60.)

    ## separation of sequences
    subdf = subdf.sort_values(["cell", "frame"])
    subdf["t_diff"] = subdf.groupby("cell")["frame"].diff()
    subdf["is_gap"] = (subdf["t_diff"] != 1) | subdf["t_diff"].isna()
    subdf["sequence_id"] = subdf.groupby("cell")["is_gap"].cumsum()

    # plot
    palette=sns.color_palette("deep", len(selected_id))
    count=0
    for cell in subdf["cell"].unique():
        cell_i=subdf[subdf["cell"]==cell]
        for seq in cell_i["sequence_id"].unique():
            seq_i=cell_i[cell_i["sequence_id"]==seq]
            sns.lineplot(data=seq_i,
                         x="xbis",
                         y="time_min",
                         color=palette[count],
                         linewidth=1,
                         markers=True,
                         marker="o",
                         markersize=8,
                         legend=False,
                         orient="y",
                         zorder=2,
                         ax=axs[0]
                     )
            
            sns.lineplot(data=seq_i,
                         x="v_dot_e1",
                         y="time_min",
                         color=palette[count],
                         linewidth=1,
                         markers=True,
                         marker="o",
                         markersize=8,
                         legend=False,
                         orient="y",
                         zorder=2,
                         ax=axs[1]
                     )
            
            subdf1=seq_i.head(1)
            sns.scatterplot(data=subdf1,
                            x="xbis",
                            y="time_min",
                            color="black",
                            marker="d",
                            edgecolor=None,
                            s=100,
                            legend=False,
                            zorder=3,
                            ax=axs[0]
                            )
            
            sns.scatterplot(data=subdf1,
                            x="v_dot_e1",
                            y="time_min",
                            color="black",
                            marker="d",
                            edgecolor=None,
                            s=100,
                            legend=False,
                            zorder=3,
                            ax=axs[1]
                            )
        count+=1
    #
      
    for ax in axs:
        ax.set_yticks([0., .25*T, .5*T, .75*T, T])
        ax.set_xlim(spec[ax][1])
        ax.set_ylim(bottom=0., top=T)
        ax.set_xlabel(xlabel=spec[ax][0], fontsize=20)
        ax.set_ylabel(ylabel="Time (min)", fontsize=20)
    #
#

def plotTimeAsFollower(df, ax, report=False):
    group=df.groupby("cell")
    percentage_follower=100.*group["follower"].sum().div(group.size())
    stats=percentage_follower.describe()
    
    # plot histogram of follower
    sns.histplot(data=percentage_follower,
                  stat="density",
                  ax=ax,
                  facecolor="#FFE0C4",
                  zorder=1)
    #

    ax.axvline(x=stats["mean"], linewidth=4, color="mediumorchid", zorder=2)

    ax.set_xlim(left=0., right=100.)
    ax.set_xticks([0., 25., 50., 75., 100.])
    ax.set_xlabel("Time spent following the wave (%)", fontsize=15)
    ax.set_ylabel("Cell count (Density)", fontsize=20)
    

    if report:
        print("\n-- report on followers")
        print(stats.T, "\n")
    #
#

def plotAlignementCurvaturePhase(df, ax, Nt=1):

    df=df.copy().reset_index(drop=True)
    subdf_start=df.groupby("cell").head(Nt)
    subdf_left=df.drop(subdf_start.index)


    ###### plot for non-followers
    ####sns.scatterplot(data=df1[df1["follower"]==False],
    ####              x="xbis",
    ####              y="v_dot_e1",
    ####              color="lightcoral",
    ####              edgecolor="black",
    ####              alpha=.5,
    ####              ax=axs[0],
    ####              legend=False,
    ####              zorder=1)
    #####

    ## plot for followers
    sns.scatterplot(data=subdf_left,
              x="xbis",
              y="v_dot_e1",
              color="cyan",
              edgecolor="black",
              alpha=.5,
              ax=ax,
              legend=False,
              zorder=1)
    #

    ## plot for followers, when they begin to follow
    sns.scatterplot(data=subdf_start,
                    x="xbis",
                    y="v_dot_e1",
                    color="black",
                    marker="d",
                    alpha=.5,
                    ax=ax,
                    legend=False,
                    zorder=2)

    ax.set_xlim(left=0., right=.5)
    ax.set_ylim(bottom=-1., top=1.)
    ax.set_xticks([0., 0.25, 0.5])
    ax.set_yticks([-1., -.5, 0., .5, 1.])

    ax.set_xlabel(r"$x/\lambda$", fontsize=20)
    ax.set_ylabel(r"$\hat{v}\cdot \vec{e}_1$", fontsize=20)
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
        ax.plot(subset["meanCurvature_hat"], subset["v_dot_e1"], color="lightgray", linewidth=1, zorder=1)

# Superposer les points colorés selon le temps
    scatter = ax.scatter(subdf["meanCurvature_hat"], subdf["v_dot_e1"], c=subdf["frame"], 
                         cmap="viridis", s=64, zorder=2, edgecolors="white", linewidth=0.5)
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
    #ax.set_ylabel(r"$\hat{v}\cdot \vec{e}_1$", fontsize=20)
#

def plotCurvatureVsXaxis(metadata, ax):

    # surface parameters
    A=1.           #metadata["Amplitude"]
    L=1.    #metadata["Wavelength"]

    # X axis
    x = np.linspace(0., L, 200)
    k=2.*np.pi/L

    # functions definition
    def h(x):
        return A*np.cos(k*x)


    # curvature definition
    def H(x):
        hx=-A*k*np.sin(k*x)
        hxx=-A*k**2.*np.cos(k*x)
        return hxx*(1.+hx**2.)**(-1.5)
    
    def gH(x):
        hx=-A*k*np.sin(k*x)
        hxx=-A*k**2.*np.cos(k*x)
        hxxx=A*k**3.*np.sin(k*x)
        return .5*(hxxx*(1.+hx**2.)**(-1.5) - 3.*hx*hxx**2.*(1.+hx**2.)**(-2.5))

    def contr(x):
        hx=-A*k*np.sin(k*x)
        g=1./np.sqrt(1.+hx**2.0)
        return [g, hx]
    #


    # Créer les segments de ligne
    points = np.array([x, h(x)]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    
    # Créer la LineCollection avec la cmap
    norm = plt.Normalize(-1., 1.)
    lc = LineCollection(segments, cmap="viridis")
    lc.set_array(H(x)/(A*k**2.))  # Assigner les valeurs de couleur
    lc.set_linewidth(6)
    
    # Ajouter la collection
    line = ax.add_collection(lc)
    ax.autoscale()
    
    # Ajouter une colorbar
    cax = inset_axes(ax,
                     width="5%",
                     height="45%",
                     loc="upper center",
                     borderpad=1.)
    
    cbar = plt.colorbar(line, cax=cax, orientation="vertical")
    cbar.set_label(r"Normalized $H$", fontsize=15)
    cbar.set_ticks([-1, -0.5, 0, 0.5, 1])
    cax.tick_params(labelsize=10)
    
    ax.set_xlabel(r"$x/\lambda$", fontsize=20)
    ax.set_ylabel(r"$h(x)/a$", fontsize=20)

    poss=[.05, 0.1, 0.15, 0.2, 0.25, .3, .35, .4, .45, .55, .6, .65, .7, .75, .8, .85, .9, .95]
    for pos in poss:
        [g, hx]=contr(pos)

        alpha = 0.06
        norm = np.sqrt(1.**2 + hx**2)
        dx = alpha * np.sign(gH(pos)) / norm
        dy = alpha * np.sign(gH(pos))*hx / norm
        
        ax.arrow(x=pos,
                 y=0.15+h(pos),
                 dx=dx,
                 dy=dy,
                 color="indianred",
                 width=.007,
                 zorder=2,
                 length_includes_head=False,
                 head_width=.02,
                 head_length=.01)
        
    
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1])
    ax.set_yticks([-1., -.5, 0., .5, 1.])

    # arrow for k
    ax.arrow(.72, 0.2, 0.2, 0., width=.04, head_length=.08, color='grey', length_includes_head=True, transform=ax.transAxes)
    ax.text(0.8, 0.3, 
                r"$\bf k$",
                color='grey', rotation = 'horizontal', fontsize=15, fontweight='bold', va='center', ha='center', transform=ax.transAxes
                )

    import matplotlib.lines as mlines
    arrow_patch = mlines.Line2D([], [], color="indianred", marker=">", markersize=10, linestyle="-", linewidth=2)
    handles=[arrow_patch]
    ax.legend(handles=handles, labels=[r"$\bf e_1$"], loc="lower left", fontsize=14)
#

def plotSpeedInPhase(metadata, df, axs, Nt=1):
    spec={axs[0]:(r"$x/\lambda$", (0., .5)), axs[1]:(r"$\hat{v}\cdot \vec{e}_1$", (-1., 1.))}
    
    ## plot Speed Vs Curvature
    sns.scatterplot(data=df,
                x="xbis",
                y="speed",
                color="cyan",
                edgecolor="black",
                alpha=.5,
                ax=axs[0],
                legend=False,
                zorder=1)
    #
    
    ## plot Speed Vs Alignment 
    sns.scatterplot(data=df,
                    x="v_dot_e1",
                    y="speed",
                    color="cyan",
                    edgecolor="black",
                    alpha=.5,
                    ax=axs[1],
                    legend=False,
                    zorder=1)
    #

    subdf=df.groupby("cell").head(Nt)
    sns.scatterplot(data=subdf,
                    x="xbis",
                    y="speed",
                    color="black",
                    marker="d",
                    alpha=.5,
                    ax=axs[0],
                    legend=False,
                    zorder=2)
    
    sns.scatterplot(data=subdf,
                x="v_dot_e1",
                y="speed",
                color="black",
                marker="d",
                alpha=.5,
                ax=axs[1],
                legend=False,
                zorder=2)

    
    for ax in axs:
        ax.axhline(y=metadata["Wave Celerity"], color="black", linestyle="--")
        ax.set_xlim(spec[ax][1])

        ax.set_xlabel(spec[ax][0], fontsize=20)
        ax.set_ylabel(r"Cell Speed (µm/s)", fontsize=20)
    #
#

def prepareData(metadf, df):

    ## passage du repère js au repère python
    df[["y", "theta"]] *= -1.

    ## compute vx, vy and speed
    df=computeCellVelocity(df, dt=metadf["Time Interval"], report=False)

    ## compute wave vector
    df=computeWaveVector(df)

    ## normalize mean curvature
    df=normalizeCurvature(df,
                        surface=metadf["Surface"],
                        a=metadf["Amplitude"],
                        L=metadf["Wavelength"],
                        report=False
                        )
    
    ## height
    df=normalizeHeight(metadf, df, report=False)

    # sequences of followers
    df, df_followers=isolateFollowers(df,
                                      opening_angle=60.,
                                      min_length=10,
                                      report=False)
    
    return df, df_followers
#

def plotTrajectoryInPhase(ax,
                          file="dt=1.0",
                          threshold=50,
                          central_metric="median",
                          interpolate=True,                   ## trajectories are interpolated or not (on a hundred of points)
                          time_interpolation_method=True,     ## if interpolated, time is unit (between [0, 1], True) or physical time (False)
                          rand=False):                        ## a few trajectory are plotted along with the global trajectory
    
    path=NDSC_assets + "\\" + file + "\\" + "data.txt"
    metadf, df = loadTrajectoryData(path, must_be_int)
    metadf, df = prepareData(metadf, df)
    
    subdf=df[["frame", "time", "cell", "xbis", "v_dot_e1"]].copy()
    ## separation of sequences
    subdf = subdf.sort_values(["cell", "frame"])
    subdf["t_diff"] = subdf.groupby("cell")["frame"].diff()
    subdf["is_gap"] = (subdf["t_diff"] != 1) | subdf["t_diff"].isna()
    subdf["sequence_id"] = subdf.groupby("cell")["is_gap"].cumsum()
    subdf=subdf.drop(["t_diff", "is_gap"], axis=1)

    ## normalize time, each sequence begins at 0
    pos0=subdf.columns.get_loc("time")
    subdf.insert(pos0+1, "time_set0", subdf.groupby(["cell", "sequence_id"])["time"].transform(lambda x: x - x.min()))
    pos0=subdf.columns.get_loc("frame")
    subdf.insert(pos0+1, "frame_set0", subdf.groupby(["cell", "sequence_id"])["frame"].transform(lambda x: x - x.min()))

    ## normalize time, EACH sequence is between [0, 1] !! cela veut dire que lorsque t_unit=0.5 : on compare des trajs à 25s avec des trajs à 50s
    pos0=subdf.columns.get_loc("time_set0")
    subdf.insert(pos0+1, "time_unit", subdf.groupby(["cell", "sequence_id"])["time"].transform(lambda x: (x - x.min()) / (x.max() - x.min())))

    ## normalize time compared to the max Time. Sequences belong to [0, 1] but can be shorter
    pos0=subdf.columns.get_loc("time_unit")
    Tmax, Tmin=subdf["time"].max(), subdf["time"].min()
    subdf.insert(pos0+1, "time_norm", subdf.groupby(["cell", "sequence_id"])["time"].transform(lambda x: (x - Tmin) / (Tmax - Tmin)))
    
    ## for convenience, each trajectory is interpolated (1D) on a hundred of points
    if interpolate:
        n_points=100
        if time_interpolation_method in ["unit", "norm"]:
            time_bins = np.linspace(0., 1., n_points)
        else:
            max_time = subdf.groupby("cell")["time"].max().max()
            time_bins = np.linspace(0., max_time, n_points)

    
        from scipy.interpolate import interp1d
        def interpolateTrajectory(group):
            if len(group) < 2:
                return None
        
            if time_interpolation_method=="unit":
                c="time_unit"
                fill_value="extrapolate"
            if time_interpolation_method=="norm":
                c="time_norm"
                fill_value="extrapolate"
            if time_interpolation_method=="phys":
                c="time"
                fill_value=np.nan
            #
    
            group = group.sort_values(c)
            # Interpolation functions
            fA = interp1d(group[c], group["xbis"], kind="linear", bounds_error=False, fill_value=fill_value)
            fB = interp1d(group[c], group["v_dot_e1"], kind="linear", bounds_error=False, fill_value=fill_value)
    
            return pd.DataFrame({
                "time_unit": time_bins,
                "xbis": fA(time_bins),
                "v_dot_e1": fB(time_bins)
            })
        #

        # Interpolation
        interpolated = subdf.groupby(["cell", "sequence_id"]).apply(interpolateTrajectory)
        #
    #

    if central_metric in ["median", "mean"]:
        # isolate frames with enough data
        frame_to_keep=subdf.groupby("time_set0")[["xbis", "v_dot_e1"]].count()
        frame_to_keep=frame_to_keep.index[frame_to_keep["xbis"] >= threshold].tolist()
        subdf=subdf[subdf["time_set0"].isin(frame_to_keep)]
        to_plot=subdf.groupby("time_set0")[["xbis", "v_dot_e1"]].agg(["mean", "median", "std"])
        to_plot = to_plot.reset_index()

        # preparation color line
        x = to_plot[("xbis", central_metric)].values
        y = to_plot[("v_dot_e1", central_metric)].values
        time = to_plot["time_set0"].values / 60.
    #

    if central_metric in ["medoid_unit"]:
        # time is unitary (between 0 and 1) and trajectory have a constant length of 100 points (interpolated)
        
        #
        # Format 1 line = 1 cell, columns = [A(t=1), B(t=1), A(t=2), B(t=2), ...]
        trajectories_matrix = []
        idx_list = []

        for idx, group in subdf.groupby(level=["cell", "sequence_id"]):
            traj_vector = np.concatenate([group["xbis"].values, group["v_dot_e1"].values])
            trajectories_matrix.append(traj_vector)
            idx_list.append(idx)
        #
        trajectories_matrix = np.array(trajectories_matrix)

        # compute distance for each pair of point
        distances = cdist(trajectories_matrix, trajectories_matrix, metric="euclidean")

        # Find medoid
        medoid_idx = np.argmin(distances.sum(axis=1))
        medoid_id = idx_list[medoid_idx]
    
        # medoid trajectory
        medoid_traj = interpolated.loc[medoid_id].reset_index(drop=True)
        
        # preparation color line
        x = medoid_traj["xbis"].values
        y = medoid_traj["v_dot_e1"].values
        time = medoid_traj["time_unit"].values
    #
    
    if central_metric in ["medoid_dba"]:
        from tslearn.metrics import dtw
        #from tslearn.barycenters import softdtw_barycenter

        def dtw_AB(traj1, traj2):
            # traj1, traj2 : arrays shape (T1, 2) et (T2, 2)
            return dtw(traj1, traj2)  # distance DTW multivariée
        #

        def compute_medoid_dtw(trajs):
            n = len(trajs)
            dist_mat = np.zeros((n, n))

            # matrice des distances (symétrique)
            for i in range(n):
                for j in range(i + 1, n):
                    d = dtw_AB(trajs[i], trajs[j])
                    dist_mat[i, j] = d
                    dist_mat[j, i] = d

            # somme des distances pour chaque trajectoire
            sum_dists = dist_mat.sum(axis=1)
        
            # indice du médoïde
            medoid_idx = np.argmin(sum_dists)
            return medoid_idx, trajs[medoid_idx]
        #

        ## formatage
        trajectories, trajectory_ids = [], []
        for (idx, sequence), group in subdf.groupby(["cell", "sequence_id"]):
            if len(group)>30:  ## to get a long enough medoid
                traj=group[["xbis", "v_dot_e1"]].values
                trajectories.append(traj)
                trajectory_ids.append((idx, sequence))

        print("Number of trajectories:", len(trajectories))
        medoid_idx, medoid_traj = compute_medoid_dtw(trajectories) 

        ### marche pas bien (car c'est une moyenne pas un médoïde)
        #medoid_traj_test=softdtw_barycenter(trajectories[:100], max_iter=5, gamma=1.)
        #ax.plot(medoid_traj_test[:,0], medoid_traj_test[:,1], "-o", color="green", zorder=10)

        ax.plot(medoid_traj[:,0], medoid_traj[:,1], "-", color="black", linewidth=5, zorder=10)
        ax.scatter(medoid_traj[0,0], medoid_traj[0,1], marker="d", color="black", facecolor=None, s=150, zorder=11)
    #
    if central_metric not in ["medoid_dba", ""]:
        # preparation to plot
        points = np.array([x, y]).T.reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)

        norm = plt.Normalize(time.min(), time.max())
        lc = LineCollection(segments, cmap="magma", norm=norm)
        lc.set_array(time)
        lc.set_linewidth(5)
    
        line=ax.add_collection(lc)
        ax.autoscale()   
    
        cax = inset_axes(ax,
                         width="5%",
                         height="45%",
                         loc="bottom left",
                         borderpad=0.5)
        
        cbar = plt.colorbar(line, cax=cax, orientation="vertical")
        cbar.set_label(r"Time", fontsize=15)
        cax.tick_params(labelsize=10)
    
        ax.scatter(x, y, c=time, cmap="magma", s=50, zorder=10, edgecolors="black", linewidth=0.5)
    
    # Tracer quelques trajectoires individuelles (semi-transparentes)
    # 1. select them
    sample_idx=[187, 20, 140, 222, 185]
    if rand:
        sample_idx = random.sample(list(subdf["cell"].unique()), 2)
    sample_traj=subdf[subdf["cell"].isin(sample_idx)]

    # 2. prepare colorbar
    norm = plt.Normalize(sample_traj["time_set0"].min(), sample_traj["time_set0"].max())

    for idx in sample_idx:
        traj = sample_traj[sample_traj["cell"] == idx]
        traj=traj[traj["sequence_id"]==1]

        # preparation color line
        x1 = traj["xbis"].values
        y1 = traj["v_dot_e1"].values
        time1 = traj["time_set0"].values / 60.

        points1 = np.array([x1, y1]).T.reshape(-1, 1, 2)
        segments1 = np.concatenate([points1[:-1], points1[1:]], axis=1)

        #norm = plt.Normalize(time1.min(), time1.max())
        lc1 = LineCollection(segments1, cmap="magma", norm=norm, alpha=0.4)
        lc1.set_array(time1)
        lc1.set_linewidth(2)

        line=ax.add_collection(lc1)
        ax.autoscale()

        if central_metric in ["medoid_dba"] and idx==sample_idx[0]:
            cax = inset_axes(ax,
                         width="5%",
                         height="45%",
                         loc="center left",
                         borderpad=0.5)
        
            cbar = plt.colorbar(line, cax=cax, orientation="vertical")
            cbar.set_label(r"Time (min)", fontsize=15)
            cax.tick_params(labelsize=10)

        ax.scatter(x1, y1, c=time1, 
                   cmap="magma", 
                   alpha=0.4, 
                   s=20, 
                   zorder=3, 
                   edgecolors="black", 
                   linewidth=0.5)
        
        # initial condition
        start = subdf[subdf["cell"] == idx].iloc[0]
        ax.scatter(start["xbis"], start["v_dot_e1"], 
                   s=100, 
                   marker="d", 
                   alpha=.4,
                   color="black",  
                   linewidth=2)
    
    ax.set_xlim(left=0., right=0.5)
    ax.set_ylim(bottom=-1., top=1.)
    ax.set_yticks([-1., -.5, 0., .5, 1.])
    ax.set_xlabel(r"$x/\lambda$", fontsize=20)
    ax.set_ylabel(r"$\hat{v}\cdot \vec{e}_1$", fontsize=20)
    
#

def plotSingleTrajInPhase(ax, file="dt=1.0", threshold=20, rand=False, N=5):
    
    path=NDSC_assets + "\\" + file + "\\" + "data.txt"
    metadf, df = loadTrajectoryData(path, must_be_int)
    metadf, df = prepareData(metadf, df)

    subdf=df[["frame", "time", "cell", "xbis", "v_dot_e1"]].copy()
    ## separation of sequences
    subdf = subdf.sort_values(["cell", "frame"])
    subdf["t_diff"] = subdf.groupby("cell")["frame"].diff()
    subdf["is_gap"] = (subdf["t_diff"] != 1) | subdf["t_diff"].isna()
    subdf["sequence_id"] = subdf.groupby("cell")["is_gap"].cumsum()
    subdf=subdf.drop(["t_diff", "is_gap"], axis=1)

    # remove too short sequences
    subdf = subdf[subdf.groupby(["cell", "sequence_id"])["cell"].transform("size") >= threshold]

    ## normalize time, each sequence begins at 0
    pos0=subdf.columns.get_loc("time")
    subdf.insert(pos0+1, "time_set0", subdf.groupby(["cell", "sequence_id"])["time"].transform(lambda x: x - x.min()))
    pos0=subdf.columns.get_loc("frame")
    subdf.insert(pos0+1, "frame_set0", subdf.groupby(["cell", "sequence_id"])["frame"].transform(lambda x: x - x.min()))

    ## normalize time compared to the max Time. Sequences belong to [0, 1] but can be shorter
    pos0=subdf.columns.get_loc("time_set0")
    Tmax, Tmin=subdf["time_set0"].max(), subdf["time_set0"].min()
    subdf.insert(pos0+1, "time_norm", subdf.groupby(["cell", "sequence_id"])["time_set0"].transform(lambda x: (x - Tmin) / (Tmax - Tmin)))

    ## count sequence from 1
    pos0=subdf.columns.get_loc("sequence_id")
    Tmax, Tmin=subdf["time_set0"].max(), subdf["time_set0"].min()
    subdf.insert(pos0+1, "sequence_id0", subdf.groupby(["cell", "sequence_id"])["sequence_id"].transform(lambda x: x-x.min()+1))


    # filter cells and sequences index
    sample_cell_idx=[216, 63, 139, 249, 119, 209, 279, 259, 176]
    if rand:
        N=min(N, len(subdf))
        sample_cell_idx = random.sample(list(subdf["cell"].unique()), N)
        print(sample_cell_idx)
    sample_traj=subdf[subdf["cell"].isin(sample_cell_idx)]
    
    first_B_per_A = sample_traj.groupby("cell")["sequence_id"].first()
    sample_traj = sample_traj.merge(first_B_per_A.reset_index(), on=["cell", "sequence_id"])

    # initial condition
    start=sample_traj.groupby("cell").head(1)
    ax.scatter(data=start,
               x="xbis",
               y="v_dot_e1",
               marker="d",
               facecolor="black",
               alpha=0.7,
               s=100, 
               zorder=3, 
               edgecolors="black", 
               linewidth=0.5
               )
    
    for idx in sample_traj["cell"].unique():
        traj = sample_traj[sample_traj["cell"] == idx]
        x = traj["xbis"].values
        y = traj["v_dot_e1"].values
        time = traj["time_set0"].values / 60.

        points = np.array([x, y]).T.reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)

        norm = plt.Normalize(time.min(), time.max())
        lc = LineCollection(segments, cmap="magma", norm=norm, alpha=.6)
        lc.set_array(time)
        lc.set_linewidth(2)

        if idx !=sample_traj["cell"].unique()[0]:
            ax.add_collection(lc)
            ax.autoscale()
        else:
            line=ax.add_collection(lc)
            ax.autoscale()
            
            cax = inset_axes(ax, width="5%", height="45%", loc="lower left", borderpad=0.5)
            norm=plt.Normalize(sample_traj["time_set0"].min(), sample_traj["time_set0"].max())
            cbar = plt.colorbar(line, cax=cax, orientation="vertical", norm=norm)
            cbar.set_label(r"Time", fontsize=15)
            cbar.set_ticks([])
    #
    ax.set_xlim(left=0., right=0.5)
    ax.set_ylim(bottom=-1., top=1.)
    ax.set_yticks([-1., -.5, 0., .5, 1.])
    ax.set_xlabel(r"$x/\lambda$", fontsize=20)
    ax.set_ylabel(r"$\hat{v}\cdot \vec{e}_1$", fontsize=20)
#







    





# ================================================================================================================
# ================================================================================================================
#                                                    User interface
# ================================================================================================================
# ================================================================================================================

# study
substrate_regime="dyna"  # dyna or static
study="spatial_curv"
file="q=1.0"

# assets
root=os.getcwd()
NDSC_assets=os.path.join(root, substrate_regime, "assets", study)
DSC_assets=os.path.join(root, substrate_regime, "assets", study)

# list of metadata which must be integers
must_be_int = ["Surface Index", "Number of Cells"]

assumptions=["NDSC", "DSC"]
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
NbLignes, NbColonnes=3, 3
Largeur, Hauteur = 20, 20
Largeur_ins, Hauteur_ins = 1.5, 1.5
fig1 = plt.figure(num="curvature_traj", figsize=(Largeur, Hauteur), constrained_layout=False)
gs1 = GridSpec(NbLignes, NbColonnes, figure=fig1, wspace=0.5, hspace=0.3)

## font size
plt.rcParams.update({"font.size": 15})

ax1 = fig1.add_subplot(gs1[0:2, 0])
ax2 = fig1.add_subplot(gs1[0:2, 1])
ax3 = fig1.add_subplot(gs1[0, 2])
ax4 = fig1.add_subplot(gs1[1, 2])

ax5 = fig1.add_subplot(gs1[2, 0])
ax6 = fig1.add_subplot(gs1[2, 1])
ax7 = fig1.add_subplot(gs1[2, 2])
#ax8 = fig1.add_subplot(gs1[1, 4])

ax6_ins = inset_axes(ax6, width=Largeur_ins, height=Hauteur_ins, loc="upper right")

## add labels
def add_label(ax, txt, x=-0.05, y=1.15):
    ax.text(x, y, txt, transform=ax.transAxes, fontsize=20, fontweight="bold", va="top", ha="right")
#

add_label(ax1, "A", y=1.07)
add_label(ax2, "B", y=1.07)
add_label(ax3, "C")
add_label(ax4, "D")
add_label(ax5, "E")
add_label(ax6, "F")
add_label(ax7, "G")
#add_label(ax8, "H")




show=False
save=True
path_to_save=os.path.join(os.getcwd(), substrate_regime, "curv_analysis.svg")

# ================================================================================================================
# ================================================================================================================
#                                                    MAIN function
# ================================================================================================================
# ================================================================================================================

## objectif : isoler les moments où la cellule suit la vague
##
## une cellule suit la vague si (les vecteurs sont normalisés) :
##    -- k.v>cos(a), avec a choisit préalablement
## 
## on regarde alors :
##    -- est-elle dans le concave (H>0) ? sur le convexe (H<0) ? 
##    -- sur quel pan de vague est-elle ? suit-elle la vague (a) ou alors est-elle portée par la vague (b) ?
##       (a) e1.v<0   (b) e1.v>0


print("--------- open data")
path=os.path.join(NDSC_assets, file, "data.txt")

# open file
metadata, data=loadTrajectoryData(path, must_be_int)

# format_data
data, data_followers=prepareData(metadata, data)



print("--------- plot")

plotCurvatureVsXaxis(metadata, ax=ax3)
plotTimeAsFollower(data, ax=ax5, report=False)
plotVsTime(metadata, data_followers, axs=[ax1, ax2], rand=False, N=1)
#plotAlignementCurvaturePhase(data_followers, Nt=1, ax=ax7)
#plotSpeedInPhase(metadata, data_followers, axs=[ax6, ax8])
#plotTrajectoryInPhase(ax4,
#                      file="dt=1.0",
#                      threshold=50,
#                      central_metric="",                        ## medoid_unit, medoid_dba, mean, median pr empty string
#                      interpolate=False,                        ## trajectories are interpolated or not (on a hundred of points)
#                      time_interpolation_method="physical",     ## if interpolated, time is unit (between [0, 1], True) or physical time (False)
#                      rand=True)

plotSingleTrajInPhase(ax4, file="dt=1.0", rand=False, N=1)

print("--------- q analysis")

def searchFilesToAnalyse(root):
    paths_to_files = []
    for root_folder, sub_folder, files in os.walk(root):
        for file in files:
            if not re.search("dt", root_folder):
                paths_to_files.append(os.path.join(root_folder, file))
    return paths_to_files
#

files=searchFilesToAnalyse(NDSC_assets)
Mean, Count, Q, STD=[], [], [], []
Mean_perst, STD_perst=[], []
for file in files:
    # open file
    metadata, data=loadTrajectoryData(file, must_be_int)
    
    # format_data
    data, data_followers=prepareData(metadata, data)

    # get q value
    Q.append(float(re.findall(r"q=(\d*\.*\d*)", file)[0]))

    ## pourcentage time spent to follow
    group=data.groupby("cell")
    percentage_follower=100.*group["follower"].sum().div(group.size())
    stats=percentage_follower.describe()
    Mean.append(stats["mean"])
    STD.append(stats["std"])

    subdata=data_followers[["frame", "time", "cell", "xbis", "v_dot_e1"]].copy()
    ## separation of sequences
    subdata = subdata.sort_values(["cell", "frame"])
    subdata["t_diff"] = subdata.groupby("cell")["frame"].diff()
    subdata["is_gap"] = (subdata["t_diff"] != 1) | subdata["t_diff"].isna()
    subdata["sequence_id"] = subdata.groupby("cell")["is_gap"].cumsum()
    subdata=subdata.drop(["t_diff", "is_gap"], axis=1)

    # remove too short sequences
    subdata = subdata[subdata.groupby(["cell", "sequence_id"])["cell"].transform("size") >= 2]
    Perst=subdata.groupby(["cell", "sequence_id"]).size()
    Mean_perst.append(Perst.mean())
    STD_perst.append(Perst.std())
    Count.append(len(subdata.groupby(["cell", "sequence_id"]).size()))



color="mediumorchid"
Mean_perst=[15.*t for t in Mean_perst]
STD_perst=[15.*t for t in STD_perst]
ax6_ins.plot(Q, Mean_perst, marker="o", color=color)
m_plus = [Mean_perst[i] + STD_perst[i] for i in range(len(Mean_perst))]
m_minus = [Mean_perst[i] - STD_perst[i] for i in range(len(Mean_perst))]
ax6_ins.plot((Q,Q), (m_plus,m_minus),'-',color=color)

ax6_ins.set_ylim(bottom=0., top=750.)
ax6_ins.tick_params(axis="y", labelcolor=color)
ax6_ins.set_xlabel(r"$q$ $(µm^2/s)$")
ax6_ins.text(0.5, 0.75, "c=0.4", transform=ax6_ins.transAxes, fontsize=15, fontweight='bold', va='top', ha='right')

ax61_ins=ax6_ins.twinx()
ax61_ins.set_yticklabels([])
ax61_ins.set_ylim(bottom=0., top=500.)
color="steelblue"
ax61_ins.plot(Q, Count, marker="o", color=color)

ax61_ins.set_ylim(bottom=0.)
ax61_ins.tick_params(axis="y", labelcolor=color)


path_root=os.path.join(root, substrate_regime, "assets", "NDSC_CurvatureCorrection_RightSinus_20min_VariousCelerity", "L=20.0")
celerities=[0.1, 0.3, 0.4, 0.6, 0.8, 1., 1.2, 1.4, 1.6]

Mean, Count, C, STD=[], [], [], []
Mean_perst, STD_perst=[], []
for c in celerities:
    print(f"... c={c}µm/s")

    path=os.path.join(path_root, f"c={c}", "data.txt")

    metadata, data=loadTrajectoryData(path, must_be_int)

    # format_data
    data, data_followers=prepareData(metadata, data)

    # get c value
    C.append(float(re.findall(r"c=(\d*\.*\d*)", path)[0]))

    ## pourcentage time spent to follow
    group=data.groupby("cell")
    percentage_follower=100.*group["follower"].sum().div(group.size())
    stats=percentage_follower.describe()
    Mean.append(stats["mean"])
    STD.append(stats["std"])

    subdata=data_followers[["frame", "time", "cell", "xbis", "v_dot_e1"]].copy()
    ## separation of sequences
    subdata = subdata.sort_values(["cell", "frame"])
    subdata["t_diff"] = subdata.groupby("cell")["frame"].diff()
    subdata["is_gap"] = (subdata["t_diff"] != 1) | subdata["t_diff"].isna()
    subdata["sequence_id"] = subdata.groupby("cell")["is_gap"].cumsum()
    subdata=subdata.drop(["t_diff", "is_gap"], axis=1)

    # remove too short sequences
    subdata = subdata[subdata.groupby(["cell", "sequence_id"])["cell"].transform("size") >= 2]
    Perst=subdata.groupby(["cell", "sequence_id"]).size()
    Mean_perst.append(Perst.mean())
    STD_perst.append(Perst.std())
    Count.append(len(subdata.groupby(["cell", "sequence_id"]).size()))


## insert ax6
color="mediumorchid"
Mean_perst=[15.*t for t in Mean_perst]
STD_perst=[15.*t for t in STD_perst]
ax6.plot(C, Mean_perst, marker="o", color=color)
m_plus = [Mean_perst[i] + STD_perst[i] for i in range(len(Mean_perst))]
m_minus = [Mean_perst[i] - STD_perst[i] for i in range(len(Mean_perst))]
ax6.plot((C,C), (m_plus,m_minus),'-',color=color)

ax6.set_ylim(bottom=0., top=750.)
ax6.tick_params(axis="y", labelcolor=color)
ax6.set_xlabel(r"Wave celerity (µm/s)")

ax6.text(0.8, 0.7, r"$q$=1", transform=ax6_ins.transAxes, fontsize=12, fontweight='bold', va='top', ha='right')
ax6.set_ylabel("Wave-following phase duration (s)", color=color)
ax61=ax6.twinx()
ax61.set_ylim(bottom=0., top=500.)
color="steelblue"
ax61.set_ylabel("Event count", color=color)
ax61.plot(C, Count, marker="o", color=color)
ax61.set_ylim(bottom=0.)
ax61.tick_params(axis="y", labelcolor=color)

path_root=os.path.join(root, substrate_regime, "assets", "NDSC_CurvatureCorrection_RightSinus_20min_VariousCelerity", "L=40.0")
celerities=[0.1, 0.3, 0.4, 0.6, 0.8, 1., 1.2, 1.4, 1.6]

Mean, Count, C, STD=[], [], [], []
Mean_perst, STD_perst=[], []
for c in celerities:
    print(f"... c={c}µm/s")

    path=os.path.join(path_root, f"c={c}", "data.txt")

    metadata, data=loadTrajectoryData(path, must_be_int)

    # format_data
    data, data_followers=prepareData(metadata, data)

    # get c value
    C.append(float(re.findall(r"c=(\d*\.*\d*)", path)[0]))

    ## pourcentage time spent to follow
    group=data.groupby("cell")
    percentage_follower=100.*group["follower"].sum().div(group.size())
    stats=percentage_follower.describe()
    Mean.append(stats["mean"])
    STD.append(stats["std"])

    subdata=data_followers[["frame", "time", "cell", "xbis", "v_dot_e1"]].copy()
    ## separation of sequences
    subdata = subdata.sort_values(["cell", "frame"])
    subdata["t_diff"] = subdata.groupby("cell")["frame"].diff()
    subdata["is_gap"] = (subdata["t_diff"] != 1) | subdata["t_diff"].isna()
    subdata["sequence_id"] = subdata.groupby("cell")["is_gap"].cumsum()
    subdata=subdata.drop(["t_diff", "is_gap"], axis=1)

    # remove too short sequences
    subdata = subdata[subdata.groupby(["cell", "sequence_id"])["cell"].transform("size") >= 2]
    Perst=subdata.groupby(["cell", "sequence_id"]).size()
    Mean_perst.append(Perst.mean())
    STD_perst.append(Perst.std())
    Count.append(len(subdata.groupby(["cell", "sequence_id"]).size()))



## insert ax6

color="mediumorchid"
Mean_perst=[15.*t for t in Mean_perst]
STD_perst=[15.*t for t in STD_perst]
ax7.plot(C, Mean_perst, marker="o", color=color)
m_plus = [Mean_perst[i] + STD_perst[i] for i in range(len(Mean_perst))]
m_minus = [Mean_perst[i] - STD_perst[i] for i in range(len(Mean_perst))]
ax7.plot((C,C), (m_plus,m_minus),'-',color=color)
ax7.set_ylabel("Wave-following phase duration (s)", color=color)

ax7.set_ylim(bottom=0., top=750.)
ax7.tick_params(axis="y", labelcolor=color)
ax7.set_xlabel(r"Wave celerity (µm/s)")
#ax7.set_yticklabels([])
ax7.text(0.8, 0.7, r"$q$=1", transform=ax6_ins.transAxes, fontsize=12, fontweight='bold', va='top', ha='right')
ax71=ax7.twinx()
ax71.set_ylim(bottom=0., top=500.)
color="steelblue"

ax71.plot(C, Count, marker="o", color=color)
ax71.set_ylim(bottom=0.)
ax71.tick_params(axis="y", labelcolor=color)
ax71.set_ylabel("Event count", color=color)













print("\n\n----------------------------------")
if show:
    print("... figure is showed")
    plt.show()

if save:
    print("... exporting figure")
    plt.savefig(path_to_save)
    print("... figure is exported")


print("\nPlot Figure 4 END OF PROCESS5\n---------------------------------------\n")













