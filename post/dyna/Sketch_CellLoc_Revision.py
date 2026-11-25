# ================================================================================================================
# ================================================================================================================
#                                                    packages importation
# ================================================================================================================
# ================================================================================================================

# gestion
import json
import os

# data et maths
import pandas as pd
import numpy as np
from scipy import stats
import random

# images
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.gridspec import GridSpec
from matplotlib.collections import LineCollection
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


# ================================================================================================================
# ================================================================================================================
#                                                    intermediate functions
# ================================================================================================================
# ================================================================================================================

#palette
no_follow="#B8D48A"
follow="#057428"
force="#DE5217"


def loadTrajectoryData(filename, comment="#", sep=",", eg="="):
    must_be_int = ["Surface Index", "Number of Cells"]
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
                        if key in must_be_int:                    # metadata is an integer
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

def computeCellVelocity(df, dt=15., report=False):
    # computes cell planar velocity
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

def isolateFollowers(df, opening_angle=60., min_length=10, report=False):
    # the opening angle is expected in degrees
    threshold=np.cos(opening_angle*np.pi/180.)

    pos0 = df.columns.get_loc("ky")
    # compute scalar product between cell velocity (normalized) and wave vector (normalized)
    df.insert(pos0+1, "v_dot_k", df["Vx_hat"]*df["kx"]+df["Vy_hat"]*df["ky"])
    df.insert(pos0+2, "follower", df["v_dot_k"]>=threshold)
    
    # the cell location is positionned on the wave (x axis)
    df.insert(pos0+3, "x_norm", (1./(2.*np.pi))*np.arccos(df["z_hat"]))
    df.loc[df["v_dot_e1"]<0., "x_norm"] += 2.*(.5-df["x_norm"])

    # isolation of the follower only
    subdf=df[df["follower"]]
    subdf = subdf.groupby("cell", group_keys=False).apply(lambda x: filterConsecutiveSequences(x, min_length=min_length))

    if report:
        print("\n-- report on v.k")
        print("threshold (deg):", opening_angle, round(threshold, 3))
        print(df.head(), "\n")
        print(subdf.head(), "\n")

        #sns.scatterplot(data=subdf, x="x", y="y", hue="cell", palette="deep", edgecolor="black", legend=False)
        #plt.show()
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

    if report:
        print("\n-- report on height")
        print("Amplitude (µm):", round(A, 3))
        print(df[["z", "z_hat", "xbis"]].describe().T, "\n")
    #

    return df
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

def computeForce(df, report=False):
    pos0 = df.columns.get_loc("stoch_norm")

    pd.set_option('display.precision', 15)
    df.insert(pos0+1, "f_pass_norm/stoch_norm", df["f_pass_norm"] / df["stoch_norm"])
    df.insert(pos0+2, "pers_norm/f_pass_norm", df["pers_norm"] / df["f_pass_norm"])
    df.insert(pos0+2, "pers_norm/stoch_norm", df["pers_norm"] / df["stoch_norm"])
    df.insert(pos0+3, "f_passY/stochY", df["f_passY"] / df["stochY"])

    if report:
        print("\n-- report on force")
        print(df[["f_pass_norm", "pers_norm", "stoch_norm", "f_pass_norm/stoch_norm", "pers_norm/f_pass_norm", "pers_norm/stoch_norm"]].describe().T, "\n")
        print(df[["f_pass_norm", "pers_norm", "stoch_norm", "f_pass_norm/stoch_norm", "pers_norm/f_pass_norm", "pers_norm/stoch_norm"]].head(), "\n")
    #

    return df


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

    ## rapport entre les forces
    df=computeForce(df, report=False)

    ## sequences of followers
    df, df_followers=isolateFollowers(df,
                                      opening_angle=90.,
                                      min_length=10,
                                      report=False)

    
    
    return df, df_followers
#

def cles_safe(x, y):
    x = np.asarray(x)
    y = np.asarray(y)

    # Tri pour approche efficace
    x = np.sort(x)
    y = np.sort(y)

    # Comptage efficace via recherche binaire
    idx = np.searchsorted(y, x, side='right')
    cles=idx.mean() / len(y)

    # CLES is symmetric. It is chosen to return the value <=0.5
    if cles<=.5:
        return cles
    else:
        return 1.-cles
#

def performStat(metadf, df, df_f):
    subdf_follo=df_f.copy()
    subdf_NOfollo=df[df["follower"]==False]
    stat_stored=[]

    ### 
    print("=========== Statistical test between Follower and NON Follower populations ===========")

    print("--- norme de F_pass ---")
    stat=stats.mannwhitneyu(subdf_follo["f_pass_norm"], subdf_NOfollo["f_pass_norm"])
    cles=cles_safe(subdf_follo["f_pass_norm"], subdf_NOfollo["f_pass_norm"])
    stat_stored.append([stat.pvalue, cles])
    print("Mann-Whitney U test:", stat)
    print("CLM effect size =", cles)

    print("\n--- norme de persistance ---")
    stat=stats.mannwhitneyu(subdf_follo["pers_norm"], subdf_NOfollo["pers_norm"])
    cles=cles_safe(subdf_follo["pers_norm"], subdf_NOfollo["pers_norm"])
    stat_stored.append([stat.pvalue, cles])
    print("Mann-Whitney U test:", stat)
    print("CLM effect size =", cles)

    print("\n--- norme de stoch ---")
    stat=stats.mannwhitneyu(subdf_follo["stoch_norm"], subdf_NOfollo["stoch_norm"])
    cles=cles_safe(subdf_follo["stoch_norm"], subdf_NOfollo["stoch_norm"])
    stat_stored.append([stat.pvalue, cles])
    print("Mann-Whitney U test:", stat)
    print("CLM effect size =", cles)

    print("\n--- rapport f_pass_norm / stoch_norm ---")
    stat=stats.mannwhitneyu(subdf_follo["f_pass_norm/stoch_norm"], subdf_NOfollo["f_pass_norm/stoch_norm"])
    cles=cles_safe(subdf_follo["f_pass_norm/stoch_norm"], subdf_NOfollo["f_pass_norm/stoch_norm"])
    stat_stored.append([stat.pvalue, cles])
    print("Mann-Whitney U test:", stat)
    print("CLM effect size =", cles)

    return stat_stored
#

def plotVsTime(metadata, df, df_f, axs, N=10, rand=False):
 
    # twin axes def
    ax2=axs[0].twinx()
    ax4=axs[1].twinx()
    ax6=axs[2].twinx()

    axs[0].set_zorder(2)
    ax2.set_zorder(1)
    axs[0].patch.set_visible(False)
    axs[1].set_zorder(2)
    ax4.set_zorder(1)
    axs[1].patch.set_visible(False)
    axs[2].set_zorder(2)
    ax6.set_zorder(1)
    axs[2].patch.set_visible(False)

    # traj to plot
    selected_id=[285]#[285, 24, 123]
    if rand:
        N=min(N, len(df["cell"].unique()))
        selected_id = random.sample(list(df["cell"].unique()), N)
        print("plotted cells:", selected_id)

    subdf=df[df["cell"].isin(selected_id)]
    subdf_f=df_f[df_f["cell"].isin(selected_id)]
    

    # convert seconds to minutes
    T=metadata["Simulation Duration"]/60.
    pos0 = subdf.columns.get_loc("time")
    subdf.insert(pos0+1, "time_min", subdf["time"] / 60.)
    pos0 = subdf_f.columns.get_loc("time")
    subdf_f.insert(pos0+1, "time_min", subdf_f["time"] / 60.)

    ## separation of sequences
    subdf_f = subdf_f.sort_values(["cell", "frame"])
    subdf_f["t_diff"] = subdf_f.groupby("cell")["frame"].diff()
    subdf_f["is_gap"] = (subdf_f["t_diff"] != 1) | subdf_f["t_diff"].isna()
    subdf_f["sequence_id"] = subdf_f.groupby("cell")["is_gap"].cumsum()
 

    # plot cell traj
    for ax in axs:
       # wave-following segments
       for _, g in subdf_f.groupby(["cell", "sequence_id"], sort=False):
           sns.lineplot(data=g,
                   y="x_norm",
                   x="time_min",
                   color=follow,
                   linewidth=3,
                   markers=True,
                   marker="o",
                   markersize=6,
                   markeredgecolor="none",
                   legend=False,
                   zorder=10,
                   ax=ax)
       
       # plot traj for entire simu
       sns.lineplot(data=subdf,
                   y="x_norm",
                   x="time_min",
                   color=no_follow,
                   linewidth=2,
                   markers=True,
                   marker="o",
                   markersize=6,
                   markeredgecolor="none",
                   legend=False,
                   zorder=5,
                   ax=ax)
       
       # axes specs
       ax.set_xticks([0., .25*T, .5*T, .75*T, T])
       ax.set_xlim(left=0., right=T)
       ax.set_ylabel(ylabel=r"$x/\lambda$", fontsize=20)
       ax.set_xlabel(xlabel="Time (min)", fontsize=20)
    #

    # norm of f_pass
    sns.lineplot(data=subdf,
                y="f_pass_norm",
                x="time_min",
                color=force,
                linewidth=3,
                linestyle="--",
                markers=True,
                marker="o",
                markersize=6,
                markeredgecolor="none",
                legend=False,
                ax=ax2)

    ## norm of persistence term
    sns.lineplot(data=subdf,
                y="pers_norm",
                x="time_min",
                color=force,
                linewidth=3,
                linestyle="--",
                markers=True,
                marker="o",
                markersize=6,
                markeredgecolor="none",
                legend=False,
                ax=ax4)

    # rapport des normes
    sns.lineplot(data=subdf,
                y="f_pass_norm/stoch_norm",
                x="time_min",
                color=force,
                linewidth=3,
                linestyle="--",
                markers=True,
                marker="o",
                markersize=6,
                markeredgecolor="none",
                legend=False,
                ax=ax6)

    ax6.set_ylim(bottom=-1., top=10.)
    ax2.set_ylabel(ylabel=r"$\bf || f_{pass} ||$ (µm/s²)", fontsize=20, color=force)
    ax4.set_ylabel(ylabel=r"$\bf || \Pi^{-1} \cdot (v-v_{ext}) ||$ (µm/s²)", fontsize=20, color=force)
    ax6.set_ylabel(ylabel=r"$\bf || f_{pass} || / || \sigma \cdot \eta ||$ (µm/s²)", fontsize=20, color=force)
    ax2.tick_params(axis="y", labelcolor=force)
    ax4.tick_params(axis="y", labelcolor=force)
    ax6.tick_params(axis="y", labelcolor=force)
#

def plotForceDistributions(metadf, df, df_f, stats, axs):
    
    order={axs[0]:("f_pass_norm", r"$\bf || f_{pass} ||$ (µm/s²)", stats[0]),
           axs[1]:("pers_norm", r"$\bf || \Pi^{-1} \cdot (v-v_{ext}) ||$ (µm/s²)", stats[1]),
           axs[2]:("f_pass_norm/stoch_norm", r"$\bf || f_{pass} || / || \sigma \cdot \eta ||$ (µm/s²)", stats[3]) ## on garde stat 3 car on a la stat sur les distributions de stoch
           }
    
    def add_stat(ax, stat):
        symb="NS"
        if stat[0]<0.05:
            symb="*"
        if stat[0]<0.001:
            symb="**"
        
        ax.text(x=.5,
                y=.9, 
                s=symb + "\n" + r"$CLES$" + r"=" + str(round(stat[1],2)),
                transform=ax.transAxes,
                fontsize=13,
                fontweight='bold',
                va='center',
                ha='center')
    #
#

    for ax in axs:
        sns.boxplot(data=df,
                    y=order[ax][0],
                    hue="follower",
                    whis=3.,
                    palette=(no_follow, follow),
                    legend=True,
                    ax=ax
                    )

        ax.set_ylabel(ylabel=order[ax][1], fontsize=20)
        ax.legend(title="Following ?", fontsize= 10, loc=2)
        if order[ax][0]=="f_pass_norm/stoch_norm":
            ax.set_ylim(bottom=-.5, top=10.)
        
        add_stat(ax, order[ax][2])
#

def plotHeatmap(metadf, df, df_f, axs):

    # plot heatmap VS time
    data_to_plot=df[df["time"] % 8 == 0].copy()  ## 8 car toutes les deux minutes
    data_to_plot["time"]=(data_to_plot["time"]/60).astype(int)

    sns.violinplot(data_to_plot,
                   x="time",
                   y="x_norm",
                   hue="follower",
                   density_norm="width",
                   inner="quart",
                   split=True,
                   gap=.1,
                   palette=(no_follow, follow),
                   cut=0,
                   legend=False,
                   ax=axs[0])

    axs[0].set_ylabel(ylabel=r"$x/\lambda$", fontsize=20)
    axs[0].set_xlabel(xlabel="Time (min)", fontsize=20)

    #plot localisation global
    sns.histplot(data=data_to_plot[data_to_plot["follower"]==True],
                 x="x_norm",
                 stat="density",
                 color=follow,
                 ax=axs[1])

    axs[1].set_xlabel(xlabel=r"$x/\lambda$", fontsize=20)
    axs[1].set_ylabel(ylabel="Density", fontsize=20)

    sns.histplot(data=data_to_plot[data_to_plot["follower"]==False],
                 x="x_norm",
                 stat="density",
                 color=no_follow,
                 ax=axs[2])

    axs[2].set_xlabel(xlabel=r"$x/\lambda$", fontsize=20)
    axs[2].set_ylabel(ylabel="Density", fontsize=20)





# ================================================================================================================
# ================================================================================================================
#                                                    UI
# ================================================================================================================
# ================================================================================================================

# study
substrate_regime="dyna"  # dyna or static
assumption="NDSC"

# assets
root=os.getcwd()
NDSC_assets=os.path.join(root, substrate_regime, "assets", assumption)




# ================================================================================================================
# ================================================================================================================
#                                                    Open figure
# ================================================================================================================
# ================================================================================================================
NbLignes, NbColonnes=4, 4
Largeur, Hauteur = 20, 20
Largeur_ins, Hauteur_ins = 1.5, 1.5
fig1 = plt.figure(num="curvature_traj", figsize=(Largeur, Hauteur), constrained_layout=True)
gs1 = GridSpec(NbLignes, NbColonnes, figure=fig1)  #wspace=0.5, hspace=0.3

## font size
plt.rcParams.update({"font.size": 15})

ax1 = fig1.add_subplot(gs1[1, 1:3])
ax3 = fig1.add_subplot(gs1[2, 1:3])
ax5 = fig1.add_subplot(gs1[3, 1:3])
ax7 = fig1.add_subplot(gs1[0, 2:])

ax9 = fig1.add_subplot(gs1[3, 0])
ax10 = fig1.add_subplot(gs1[1, 2])
ax11 = fig1.add_subplot(gs1[2, 2])
ax12 = fig1.add_subplot(gs1[0, 1])

ax13 = fig1.add_subplot(gs1[3, 0])
ax14 = fig1.add_subplot(gs1[1, 0])
ax15 = fig1.add_subplot(gs1[2, 0])
ax16 = fig1.add_subplot(gs1[3, 0])

ax12_ins=inset_axes(ax12, width=Largeur_ins, height=Hauteur_ins, loc="upper right")


#

## add labels
def add_label(ax, txt, x=-0.05, y=1.15):
    ax.text(x, y, txt, transform=ax.transAxes, fontsize=20, fontweight="bold", va="top", ha="right")
#

add_label(ax1, "a")
add_label(ax3, "b")
add_label(ax5, "c")

add_label(ax9, "d")
add_label(ax10, "e")
add_label(ax11, "f")
add_label(ax12, "g")

add_label(ax13, "h")
add_label(ax14, "i")
add_label(ax15, "j")
add_label(ax16, "k")


show=False
save=True
path_to_save=os.path.join(os.getcwd(), substrate_regime, "curv_analysis.svg")

# ================================================================================================================
# ================================================================================================================
#                                                    MAIN
# ================================================================================================================
# ================================================================================================================

print("--------- open data")
path=os.path.join(NDSC_assets, "L=20", "c=0.4","data.txt")

# open file
metadata, data=loadTrajectoryData(path)

# format_data
data, data_followers=prepareData(metadata, data)


plotCurvatureVsXaxis(metadata, ax=ax13)
plotVsTime(metadata, data, data_followers, axs=[ax1, ax3, ax5], rand=False, N=1)
stats=performStat(metadata, data, data_followers)
plotForceDistributions(metadata, data, data_followers, stats, axs=[ax9, ax10, ax11])

plotHeatmap(metadata, data, data_followers, axs=[ax7, ax12, ax12_ins])




print("\n\n----------------------------------")
if show:
    print("... figure is showed")
    plt.show()

if save:
    print("... exporting figure")
    plt.savefig(path_to_save)
    print("... figure is exported")


print("\nPlot Figure 4 END OF PROCESS5\n---------------------------------------\n")


