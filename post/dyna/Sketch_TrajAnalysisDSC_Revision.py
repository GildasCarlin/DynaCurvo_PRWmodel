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
from itertools import combinations, pairwise

# images
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.gridspec import GridSpec
from matplotlib.collections import LineCollection
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.patches as mpatches


# ================================================================================================================
# ================================================================================================================
#                                                    intermediate functions
# ================================================================================================================
# ================================================================================================================

#palettes
cT_ini="#3F2773"  # Time t_ini, t=10 min, t=20 min
cT10="#1F9E89"
cT20="#F2D027"

cNDSC="#E1BE6A"   # assumption
cDSC="#40B0A6"
cSubSpeed="darkgray"
cCellSpeed="#E2DED4"

cTdynaLow="#80FF5B"
cTdynaHigh="#00E4E6"


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

def formatTraj(df, center=False, convert=False):

    # theta column
    df=df.drop("theta", axis=1)
    
    # trajectory are centered at origin.
    if center: df[['x', 'y']] -= df.groupby('cell')[['x', 'y']].transform('first')

    # js frame vs python frame: the y-axis is reversed
    df['y'] *= -1.0

    # convert µm to mm
    if convert: df[['x', 'y']] /= 1000.0

    return df
#

def plotDist(data, quantity, stats, ax, whis = 1.5, plot_legend=True, palette=None):

    # plot distributions
    #-----------------------------------------------------
    sns.boxplot(data=data, 
                   x='celerity', 
                   y=quantity,
                   hue="assumption", 
                   linewidth=1, 
                   palette=palette, 
                   ax=ax,
                   whis = whis,
                   legend=plot_legend,
                   zorder=10
                   )
    #

    # axis specificity
    if quantity=='dx':
        ax.axhline(y=.5, linewidth = 2, color='#8C0335', linestyle='--', zorder=1)
        ax.set_ylim(bottom=0., top=1.07)
    
    if quantity in ['VzSub']:
        celerities=sorted(data["celerity"].unique())
        X, Y = [], []
        for i, cel in enumerate(celerities):
            Y.append(10.*2.*np.pi/20. * cel)
            X.append(i)
        ax.scatter(X, Y, marker="x", s=400, linewidth=2, color="red", zorder=10)
        ax.scatter(X, [-1.*v for v in Y], marker="x", s=400, linewidth=2, color="red", zorder=1)
        
    if quantity in ["speedObs"]:
        ax.set_ylim(bottom=-.1, top=9.)
    
    if quantity in ["speedCell"]:
        ax.set_ylim(bottom=-.1, top=6.)

    if quantity=='r_d':
        ax.axhline(y=0., linewidth = 2, color='#8C0335', linestyle='--', zorder=1)
        ax.axhline(y=1., linewidth = 2, color='#8C0335', linestyle='--', zorder=1)
        ax.set_ylim(bottom=-2., top=2.)
    
    if quantity=='last_y':
        ax.axhline(y=0., linewidth = 2, color='#8C0335', linestyle='--', zorder=1)
        ax.set_ylim(bottom=-.2, top=0.85)
    #

    # labels fig
    YLabel={"dx":r'$D_{\mathbf{k}}$'+" at "+r'$T$=20 min',
        "speedObs":r'Cell speed $|| \mathbf{v^{obs}} ||_2$ (µm/s)',
        "speedCell":r'Cell speed $|| \mathbf{v} ||_2$ (µm/s)',
        "VzSub":'Surface velocity (µm/s)',
        "r_d":r'$f$'+" at "+r'$T$=20 min',
        "last_y":r"$B_T$"+" at "+r'$T$=20 min (mm)'
        }
    ax.set_ylabel(YLabel[quantity], fontsize=25)
    ax.set_xlabel(r'Wave celerity (µm/s)', fontsize=25)
    ax.tick_params(axis='both', which='major', labelsize=20) 
    if plot_legend:
        ax.legend(fontsize=15, title='Assumption', title_fontsize=17, loc = "best")
    #

    def add_stat(ax, stats):
        celerities = stats.keys()
        X=ax.get_xticks()
        l=0.2
        count=0
        for c in celerities:
            x_convert=X[count]
            ymin, ymax = ax.get_ylim()
            xmin, xmax = ax.get_xlim()
            y_convert = ymin + (.95-0.04) * (ymax - ymin)
            ax.hlines(y=y_convert, xmin=x_convert-l, xmax=x_convert+l, color='black', linewidth=2.)
            count+=1

            _, p, d = stats[c]
            symb="NS"
            if p<0.05:  symb="*"
            if p<0.001: symb="**"
            xc=(x_convert-xmin)/(xmax-xmin)
            ax.text(x=xc, y=.95, 
                    s=symb + "\n" + r"$CLES$" + r"=" + str(round(d,2)),
                    fontsize=8,
                    transform=ax.transAxes,
                    va='center',
                    ha='center')
                


    add_stat(ax, stats)
#

def computeCLES(x, y):
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

def performStats(df, value_col, col1="celerity", col2="assumption"):
    
    b1, b2 = sorted(df[col2].dropna().unique())
    res = {}

    for c, g in df.groupby(col1):
        X1 = g.loc[g[col2] == b1, value_col].to_numpy()
        X2 = g.loc[g[col2] == b2, value_col].to_numpy()

        u, p=stats.mannwhitneyu(X1, X2)
        cles=computeCLES(X1, X2)

        res[c]=(u, p, cles)

    return res


def plotTraj(metadf, df, ax, N=60, centered=True, palette="deep"):

    c=metadf["Wave Celerity"]
    assumption=metadf["surf_assumption"]
    
    ## random cells id to plot
    IDlist = np.random.choice(df['cell'].unique(), size=N, replace=False)
    
    # extract trajectories to plot
    data_to_plot = df[df['cell'].isin(IDlist)].copy()
 
    sns.scatterplot(data=data_to_plot, 
                    x='x', 
                    y='y', 
                    hue='cell', 
                    ax=ax, 
                    palette=palette, 
                    marker='o', 
                    edgecolor='none'
                    )
    
    ax.get_legend().set_visible(False)
    ax.set_xlabel('x (mm)', fontsize=22)
    ax.set_ylabel('y (mm)', fontsize=22)
    ax.tick_params(axis='both', which='major', labelsize=15)

    xticks=[-1., -0.5, 0., .5, 1.]
    yticks=[-0.25, 0., .25, .5, .75]
    ax.set_yticks(yticks) 
    ax.set_xticks(xticks)
    ax.set_xlim(xticks[0], xticks[-1])
    ax.set_ylim(yticks[0], yticks[-1])

    hyp="DSC(-)"
    if assumption: hyp="DSC(+)"
    
    ax.text(.5, 1.05, r"{} cells, $c$={} µm/s".format(hyp, c), 
        fontsize=20, va='center', ha='center', transform=ax.transAxes)
    
    
    carrow, cdot="red", "deeppink"
    if not assumption: carrow, cdot="gray", "black"
    if c == 0.1:
        ax.arrow(0.1, .3, 0., 0.4, width=.03, head_length=.07, color=carrow, transform=ax.transAxes, length_includes_head=True)
        ax.text(0.15, 0.5, "Wave",
                color=carrow, rotation = 'vertical', fontsize=12, fontweight='bold', va='center', ha='center', transform=ax.transAxes)
        ax.arrow(.2, 0.9, -0.15, 0., width=.03, head_length=.07, color=carrow, length_includes_head=True, transform=ax.transAxes)
        ax.arrow(.2, 0.9, 0.15, 0., width=.03, head_length=.07, color=carrow, length_includes_head=True, transform=ax.transAxes)
        ax.text(0.2, 0.81, "Grooves",
                color=carrow, rotation = 'horizontal', fontsize=12, fontweight='bold', va='center', ha='center', transform=ax.transAxes)

    if centered:
        ax.scatter([0.0], [0.0], s=200, color=cdot)
    else:
        rect=mpatches.Rectangle([-0.2, -0.2], width=0.4, height=0.4, ec=cdot, ls="--", alpha=1., fill=0, lw=4, zorder=2)
        ax.add_patch(rect)
#

def plotBaryc(metadf, df, list_frame, ax, palette):

    # select axe to plot
    c=metadf["Wave Celerity"]
    dt = metadf["Time Interval"]

    pop = df[df['frame'].isin(list_frame)]
    sns.kdeplot(data=pop,
                x='x', 
                y='y',
                hue='frame',
                ax=ax, 
                palette=palette,   
                alpha=.8,
                levels=5,
                thresh=.01,
                zorder=1,
                legend=False
                )
    #

    # plot barycenter
    baryc = pop.groupby('frame').mean()
    ax.scatter(baryc['x'], baryc['y'], marker='o', s=150, c=palette, alpha = 1., edgecolors='black', zorder=2)

    if c in [0.0, 0.1, 2.0]:
        items = []
        for i in range(len(list_frame)):
            items.append(mpatches.Patch(facecolor=palette[i], label=float(list_frame[i])*dt/60., edgecolor='black'))
        #
        ax.legend(handles = items, fontsize=11, title='t (min)', title_fontsize=10)
    #

    ax.set_xlabel('')
    ax.set_ylabel('') 

    xticks=[-1., -.5, 0., .5, 1.]
    yticks=[-0.25, 0., .25, .5, .75]
    ax.set_yticks(yticks) 
    ax.set_xticks(xticks)
    ax.set_xlim(xticks[0], xticks[-1])
    ax.set_ylim(yticks[0], yticks[-1])
    ax.tick_params(axis='both', which='major', labelsize=9)
    ax.patch.set_alpha(0.5)
#



# ================================================================================================================
# ================================================================================================================
#                                                    UI
# ================================================================================================================
# ================================================================================================================

# study
substrate_regime="dyna"  # dyna or static

# assets
root=os.getcwd()
NDSC_assets=os.path.join(root, substrate_regime, "assets", "NDSC")
NDSC_speed_assets=os.path.join(root, substrate_regime, "assets", "NDSC_speed")
DSC_assets=os.path.join(root, substrate_regime, "assets", "DSC")
DSC_speed_assets=os.path.join(root, substrate_regime, "assets", "DSC_speed")

# cases of study
WL=[20.]
CL=[0., .1, .4, .6, 2.]



# ================================================================================================================
# ================================================================================================================
#                                                    Open figure
# ================================================================================================================
# ================================================================================================================
NbLignes, NbColonnes=3, 4
Hauteur, Largeur,  = 15, 20
Largeur_ins, Hauteur_ins = 1.7, 1.7
fig1 = plt.figure(num="TrajAnalysis_DSC", figsize=(Largeur, Hauteur), constrained_layout=True)
gs1 = GridSpec(NbLignes, NbColonnes, figure=fig1) 

## font size
plt.rcParams.update({"font.size": 15})

## axes pour plot les trajectoires / barycentres en insert
ax_traj0 = fig1.add_subplot(gs1[0, 0])
ax_traj1 = fig1.add_subplot(gs1[0, 1])
ax_traj2 = fig1.add_subplot(gs1[0, 2])
ax_traj3 = fig1.add_subplot(gs1[0, 3])

ax_traj4 = fig1.add_subplot(gs1[1, 0])
ax_traj5 = fig1.add_subplot(gs1[1, 1])
ax_traj6 = fig1.add_subplot(gs1[1, 2])
ax_traj7 = fig1.add_subplot(gs1[1, 3])

## axes d'insert de barycentres
ax_ins0 = inset_axes(ax_traj4, width=Largeur_ins, height=Hauteur_ins)
ax_ins1 = inset_axes(ax_traj5, width=Largeur_ins, height=Hauteur_ins)
ax_ins2 = inset_axes(ax_traj6, width=Largeur_ins, height=Hauteur_ins)
ax_ins3 = inset_axes(ax_traj7, width=Largeur_ins, height=Hauteur_ins)

## axes pour plot les autres indicateurs (Dx, velocity etc)
ax1 = fig1.add_subplot(gs1[2, 0])
ax2 = fig1.add_subplot(gs1[2, 1])
ax3 = fig1.add_subplot(gs1[2, 2])
ax4 = fig1.add_subplot(gs1[2, 3])


## add labels
def add_label(ax, txt, x=-0.05, y=1.15):
    ax.text(x, y, txt, transform=ax.transAxes, fontsize=20, fontweight="bold", va="top", ha="right")
#

add_label(ax_traj0, 'a.1')
add_label(ax_traj1, 'a.2')
add_label(ax_traj2, 'a.3')
add_label(ax_traj3, 'a.4')

add_label(ax_traj4, 'b.1')
add_label(ax_traj5, 'b.2')
add_label(ax_traj6, 'b.3')
add_label(ax_traj7, 'b.4')

add_label(ax1, 'c')
add_label(ax2, 'd')
add_label(ax3, 'e')
add_label(ax4, 'f')


show=False
save=True
path_to_save=os.path.join(root, substrate_regime, "trajAnalysis_DSC.svg")

# ================================================================================================================
# ================================================================================================================
#                                                    MAIN function
# ================================================================================================================
# ================================================================================================================


# plot distribution ------------------------------------------------------
print('\n------------------------------------ get general indicators')
print('... read data')

# open NDSC data
assumption = "DSC(-)"
last_NDSC = pd.read_csv(os.path.join(NDSC_assets, "last.csv"))
last_NDSC = last_NDSC[last_NDSC['celerity'].isin(CL) & last_NDSC['wavelength'].isin(WL)]
last_NDSC.insert(0, 'assumption', assumption)

dx_NDSC = pd.read_csv(os.path.join(NDSC_assets, "dx.csv"))
dx_NDSC = dx_NDSC[dx_NDSC['celerity'].isin(CL) & dx_NDSC['wavelength'].isin(WL)]
dx_NDSC.insert(0, 'assumption', assumption)

# open DSC data
assumption = "DSC(+)"
last_DSC = pd.read_csv(os.path.join(DSC_assets, "last.csv"))
last_DSC = last_DSC[last_DSC['celerity'].isin(CL) & last_DSC['wavelength'].isin(WL)]
last_DSC.insert(0, 'assumption', assumption)

dx_DSC = pd.read_csv(os.path.join(DSC_assets, "dx.csv"))
dx_DSC = dx_DSC[dx_DSC['celerity'].isin(CL) & dx_DSC['wavelength'].isin(WL)]
dx_DSC.insert(0, 'assumption', assumption)

# concat
last = pd.concat([last_NDSC, last_DSC], ignore_index=True)
dx =   pd.concat([dx_NDSC, dx_DSC], ignore_index=True)

# convert µm to mm (barycenter position)
last['last_y'] /= 1000.


print('... plot---last')
plotDist(data=       last,
         quantity=   "last_y",
         stats=      performStats(last, "last_y"),
         ax=         ax1,
         palette=    (cNDSC, cDSC),
         plot_legend=False)

subdata_last=last[last['celerity'].isin(CL[1:])] # rd is not defined for c=0µm/s so it is removed
print('... plot r_d')
plotDist(data=       subdata_last,
         quantity=   'r_d',
         stats=      performStats(subdata_last, "r_d"),
         ax=         ax2,
         palette=    (cNDSC, cDSC),
         plot_legend=True)

print('---------------------- plot dx')
plotDist(data=       dx,
         quantity=   'dx',
         stats=      performStats(dx, "dx"),
         ax=         ax3,
         palette=    (cNDSC, cDSC),
         plot_legend=True)



# plot trajectories ------------------------------------------------------

mapping_CL_to_axes = {CL[1]:(ax_traj0, ax_traj4, ax_ins0),
                      CL[2]:(ax_traj1, ax_traj5, ax_ins1),
                      CL[3]:(ax_traj2, ax_traj6, ax_ins2), 
                      CL[4]:(ax_traj3, ax_traj7, ax_ins3)}

for celerity in CL[1:]:

    print('\n---------------------- get traj for c =',str(celerity))
    metadata_NDSC, data_NDSC=loadTrajectoryData(os.path.join(NDSC_assets, "L=20","c="+str(celerity), "data.txt"))
    metadata_DSC, data_DSC=loadTrajectoryData(os.path.join(DSC_assets, "L=20","c="+str(celerity), "data.txt"))

    print('---------------------- format traj')
    data_NDSC=formatTraj(data_NDSC, center=True, convert=True)
    data_DSC=formatTraj(data_DSC, center=True, convert=True)

    print("---------------------- plot traj")
    axs=mapping_CL_to_axes[celerity]
    plotTraj(metadf=metadata_NDSC, 
              df=    data_NDSC, 
              ax=      axs[0],
              palette="gray",
              N=       50
              )

    plotTraj(metadf=metadata_DSC, 
              df=    data_DSC, 
              ax=      axs[1],
              N=       50
              )

    print("---------------------- plot barycenter")
    plotBaryc(metadf=  metadata_DSC, 
               df=      data_DSC, 
               list_frame=[10, 40, 80],
               ax=        axs[2],
               palette=   [cT_ini, cT10, cT20]
               )
#



# plot Dk vs time ------------------------------------------------------
# either DSC(+) or DSC(-)
# L=20 µm and c=0.4 µm/s

def getTemporalDx(metadf, df, start=8, step=8):
    
    list_frames=df["time"].unique()[start::step]
    Dx_by_time=[]

    for t in list_frames:
        subdf=df[df["time"]<=t]
        Dx_by_cell = {}
        for i in subdf['cell'].unique(): # calculates Dx for each cell
            Sx, Sy = subdf.loc[subdf['cell']==i,['x','y']].diff().abs().sum(axis=0)
            Dx_by_cell[i]=Sx/(Sx+Sy)
        #
        Dx_by_time.append(Dx_by_cell)

    res=pd.DataFrame(Dx_by_time)
    res=res.stack().reset_index()
    res.columns=["time", "cell", "dx"]
    
    mapping = dict(zip(res["time"].unique(), list_frames))
    res["time"] = res["time"].map(mapping)
    res["time"] /= 60.  # conversion sec min
    return res
# 



# new assets

NDSC_assets=os.path.join(root, substrate_regime, "assets", "NDSC_taudyna")
DSC_assets =os.path.join(root, substrate_regime, "assets", "DSC_taudyna")
taux=[30.0, 120.0, 480.0]

# boucle sur les NDSC

print("\n... Taudyna NDSC")
metadata_NDSC, data_NDSC=loadTrajectoryData(os.path.join(NDSC_assets, "L=20", "c=0.4", "taudyna=120.0", "data.txt"))
res=getTemporalDx(metadata_NDSC, data_NDSC)
sns.lineplot(data=res,
             x="time",
             y="dx",
             err_style="bars",
             errorbar=("sd", 1),
             color=cNDSC,
             linewidth=3,
             ax=ax4)
#


palette={30.0:cTdynaLow,
         120.0:cDSC,
         480.0:cTdynaHigh}
items = []
for tau in taux:
    print("... Taudyna DSC {}".format(tau))

    items.append(mpatches.Patch(facecolor=palette[tau], label="{}".format(round(tau/60.,2)), edgecolor='black'))
    
    metadata_DSC, data_DSC=loadTrajectoryData(os.path.join(DSC_assets, "L=20", "c=0.4", "taudyna={}".format(tau), "data.txt"))
    res=getTemporalDx(metadata_DSC, data_DSC)

    sns.lineplot(data=res,
                 x="time",
                 y="dx",
                 err_style="bars",
                 errorbar=("sd", 1),
                 color=palette[tau],
                 linewidth=3,
                 ax=ax4)
#

items.append(mpatches.Patch(facecolor=cNDSC, label="DSC(-)", edgecolor='black'))
ax4.text(0.3, 0.9, r"$c$=0.4 µm/s",
                color="black", rotation = 'horizontal', fontsize=12, va='center', ha='center', transform=ax4.transAxes)
ax4.axhline(y=.5, linewidth = 2, color='#8C0335', linestyle='--', zorder=1)
ax4.set_ylim(bottom=0.3, top=1.)
ax4.legend(handles = items, fontsize=11, title=r"$\tau_{dyna}$ (min)", title_fontsize=10, loc="upper right")
ax4.set_ylabel(r'$D_{\mathbf{k}}$', fontsize=25)
ax4.set_xlabel(r'Time (min)', fontsize=25)

































print("\n\n----------------------------------")
if show:
    print("... figure is showed")
    plt.show()

if save:
    print("... exporting figure")
    plt.savefig(path_to_save)
    print("... figure is exported")


print("\nPlot Figure 3 END OF PROCESS5\n---------------------------------------\n")