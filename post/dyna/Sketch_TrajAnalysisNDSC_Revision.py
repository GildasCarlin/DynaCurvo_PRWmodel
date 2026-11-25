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
import scikit_posthocs as sph
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

def plotDist(data, quantity, stats, ax, whis = 1.5, perform_stat=True, plot_legend=True, color=None):

    # plot distributions
    #-----------------------------------------------------
    sns.boxplot(data=data, 
                   x='celerity', 
                   y=quantity, 
                   linewidth=1, 
                   color=color, 
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
    YLabel={"dx":r'$D_{\mathbf{k}}$'+" at "+r"$T$=20 min",
        "speedObs":r'Cell speed $|| \mathbf{v^{obs}} ||_2$ (µm/s)',
        "speedCell":r'Cell speed $|| \mathbf{v} ||_2$ (µm/s)',
        "VzSub":'Surface velocity (µm/s)',
        "r_d":r'$f$'+" at "+r"$T$=20 min",
        "last_y":r"$B_T$"+" at "+r'$T$=20 min (mm)'
        }
    ax.set_ylabel(YLabel[quantity], fontsize=25)
    ax.set_xlabel(r'Wave celerity (µm/s)', fontsize=25)
    ax.tick_params(axis='both', which='major', labelsize=20) 
    if plot_legend:
        ax.legend(fontsize=20, title='Assumption', title_fontsize=22, loc = "best")
    #

    def add_stat(ax, stats):
        celerities = stats.columns
        X=ax.get_xticks()
        l=0.2
        count=0
        for c1, c2 in pairwise(celerities):
            x1, x2=X[count], X[count+1]
            x_convert=.5*(x1+x2)
            ymin, ymax = ax.get_ylim()
            xmin, xmax = ax.get_xlim()
            y_convert = ymin + (.95-0.04) * (ymax - ymin)
            ax.hlines(y=y_convert, xmin=x_convert-l, xmax=x_convert+l, color='black', linewidth=2.)
            count+=1

            p, d = stats.loc[c1, c2]
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
                


    if perform_stat: add_stat(ax, stats)
#



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

def performStats(df, value_col, cat_col="celerity"):
    
    groups = {cat: g[value_col].to_numpy() for cat, g in df.groupby(cat_col)}
    
    celerities = list(groups)
    res = pd.DataFrame(index=celerities, columns=celerities, dtype=object)
    
    # post hoc dunn's test. Kruskall wallis already performed in Fig2_statistical_analysis.py
    p_values = sph.posthoc_dunn([groups[c] for c in celerities], p_adjust="bonferroni")
    p_values.index = celerities
    p_values.columns = celerities

    for c1, c2 in combinations(celerities, 2):

        # mann Whitney U test and CLES effect size
        #u, p=stats.mannwhitneyu(X1, X2)
        p = p_values.loc[c1, c2]
        cles=computeCLES(groups[c1], groups[c2])

        res.loc[c1, c2]=(p, cles)
        res.loc[c2, c1]=(p, cles)

    for c in celerities:
        res.loc[c, c]=(np.nan, np.nan, np.nan)

    return res

def plotTraj(metadata, data, ax, N=60, centered=True):

    c=metadata["Wave Celerity"]
    
    ## random cells id to plot
    IDlist = np.random.choice(data['cell'].unique(), size=N, replace=False)
    
    # extract trajectories to plot
    data_to_plot = data[data['cell'].isin(IDlist)].copy()
 
    sns.scatterplot(data=data_to_plot, 
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

    xticks=[-1., -0.5, 0., .5, 1.]
    yticks=[-0.25, 0., .25, .5, .75]
    ax.set_yticks(yticks) 
    ax.set_xticks(xticks)
    ax.set_xlim(xticks[0], xticks[-1])
    ax.set_ylim(yticks[0], yticks[-1])

    ax.text(.5, 1.05, r'DSC(-) cells, $c$={} µm/s'.format(c), 
        fontsize=20, va='center', ha='center', transform=ax.transAxes)
    
    if c == 0.0:
        ax.arrow(0.1, .3, 0., 0.4, width=.03, head_length=.07, color='red', transform=ax.transAxes, length_includes_head=True)
        ax.text(0.15, 0.5, "Wave",
                color='red', rotation = 'vertical', fontsize=12, fontweight='bold', va='center', ha='center', transform=ax.transAxes)
        ax.arrow(.2, 0.9, -0.15, 0., width=.03, head_length=.07, color='red', length_includes_head=True, transform=ax.transAxes)
        ax.arrow(.2, 0.9, 0.15, 0., width=.03, head_length=.07, color='red', length_includes_head=True, transform=ax.transAxes)
        ax.text(0.2, 0.81, "Grooves",
                color='red', rotation = 'horizontal', fontsize=12, fontweight='bold', va='center', ha='center', transform=ax.transAxes)


    if centered:
        ax.scatter([0.0], [0.0], s=200, color='deeppink')
    else:
        rect=mpatches.Rectangle([-0.2, -0.2], width=0.4, height=0.4, ec='deeppink', ls="--", alpha=1., fill=0, lw=4, zorder=2)
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

# cases of study
WL=[20.]
CL=[0., .1, .4, .6, 2.]



# ================================================================================================================
# ================================================================================================================
#                                                    Open figure
# ================================================================================================================
# ================================================================================================================
NbLignes, NbColonnes=2, 5
Largeur, Hauteur = 25, 10
Largeur_ins, Hauteur_ins = 1.7, 1.7
fig1 = plt.figure(num="TrajAnalysis_NDSC", figsize=(Largeur, Hauteur), constrained_layout=True)
gs1 = GridSpec(NbLignes, NbColonnes, figure=fig1) 

## font size
plt.rcParams.update({"font.size": 15})

## axes pour plot les trajectoires / barycentres en insert
ax_traj0 = fig1.add_subplot(gs1[0, 0])
ax_traj1 = fig1.add_subplot(gs1[0, 1])
ax_traj2 = fig1.add_subplot(gs1[0, 2])
ax_traj3 = fig1.add_subplot(gs1[0, 3])
ax_traj4 = fig1.add_subplot(gs1[0, 4])

## axes d'insert de barycentres
ax_ins0 = inset_axes(ax_traj0, width=Largeur_ins, height=Hauteur_ins)
ax_ins1 = inset_axes(ax_traj1, width=Largeur_ins, height=Hauteur_ins)
ax_ins2 = inset_axes(ax_traj2, width=Largeur_ins, height=Hauteur_ins)
ax_ins3 = inset_axes(ax_traj3, width=Largeur_ins, height=Hauteur_ins)
ax_ins4 = inset_axes(ax_traj4, width=Largeur_ins, height=Hauteur_ins)

## axes pour plot les autres indicateurs (Dx, velocity etc)
ax1 = fig1.add_subplot(gs1[1, 0])
ax2 = fig1.add_subplot(gs1[1, 1])
ax3 = fig1.add_subplot(gs1[1, 2])
ax4 = fig1.add_subplot(gs1[1, 3])
ax5 = fig1.add_subplot(gs1[1, 4])


## add labels
def add_label(ax, txt, x=-0.05, y=1.15):
    ax.text(x, y, txt, transform=ax.transAxes, fontsize=20, fontweight="bold", va="top", ha="right")
#

add_label(ax_traj0, 'a.1')
add_label(ax_traj1, 'a.2')
add_label(ax_traj2, 'a.3')
add_label(ax_traj3, 'a.4')
add_label(ax_traj4, 'a.5')

add_label(ax1, 'b')
add_label(ax2, 'c')
add_label(ax3, 'd')
add_label(ax4, 'e')
add_label(ax5, 'f')


show=False
save=True
path_to_save=os.path.join(root, substrate_regime, "trajAnalysis_NDSC.svg")

# ================================================================================================================
# ================================================================================================================
#                                                    MAIN function
# ================================================================================================================
# ================================================================================================================

print('\n------------------------------------ get general indicators')
print('... read data')

# open data
assumption = "DSC(-)"
last = pd.read_csv(os.path.join(NDSC_assets, "last.csv"))
last = last[last['celerity'].isin(CL) & last['wavelength'].isin(WL)]
last.insert(0, 'assumption', assumption)

dx = pd.read_csv(os.path.join(NDSC_assets, "dx.csv"))
dx = dx[dx['celerity'].isin(CL) & dx['wavelength'].isin(WL)]
dx.insert(0, 'assumption', assumption)

vd = pd.read_csv(os.path.join(NDSC_speed_assets, "vd.csv"))
vd = vd[vd['celerity'].isin(CL) & vd['wavelength'].isin(WL)]
vd.insert(0, 'assumption', assumption)

### convert µm to mm (barycenter position)
last['last_y'] /= 1000.


print('... plot---last')
plotDist(data=       last,
         quantity=   "last_y",
         stats=      performStats(last, "last_y"),
         ax=         ax1,
         color=      cNDSC,
         plot_legend=False)


subdata_last=last[last['celerity'].isin(CL[1:])] # rd is not defined for c=0µm/s so it is removed
print('... plot r_d')
plotDist(data=       subdata_last,
         quantity=   'r_d',
         stats=      performStats(subdata_last, "r_d"),
         ax=         ax2,
         color=      cNDSC,
         plot_legend=False)

print('---------------------- plot dx')
plotDist(data=       dx,
         quantity=   'dx',
         stats=      performStats(dx, "dx"),
         ax=         ax3,
         color=      cNDSC,
         plot_legend=False)

print('---------------------- plot substrate vertical velocity')
plotDist(data=       vd.copy().sample(frac=0.03), 
          quantity=   "VzSub", 
          perform_stat=False,
          stats=None,
          ax=         ax4, 
          color=      cCellSpeed,
          plot_legend=False
          )

print('---------------------- plot cell speed')
vd_to_plot=vd.copy().sample(frac=0.03)[["celerity", "speedObs", "speedCell"]].set_index("celerity").stack().reset_index()
vd_to_plot.columns = ["celerity", "variable", "speed"]

sns.boxplot(data=vd_to_plot, 
            x='celerity', 
            y="speed",
            hue="variable",
            palette=(cSubSpeed, cNDSC),
            linewidth=1, 
            ax=ax5,
            whis = 3.,
            )
#

ax5.set_ylabel("Cell speed (µm/s)", fontsize=25)
ax5.set_xlabel(r'Wave celerity (µm/s)', fontsize=25)
ax5.tick_params(axis='both', which='major', labelsize=20) 
ax5.legend()
#legend
items = []
items.append(mpatches.Patch(facecolor=cCellSpeed, label=r"$||\mathbf{v^{obs}}||_2$", edgecolor='black'))
items.append(mpatches.Patch(facecolor=cNDSC, label=r"$||\mathbf{v}||_2$", edgecolor='black'))
ax5.legend(fontsize=20, loc=2, handles = items)

mapping_CL_to_axes = {CL[0]:(ax_traj0, ax_ins0),
                      CL[1]:(ax_traj1, ax_ins1),
                      CL[2]:(ax_traj2, ax_ins2), 
                      CL[3]:(ax_traj3, ax_ins3),
                      CL[4]:(ax_traj4, ax_ins4)}

for celerity in CL:

    print('\n---------------------- get traj for c =',str(celerity))
    metadata, data=loadTrajectoryData(os.path.join(NDSC_assets, "L=20","c="+str(celerity), "data.txt"))

    print('---------------------- format traj')
    data=formatTraj(data, center=True, convert=True)

    print("---------------------- plot traj")
    axs=mapping_CL_to_axes[celerity]
    plotTraj(metadata=metadata, 
              data=    data, 
              ax=      axs[0],
              N=       50
              )

    print("---------------------- plot barycenter")
    plotBaryc(metadf=  metadata, 
               df=      data, 
               list_frame=[10, 40, 80],
               ax=        axs[1],
               palette=   [cT_ini, cT10, cT20]
               )



print("\n\n----------------------------------")
if show:
    print("... figure is showed")
    plt.show()

if save:
    print("... exporting figure")
    plt.savefig(path_to_save)
    print("... figure is exported")


print("\nPlot Figure 2 END OF PROCESS5\n---------------------------------------\n")