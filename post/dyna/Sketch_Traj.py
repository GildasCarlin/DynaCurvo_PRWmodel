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

# images
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

# ================================================================================================================
# ================================================================================================================
#                                                    intermediate functions
# ================================================================================================================
# ================================================================================================================

#palettes
pTraj="deep"

cT_ini="#3F2773"  # Time t_ini, t=10 min, t=20 min
cT10="#1F9E89"
cT20="#F2D027"

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



    hyp="DSC(-)"
    if assumption: hyp="DSC(+)"
    
    #ax.text(.5, 1.05, r"{} cells".format(hyp, c), 
    #    fontsize=20, va='center', ha='center', transform=ax.transAxes)
    
    
    carrow="red"
    ax.arrow(0.92, .1, 0., 0.4, width=.02, head_length=.07, color=carrow, alpha=1., transform=ax.transAxes, length_includes_head=True)
    ax.text(0.97, 0.27, "0 h -> 1 h",
            color=carrow, alpha=1., rotation = 'vertical', fontsize=12, fontweight='bold', va='center', ha='center', transform=ax.transAxes)
    ax.arrow(0.7, .8, -0.4, 0., width=.02, head_length=.07, color=carrow, alpha=1., transform=ax.transAxes, length_includes_head=True)
    ax.text(0.51, 0.84, "1 h -> 2 h",
            color=carrow, alpha=1., rotation = 'horizontal', fontsize=12, fontweight='bold', va='center', ha='center', transform=ax.transAxes)

#

def plotBaryc(metadf, df, list_frame, ax, palette):

    # select axe to plot
    c=metadf["Wave Celerity"]
    dt = metadf["Time Interval"]

    # select
    subdf=df.copy().drop(labels=["z", "v_dot_e1", "meanCurvature"], axis=1)
    subdf=subdf[subdf["frame"].isin(list_frame)]

    baryc=subdf.groupby(["frame"]).mean().drop(labels=["cell"], axis=1)

    sns.scatterplot(data=baryc,
                    x='x', 
                    y='y',
                    hue='frame',
                    s=250,
                    ax=ax, 
                    palette=palette,
                    linewidth=2.,   
                    legend=False)
    #
    
    ax.autoscale()

    xticks=[-3.5, -3., -2.5, -2., -1.5, -1., -.5, 0., .5, 1., 1.5, 2., 2]
    yticks=[-.5, 0., .5, 1., 1.5, 2., 2.5, 3., 3.5]
    ax.set_yticks(yticks[::2]) 
    ax.set_xticks(xticks[::2])
    ax.set_xlim(xticks[0], xticks[-1])
    ax.set_ylim(yticks[0], yticks[-1])
    
    norm = mpl.colors.Normalize(vmin=baryc["time"].min()/60., vmax=baryc["time"].max()/60.)
    sm = mpl.cm.ScalarMappable(norm=norm, cmap="viridis")

    cax=ax.inset_axes([0.88, 0.50, 0.025, 0.45])

    cbar = plt.colorbar(sm, cax=cax, orientation="vertical")
    cbar.set_label(r"Time (min)", fontsize=15)
    cbar.ax.tick_params(labelsize=9)
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
NDSC_assets=os.path.join(root, substrate_regime, "assets", "NDSC_rotating")


# figure
# ---------------------------
fig, ax = plt.subplots(figsize=(10, 7))

## font size
plt.rcParams.update({"font.size": 15})

## add labels
def add_label(ax, txt, x=-0.05, y=1.15):
    ax.text(x, y, txt, transform=ax.transAxes, fontsize=20, fontweight="bold", va="top", ha="right")
#

add_label(ax, 'f')

show=False
save=True
path_to_save=os.path.join(root, substrate_regime, "Exotic_surfaces.svg")



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
metadata_NDSC, data_NDSC=loadTrajectoryData(os.path.join(NDSC_assets, "data.txt"))

print('---------------------- format traj')
data_NDSC=formatTraj(data_NDSC, center=True, convert=True)

print("---------------------- plot traj")
plotTraj(metadf=metadata_NDSC, 
          df=    data_NDSC, 
          ax=      ax,
          palette=pTraj,
          N=       50
          )

plotBaryc(metadf=metadata_NDSC,
          df=data_NDSC,
          list_frame=[10, 80, 160, 240, 320, 400, 480],
          ax=ax,
          palette="viridis")




print("\n\n----------------------------------")
if show:
    print("... figure is showed")
    plt.show()

if save:
    print("... exporting figure")
    plt.savefig(path_to_save)
    print("... figure is exported")


print("\nPlot Figure 5 END OF PROCESS5\n---------------------------------------\n")