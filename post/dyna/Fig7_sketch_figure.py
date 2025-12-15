# ================================================================================================================
# ================================================================================================================
#                                                    packages importation
# ================================================================================================================
# ================================================================================================================

import json

import numpy as np
#from sklearn.metrics import r2_score

import os
import regex as re
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.lines import Line2D
from matplotlib.gridspec import GridSpec
import seaborn as sns

# ================================================================================================================
# ================================================================================================================
#                                                    open figure
# ================================================================================================================
# ================================================================================================================

######## open figure

NbLignes, NbColonnes=3, 3
Largeur, Hauteur = 25, 25
Largeur_ins, Hauteur_ins = 3.4, 3.4
fig1 = plt.figure(num='role_of_q', figsize=(Largeur, Hauteur), constrained_layout=True)
gs1 = GridSpec(NbLignes, NbColonnes, figure=fig1)

## axes pour scatter les trajectoires
ax1 = fig1.add_subplot(gs1[0, 0])
ax2 = fig1.add_subplot(gs1[0, 1])

## axes pour les statistiques
ax4 = fig1.add_subplot(gs1[1, 0])
ax5 = fig1.add_subplot(gs1[1, 1])
ax3 = fig1.add_subplot(gs1[0, 2])
ax6 = fig1.add_subplot(gs1[1, 2])
ax7 = fig1.add_subplot(gs1[2, 0])
ax8 = fig1.add_subplot(gs1[2, 1])
ax9 = fig1.add_subplot(gs1[2, 2])

#ax_ins1 = inset_axes(ax2, width=Largeur_ins, height=Hauteur_ins)
#ax_ins5 = inset_axes(ax5, width=Largeur_ins, height=Hauteur_ins)

## add labels
def add_label(ax, txt, x=-0.05, y=1.1):
    ax.text(x, y, txt, transform=ax.transAxes, fontsize=20, fontweight='bold', va='top', ha='right')
#

#add_label(ax_traj1, 'A')
#add_label(ax_traj2, 'B')
#add_label(ax1, 'C')
#add_label(ax2, 'D')
#add_label(ax3, 'E')
#add_label(ax4, 'F')
#add_label(ax5, 'G')
#add_label(ax6, 'H')
#add_label(ax7, 'I')

## font size
plt.rcParams.update({'font.size': 50})


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
def get_assumption(metadata):
    c, w=metadata["curvature_assumption"], metadata["surf_assumption"]

    if c==1.0 and w==1.0:
        return 'DSC'
    
    if c==1.0 and w==0.0:
        return 'NDSC'
#


def plot_traj(metadata, data, N, ax, xticks, yticks):

    ## random cells id to plot
    IDlist = np.random.choice(data['cell'].unique(), size=N, replace=False)

    # select axe to plot
    #assumption=get_assumption(metadata)
    #ax = mapping[assumption]
   ## wl, c = metadata['Wavelength'], metadata['Wave Celerity']

    
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
    #

    ax.text(0.12, 0.55, path[len(substrate_regime)+len(study)+len("assets")+3:-1-len("data.txt")], transform=ax.transAxes, fontsize=20)
    
    ax.get_legend().set_visible(False)
    ax.set_xlabel('x (mm)', fontsize=22)
    ax.set_ylabel('y (mm)', fontsize=22)
    ax.tick_params(axis='both', which='major', labelsize=15)

    ax.set_yticks(yticks) 
    ax.set_xticks(xticks)
    ax.set_xlim(xticks[0], xticks[-1])
    ax.set_ylim(yticks[0], yticks[-1])
#

def s_dyna(x, V):
    return x/((x/V)**3.+1.)

def plot_central(data, quantity, hue, ax, palette, plot_label=True, plot_legend=True, plot_sdyna=False, x_sdyna=30., remove_ylabel=False):
    # specify options
    if quantity=='last_y':
        infos={'x_label':'Wave celerity (µm/s)', 'y_label':r'$B_t$ at $t$=20min (mm)'}
    #
    
    if quantity=='dx':
        infos={'x_label':'Wave celerity (µm/s)', 'y_label':r'$\overline{D}_{\perp}$'}

        ## dx=0.5 ==> no bias
        ax.axhline(y=.5, linewidth = 2, color='#8C0335', linestyle='--')
        ## dx belongs to [0, 1]
        ax.set_ylim(bottom=0., top=1.)
    #

    if quantity=='vd':
        infos={'x_label':'Wave celerity', 'y_label':r'$\overline{v}$ (µm/s)'}
        
        ## plot the theoretical expression of s_dyna
        if plot_sdyna:
            x=np.linspace(0., x_sdyna, 70)
            ax.plot(x, x/((x/4.)**3.+1.), color='black', linewidth=3., zorder=2)
        #
    #

    if quantity=='r_d':
        infos={'x_label':'Wave celerity (µm/s)', 'y_label':r'$\overline{r}_d$'}
        ax.axhline(y=0., linewidth = 2, color='#8C0335', linestyle='--')
    #

    ## remove y label of inset panels
    if remove_ylabel:
        ax.set_yticklabels([])
    #

    sns.scatterplot(data=data,
                    x='celerity',
                    y=quantity,
                    hue=hue,
                    ax=ax,
                    palette=palette,
                    legend=plot_legend,
                    s=150,
                    zorder=1
                    )
    #
    
    ax.tick_params(axis='both', which='major', labelsize=20)
    if plot_label:
        ax.set_xlabel(infos['x_label'], fontsize=25)
        ax.set_ylabel(infos['y_label'], fontsize=25)
    else:
        ax.set_xlabel(r'')
        ax.set_ylabel(r'')
    #

    if plot_legend:
        handles, labels = plt.gca().get_legend_handles_labels()
        mpatch1 = mpatches.Patch(facecolor=palette[1], label=r'1', edgecolor='black')
        mpatch2 = mpatches.Patch(facecolor=palette[0], label=r'0', edgecolor='black')
        handles.extend([mpatch2, mpatch1])
        ax.legend(handles=handles, fontsize=25, title=r"$q$ (µm$^2$/s)", title_fontsize=28, loc='upper right')
        ax.tick_params(axis='both', which='major', labelsize=20)
#

def plot_critical(data, y, param, r2, ax, palette, plot_legend=False, xticks = [2.0, 4.0, 6.0], alpha=.5, th=(False, .79, 'th')):
    marker={'last_y':'^', 'dx':'o', 'vd':'s'}
    ylabel={'celerity':"Critical wave celerity (µm/s)", 'vd':r"Maximum $\overline{v}$ (µm/s)"}
    
    ## plot theoretical expression
    if th[0]:
        x=np.linspace(2., 6., 3)
        ax.plot(x, th[1]*x, color='black', zorder=1, linestyle='--', linewidth=5.)
    #
    
    ## plot linear interpolation
    sns.lineplot(data=data,
                    x=param,
                    y=y,
                    hue='wavelength',
                    ax=ax,
                    palette=palette,
                    linewidth=3.,
                    legend=plot_legend,
                    zorder=2
                    )
    #

    ## plot points
    sns.scatterplot(data=data,
                    x=param,
                    y=y,
                    hue='wavelength',
                    ax=ax,
                    palette=palette,
                    linewidth=3.,
                    marker=marker[data.columns[-1]],
                    s=200,
                    legend=plot_legend,
                    zorder=3
                    )
    #

    ## reduce alpha of points
    for coll in ax.collections:
        facecolors = coll.get_facecolors()
        facecolors[:, -1] = alpha  # % de transparence
        coll.set_facecolors(facecolors)
    #


    ## write R^2 values
    for i, val in enumerate(r2):
        txt = f"$R^2$ = {val:.2f}"  
        ax.text(.7, .2 - i * 0.05, txt, transform=ax.transAxes, fontsize=21, va='top', ha='left', color=palette[i])
    #
    
    ## axis
    ax.tick_params(axis='both', which='major', labelsize=20)
    ax.set_xlabel("$V_{max}$ (µm/s)", fontsize=25)
    ax.set_ylabel(ylabel[y], fontsize=25)
    ax.set_xticks(xticks)

    ## legend
    if plot_legend:
        handles, labels = plt.gca().get_legend_handles_labels()
        mpatch1 = mpatches.Patch(facecolor=palette[0], label=r'20', edgecolor='black')
        mpatch2 = mpatches.Patch(facecolor=palette[1], label=r'40', edgecolor='black')
        mpatch3 = mpatches.Patch(facecolor=palette[2], label=r'80', edgecolor='black')
        point1 = Line2D([0], [0], label=th[2], markersize=10, linestyle='--', color='black', linewidth=3.)
        ax.legend(handles=[mpatch1, mpatch2, mpatch3, point1], fontsize=23, title='Wavelength (µm)', title_fontsize=25)
        ax.tick_params(axis='both', which='major', labelsize=20)
    #
#


# ================================================================================================================
# ================================================================================================================
#                                                    User interface
# ================================================================================================================
# ================================================================================================================

## loading data
must_be_int = ['Surface Index', 'Number of Cells']

substrate_regime="dyna"
study="nonDiagonal_rotation"




## x and y axis
xticks = [-2., 0., 2.]
yticks = [-0.5, 0., 3.5]

mapping_traj={1:ax1, 2:ax2, 3:ax3, 4:ax4, 5:ax5, 6:ax6, 7:ax7, 8:ax8}
Nb_to_plot=200


################### compute r2 score
perform_stat = False
dependant_variable='vd'

################### export figure or not
save = True
path_to_save = substrate_regime + "\\assets\\" + study + "\\" + "test.svg"


################### platte de couleur 
wl_palette = ['#601A4A',  # wavelength=20 µm
              '#EE442F',  # 40 µm
              '#63ACBE'   # 80 µm
              ]
#
q_palette = ["#986eac", "#601A4A"]
assumption_palette = ['#E1BE6A', '#40B0A6']


# ================================================================================================================
# ================================================================================================================
#                                                    MAIN function
# ================================================================================================================
# ================================================================================================================

root = os.getcwd()
list_of_paths=search_files_to_analyse(substrate_regime + "\\assets\\" + study)

count=0
for path in list_of_paths:
    count+=1

    print("\nload file nb", count)
    print(path)
    print("------------------------------")

    metadata, data=load_trajectory_data(path, must_be_int)

    # conversion µm -> mm and reverse y-axis
    data["y"] /= -1000.
    data["x"] /= 1000.

    plot_traj(metadata, data, Nb_to_plot, mapping_traj[count], xticks, yticks)
#

    



if save:
    print('\n------------------------------------ export figure')
    plt.savefig(path_to_save)