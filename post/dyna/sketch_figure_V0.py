# ================================================================================================================
# ================================================================================================================
#                                                    packages importation
# ================================================================================================================
# ================================================================================================================

import numpy as np

import os
import regex as re
import pandas as pd

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.gridspec import GridSpec
import matplotlib.patches as mpatches
import seaborn as sns

# ================================================================================================================
# ================================================================================================================
#                                                    open figure
# ================================================================================================================
# ================================================================================================================

######## open figure

NbLignes, NbColonnes=4, 4
Largeur, Hauteur = 30, 30
Largeur_ins, Hauteur_ins = 2.3, 2.3
fig1 = plt.figure(num='data_analysis_dynamic_curvotaxis', figsize=(Largeur, Hauteur), constrained_layout=True)
gs1 = GridSpec(NbLignes, NbColonnes, figure=fig1)

## axes pour plot les trajectoires / barycentres
ax_traj1 = fig1.add_subplot(gs1[0, 0])
ax_traj2 = fig1.add_subplot(gs1[0, 1])
ax_traj3 = fig1.add_subplot(gs1[0, 2])
ax_traj4 = fig1.add_subplot(gs1[0, 3])

ax_traj5 = fig1.add_subplot(gs1[1, 0])
ax_traj6 = fig1.add_subplot(gs1[1, 1])
ax_traj7 = fig1.add_subplot(gs1[1, 2])
ax_traj8 = fig1.add_subplot(gs1[1, 3])

## axe pour captude d'écran d'une simulation
ax_fig = fig1.add_subplot(gs1[2, 0])

# axes pour plot les autres indicateurs (Dx, velocity etc)
ax1 = fig1.add_subplot(gs1[2, 1])
ax2 = fig1.add_subplot(gs1[2, 2])
ax3 = fig1.add_subplot(gs1[2, 3])

ax4 = fig1.add_subplot(gs1[3, 0])
ax5 = fig1.add_subplot(gs1[3, 1])
ax6 = fig1.add_subplot(gs1[3, 2])
ax7 = fig1.add_subplot(gs1[3, 3])

## add labels
def add_label(ax, txt, x=-0.05, y=1.1):
    ax.text(x, y, txt, transform=ax.transAxes, fontsize=20, fontweight='bold', va='top', ha='right')
#

add_label(ax_traj1, 'A.1')
add_label(ax_traj2, 'A.2')
add_label(ax_traj3, 'A.3')
add_label(ax_traj4, 'A.4')

add_label(ax_traj5, 'B.1')
add_label(ax_traj6, 'B.2')
add_label(ax_traj7, 'B.3')
add_label(ax_traj8, 'B.4')

add_label(ax_fig, 'C')

add_label(ax1, 'D')
add_label(ax2, 'E')
add_label(ax3, 'F')
add_label(ax4, 'G')
add_label(ax5, 'H')
add_label(ax6, 'I')
add_label(ax7, 'J')

axins1 = inset_axes(ax_traj1, width=Largeur_ins, height=Hauteur_ins)
axins2 = inset_axes(ax_traj2, width=Largeur_ins, height=Hauteur_ins)
axins3 = inset_axes(ax_traj3, width=Largeur_ins, height=Hauteur_ins)
axins4 = inset_axes(ax_traj4, width=Largeur_ins, height=Hauteur_ins)

## font size
plt.rcParams.update({'font.size': 50})


# ================================================================================================================
# ================================================================================================================
#                                                    intermediate functions
# ================================================================================================================
# ================================================================================================================


def plot_traj(metadata, data, IDlist, mapping, xticks, yticks):

    # select axe to plot
    wl, c = metadata[0], metadata[1]
    ax = mapping[(wl, c)]
    
    sns.scatterplot(data=data[data['cell'].isin(IDlist)], 
                    x='x', 
                    y='y', 
                    hue='cell', 
                    ax=ax, 
                    palette='deep', 
                    marker='o', 
                    edgecolor='none'
                    )
    
    ax.get_legend().set_visible(False)
    ax.set_xlabel('x (µm)', fontsize=17)
    ax.set_ylabel('y (µm)', fontsize=17)

    ax.text(0.8, 1.1, 
            r'$\lambda$={}µm, $c$={}µm/s'.format(wl, c), 
            transform=ax.transAxes, 
            fontsize=20, fontweight='bold', va='top', ha='right'
            )
    #

    ax.set_yticks(yticks) 
    ax.set_xticks(xticks)
    ax.set_xlim(xticks[0], xticks[-1])
    ax.set_ylim(yticks[0], yticks[-1])
#


def plot_baryc(metadata, data, list_frame, mapping, xticks, yticks):

    # select axe to plot
    wl, c, dt = metadata[0], metadata[1], metadata[2]
    axins, plot_legend = mapping[(wl, c)]
    palette = ["#CA47F1", "#9D70FF", "#5EB3FD"]

    pop = data[data['frame'].isin(list_frame)]
    sns.kdeplot(data=pop,
                x='x', 
                y='y',
                hue='frame',
                ax=axins, 
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
    axins.scatter(baryc['x'], baryc['y'], marker='o', s=150, c=palette, alpha = 1., edgecolors='black', zorder=2)

    if plot_legend:
        # legend and labels
        handles, labels = plt.gca().get_legend_handles_labels()
        items = []
        for i in range(len(list_frame)):
            items.append(mpatches.Patch(facecolor=palette[i], label=float(list_frame[i])*dt/60., edgecolor='black'))
        #
        handles.extend(items)
        axins.legend(handles = handles, fontsize=15, title='t (min)', title_fontsize=15)
    #

    axins.set_xlabel('x (µm)', fontsize=20)
    axins.set_ylabel('y (µm)', fontsize=20) 

    axins.set_yticks(yticks) 
    axins.set_xticks(xticks)
    axins.set_xlim(xticks[0], xticks[-1])
    axins.set_ylim(yticks[0], yticks[-1])
#

def plot_bary_by_c(data, ax):
    sns.scatterplot(data=data,
                    x='celerity',
                    y='last_y',
                    hue='wavelength',
                    ax=ax,
                    palette=['#F24141', '#F29F05', '#09A603'])
    #
    
    if source=='cell_is_not_surfer':
        ax.set_xlim(left=0., right=10.)
    #
    ax.tick_params(axis='both', which='major', labelsize=20)
    ax.set_xlabel('Wave celerity (µm/s)', fontsize=25)
    ax.set_ylabel('Population barycenter (µm)', fontsize=25)
    ax.legend(fontsize=15, title='Wavelength (µm)', title_fontsize=15)
#

def plot_dist(data, quantity, c_list, ax1, plot_legend=True):

    # specify options
    if quantity=='dx':
        infos={'y_label1':r'$D_x$ distribution',
               'y_label2':r'$\overline{D}_x$',
        }

        ## dx=0.5 ==> no bias
        ax1.axhline(y=.5, linewidth = 2, color='#8C0335', linestyle='--')
        
        ## dx belongs to [0, 1]
        ax1.set_ylim(bottom=0., top=1.)
    #

    if quantity=='vd':
        infos={'y_label1':'Cell velocity (µm/s)',
               'y_label2':'Median cell velocity (µm/s)',
        }
        ax1.set_ylim(bottom=0., top=5.)
    #

    if quantity=='r_d':
        infos={'y_label1':r'$r_d$ distribution',
               'y_label2':r'$\overline{r}_d$',
        }
        ax1.axhline(y=0., linewidth = 2, color='#8C0335', linestyle='--')
    #
        
    # plot distributions
    sns.violinplot(data=data[data['celerity'].isin(c_list)], 
                   x='wavelength', 
                   y=quantity, 
                   hue='celerity',
                   bw_adjust=.5, 
                   cut=1,
                   linewidth=1, 
                   density_norm='area', 
                   palette=['#4F1167', '#2D568C', '#008C8D', '#FFE529'], 
                   ax=ax1,
                   legend=plot_legend
                   )
    #

    # labels fig 1
    ax1.set_ylabel(infos['y_label1'], fontsize=25)
    ax1.set_xlabel(r'Wavelength (µm)', fontsize=25)
    if plot_legend:
        ax1.legend(fontsize=15, title='Wave celerity (µm/s)', title_fontsize=15)
    #

    ax1.tick_params(axis='both', which='major', labelsize=20)    
#

def plot_central(data, quantity, ax1, plot_legend=True):
    # specify options
    if quantity=='dx':
        infos={'y_label1':r'$D_x$ distribution',
               'y_label2':r'$\overline{D}_x$',
        }

        ## dx=0.5 ==> no bias
        ax1.axhline(y=.5, linewidth = 2, color='#8C0335', linestyle='--')
        
        ## dx belongs to [0, 1]
        ax1.set_ylim(bottom=0., top=1.)

    if quantity=='vd':
        infos={'y_label1':'Cell velocity (µm/s)',
               'y_label2':'Median cell velocity (µm/s)',
        }

        ax1.set_ylim(top=3.)
    if quantity=='r_d':
        infos={'y_label1':r'$r_d$ distribution',
               'y_label2':r'$\overline{r}_d$',
        }

        ax1.axhline(y=0., linewidth = 2, color='#8C0335', linestyle='--')
    #

    if source=='cell_is_not_surfer':
        ax1.set_xlim(right=5.)
    #

    sns.scatterplot(data=data,
                    x='celerity',
                    y=quantity,
                    hue='wavelength',
                    ax=ax1,
                    palette=['#F24141', '#F29F05', '#09A603'],
                    legend=plot_legend
                    )
    
    ax1.tick_params(axis='both', which='major', labelsize=20)
    ax1.set_xlabel('Wave celerity (µm/s)', fontsize=25)
    ax1.set_ylabel(infos['y_label2'], fontsize=25)
    ax1.legend(fontsize=15, title='Wavelength (µm)', title_fontsize=15)
#



# ================================================================================================================
# ================================================================================================================
#                                                    User interface
# ================================================================================================================
# ================================================================================================================

## source    cell_is_surfer  /  cell_is_not_surfer  /  role_of_Vmax\\Vmax=...\\cell_is_surfer
source = 'cell_is_not_surfer'

## nb oc cell trak to plot
NtoPlot=30

## metadada
wavelength = 20.
dt = 15.

## list of wave celerity to plot distributions
#if source=='cell_is_not_surfer'
c_list = [.1, .6, 1., 10.]

if source=='cell_is_surfer':
    c_list = [.1, 1., 3., 10.]
#

## list of frames to plot
t_list = [10, 40, 80]

# mapping between simulation parameters and ax to plot traj
ax_mapping_traj = {
    (20., .1):ax_traj1,
    (20., .6):ax_traj2,
    (20., 1.):ax_traj3,
    (20., 10.):ax_traj4
}

ax_mapping_baryc = {
    (20., .1): (axins1, True),
    (20., .6): (axins2, False),
    (20., 1.): (axins3, False),
    (20., 10.):(axins4, False)
}

# export figure or not
save = True

### ticks for traj plot

### if source=='cell_is_not_surfer'
xticks = [-1000., -500., 0., 500., 1000.]
yticks = [-300., 0., 1000., 1700.]

if source=='cell_is_surfer':
    xticks = [-1900., -1000., 0., 1000., 1900.]
    yticks = [-300., 0., 1000., 2000., 3500.]
#


# ================================================================================================================
# ================================================================================================================
#                                                    MAIN function
# ================================================================================================================
# ================================================================================================================


print('\n------------------------------------ plot trajectories and barycenters')

for root_folder, sub_folder, files in os.walk(os.getcwd() + '\\' + source + '\\analysis_results\\traj'):
    for file in files:
        c = float(re.findall(r'c=(\d*.\d*).csv', file)[0])
        metadata = [wavelength, c, dt]
        traj_df = pd.read_csv(os.path.join(root_folder, file))

        print('...... w='+str(wavelength)+'µm, c='+str(c)+'µm/s')
        print('... trajectories')
        plot_traj(metadata, traj_df, np.random.choice(traj_df['cell'].unique(), size=NtoPlot, replace=False), ax_mapping_traj, xticks, yticks)

        print('... barycenters')
        plot_baryc(metadata, traj_df, t_list, ax_mapping_baryc, xticks, yticks)
    #
#

print('------------------------------------ final plots:')

## dx
print('... dx')
plot_dist(pd.read_csv(os.getcwd() + '\\' + source + '\\analysis_results\\global\\dx.csv'), 'dx', c_list, ax1=ax1)
plot_central(pd.read_csv(os.getcwd() + '\\' + source + '\\analysis_results\\global\\central_dx.csv'), 'dx', ax1=ax5)

print('... vd')
plot_dist(pd.read_csv(os.getcwd() + '\\' + source + '\\analysis_results\\global\\vd.csv'), 'vd', c_list, ax1=ax2)
plot_central(pd.read_csv(os.getcwd() + '\\' + source + '\\analysis_results\\global\\central_vd.csv'), 'vd', ax1=ax6)

print('... rd')
data = pd.read_csv(os.getcwd() + '\\' + source + '\\analysis_results\\global\\last.csv')
data_central = pd.read_csv(os.getcwd() + '\\' + source + '\\analysis_results\\global\\central_last.csv')
plot_dist(data[['wavelength', 'celerity', 'r_d']], 'r_d', c_list, ax1=ax3)
plot_central(data_central[['wavelength', 'celerity', 'r_d']], 'r_d', ax1=ax7)

print('... barycenter')
plot_bary_by_c(data_central[['wavelength', 'celerity', 'last_y']], ax = ax4)



if save:
    print('------------------------------------ export figure:')
    plt.savefig(source + '\\traj_analysis.svg')