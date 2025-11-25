# ================================================================================================================
# ================================================================================================================
#                                                    packages importation
# ================================================================================================================
# ================================================================================================================

import os
import pandas as pd

import numpy as np

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.gridspec import GridSpec
import seaborn as sns


# ================================================================================================================
# ================================================================================================================
#                                                    open figure
# ================================================================================================================
# ================================================================================================================

######## open figure

NbLignes, NbColonnes=2, 4
Largeur, Hauteur = 30, 15
Largeur_ins, Hauteur_ins = 2.5, 2.5

fig1 = plt.figure(num='effect_of_wavelength', figsize=(Largeur, Hauteur), constrained_layout=True)
gs1 = GridSpec(NbLignes, NbColonnes, figure=fig1)

# axes pour NDSC cells
ax1 = fig1.add_subplot(gs1[0, 0])
ax2 = fig1.add_subplot(gs1[0, 1])
ax3 = fig1.add_subplot(gs1[0, 2])
ax4 = fig1.add_subplot(gs1[0, 3])

# axes pour DSC cells
ax5 = fig1.add_subplot(gs1[1, 0])
ax6 = fig1.add_subplot(gs1[1, 1])
ax7 = fig1.add_subplot(gs1[1, 2])
ax8 = fig1.add_subplot(gs1[1, 3])

## axes d'insert 
ax_ins1 = inset_axes(ax1, width=Largeur_ins, height=Hauteur_ins)
ax_ins2 = ax2.inset_axes([.017, 0.06, .4, .4])
ax_ins3 = inset_axes(ax3, width=Largeur_ins, height=Hauteur_ins)
ax_ins4 = ax4.inset_axes([.59, 0.057, .4, .4])

ax_ins5 = inset_axes(ax5, width=Largeur_ins, height=Hauteur_ins)
ax_ins6 = ax6.inset_axes([0.08, 0.59, .4, .4])
ax_ins7 = inset_axes(ax7, width=Largeur_ins, height=Hauteur_ins)
ax_ins8 = ax8.inset_axes([.59, 0.57, .4, .4])

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
add_label(ax7, 'G')
add_label(ax8, 'H')

## font size
plt.rcParams.update({'font.size': 50})

# ================================================================================================================
# ================================================================================================================
#                                                    intermediate functions
# ================================================================================================================
# ================================================================================================================

def plot_bary_by_c(data, ax, palette, plot_label=True, plot_legend=True, plot_sdyna=False):

    sns.scatterplot(data=data,
                    x='celerity',
                    y='last_y',
                    hue='wavelength',
                    ax=ax,
                    palette=palette,
                    s=150,
                    legend=plot_legend
    )
    
    if ax == ax1:
        y0, y1 = ax.get_ylim()
        ax.text(-10., .5*(y0+y1), 
                r"NDSC cells $\bf f_{dyna}=0$", 
                rotation = 'vertical', fontsize=30, fontweight='bold', va='center', ha='center'
                )
        #
    #

    if ax == ax5:
        y0, y1 = ax.get_ylim()
        ax.text(-10., .5*(y0+y1), 
                r"DSC cells $\bf f_{dyna}\neq0$", 
                rotation = 'vertical', fontsize=30, fontweight='bold', va='center', ha='center'
                )
        #
    #
    
    ax.tick_params(axis='both', which='major', labelsize=20)
    if plot_label:
        ax.set_xlabel('Wave celerity (µm/s)', fontsize=25)
        ax.set_ylabel(r'$B_T$ with $T$=20min (mm)', fontsize=25)
    else:
        ax.set_xlabel(r'')
        ax.set_ylabel(r'')
    #

    if plot_legend:
        handles, labels = plt.gca().get_legend_handles_labels()
        mpatch1 = mpatches.Patch(facecolor=palette[0], label=r'20', edgecolor='black')
        mpatch2 = mpatches.Patch(facecolor=palette[1], label=r'40', edgecolor='black')
        mpatch3 = mpatches.Patch(facecolor=palette[2], label=r'80', edgecolor='black')
        handles.extend([mpatch1, mpatch2, mpatch3])
        ax.legend(handles=handles, fontsize=20, title=r"Wavelength (µm)", title_fontsize=22)
        ax.tick_params(axis='both', which='major', labelsize=20)
#

def plot_central(data, quantity, ax, palette, plot_label=True, plot_legend=True, plot_sdyna=False, x_sdyna=30., remove_ylabel=False):
    # specify options
    if quantity=='dx':
        infos={'y_label':r'$\overline{D}_{\perp}$'}

        ## dx=0.5 ==> no bias
        ax.axhline(y=.5, linewidth = 2, color='#8C0335', linestyle='--')
        
        ## dx belongs to [0, 1]
        ax.set_ylim(bottom=0., top=1.)

    if quantity=='vd':
        infos={'y_label':r'$\overline{v}$ (µm/s)'}
    #

    if quantity=='r_d':
        infos={'y_label':r'$\overline{r}_d$'}
        ax.axhline(y=0., linewidth = 2, color='#8C0335', linestyle='--')
    #

    if plot_sdyna:
        x=np.linspace(0., x_sdyna, 70)
        ax.plot(x, x/((x/4.)**3.+1.), color='black', linewidth=3.)
    #

    if remove_ylabel:
        ax.set_yticklabels([])
    #

    sns.scatterplot(data=data,
                    x='celerity',
                    y=quantity,
                    hue='wavelength',
                    ax=ax,
                    palette=palette,
                    legend=plot_legend,
                    s=150
                    )
    #
    
    ax.tick_params(axis='both', which='major', labelsize=20)
    if plot_label:
        ax.set_xlabel('Wave celerity (µm/s)', fontsize=25)
        ax.set_ylabel(infos['y_label'], fontsize=25)
    else:
        ax.set_xlabel(r'')
        ax.set_ylabel(r'')
    #

    if plot_legend:
        handles, labels = plt.gca().get_legend_handles_labels()
        mpatch1 = mpatches.Patch(facecolor=palette[0], label=r'20', edgecolor='black')
        mpatch2 = mpatches.Patch(facecolor=palette[1], label=r'40', edgecolor='black')
        mpatch3 = mpatches.Patch(facecolor=palette[2], label=r'80', edgecolor='black')
        handles.extend([mpatch1, mpatch2, mpatch3])
        ax.legend(handles=handles, fontsize=20, title=r"Wavelength (µm)", title_fontsize=22, loc=4)
        ax.tick_params(axis='both', which='major', labelsize=20)
#

# ================================================================================================================
# ================================================================================================================
#                                                    User interface
# ================================================================================================================
# ================================================================================================================

## metadada
assumptions = ['cell_is_not_surfer', 'cell_is_surfer']
wavelengths = [20., 40., 80.]
max_celerity = plot_or_not={'cell_is_not_surfer':2., 'cell_is_surfer':5.}
dt = 15.

# mapping between assumption and quantity and axes
mapping={('cell_is_not_surfer', 'b'):ax1,
         ('cell_is_not_surfer', 'dx'):ax2,
         ('cell_is_not_surfer', 'vd'):ax3,
         ('cell_is_not_surfer', 'r_d'):ax4,
         ('cell_is_surfer', 'b'):ax5,
         ('cell_is_surfer', 'dx'):ax6,
         ('cell_is_surfer', 'vd'):ax7,
         ('cell_is_surfer', 'r_d'):ax8
         }

mapping_ins={('cell_is_not_surfer', 'b'):ax_ins1,
             ('cell_is_not_surfer', 'dx'):ax_ins2,
             ('cell_is_not_surfer', 'vd'):ax_ins3,
             ('cell_is_not_surfer', 'r_d'):ax_ins4,
             ('cell_is_surfer', 'b'):ax_ins5,
             ('cell_is_surfer', 'dx'):ax_ins6,
             ('cell_is_surfer', 'vd'):ax_ins7,
             ('cell_is_surfer', 'r_d'):ax_ins8
             }
#

# plot sdyna theorical expression
plot_or_not={'cell_is_not_surfer':False, 'cell_is_surfer':True}

# export figure or not
save = True
path_to_save = 'wavelength_effect.svg'


## platte de couleur 
wl_palette = ['#601A4A',  # wavelength=20 µm
              '#EE442F',  # 40 µm
              '#63ACBE'   # 80 µm
              ] 
#

# ================================================================================================================
# ================================================================================================================
#                                                    MAIN function
# ================================================================================================================
# ================================================================================================================

root = os.getcwd()

for assumption in assumptions:
    print('\n----------------------------------- treatment of ' + assumption)
    print('... open data')
    path_to_last = root + '\\' + assumption + '\\analysis_results\\global\\central_last.csv'
    path_to_dx =   root + '\\' + assumption + '\\analysis_results\\global\\central_dx.csv'
    path_to_vd =   root + '\\' + assumption + '\\analysis_results\\global\\central_vd.csv'

    data_last = pd.read_csv(path_to_last)
    data_dx =   pd.read_csv(path_to_dx)
    data_vd =   pd.read_csv(path_to_vd)

    print('... sort by wavelength')
    data_last = data_last[data_last['wavelength'].isin(wavelengths)]
    data_dx =   data_dx[data_dx['wavelength'].isin(wavelengths)]
    data_vd =   data_vd[data_vd['wavelength'].isin(wavelengths)]

    # µm -->  mm
    data_last['last_y'] /= 1000. 

    if assumption=='cell_is_surfer':
        for wl in data_last['wavelength'].unique():
            f=data_last[data_last['wavelength']==wl]['r_d']
            dx=data_last[data_last['wavelength']==wl]['celerity']
            grad=np.gradient(f, dx)
            idx=np.where(grad==grad.min())[0]
            print('... min of gradient of rd\n', data_last.iloc[idx])
            

    print('... plot general indicators')
    plot_bary_by_c(data_last, 
                   ax=mapping[(assumption, 'b')], 
                   palette=wl_palette,
                   plot_legend=False 
                   )
    
    plot_central(data_dx, 'dx', 
                 ax=mapping[(assumption, 'dx')], 
                 palette=wl_palette, 
                 plot_legend=True,
                 )
    
    plot_central(data_vd, 'vd', 
                 ax=mapping[(assumption, 'vd')], 
                 palette=wl_palette, 
                 plot_legend=False,
                 plot_sdyna=plot_or_not[assumption]
                 )
    
    plot_central(data_last, 'r_d', 
                 ax=mapping[(assumption, 'r_d')], 
                 palette=wl_palette, 
                 plot_legend=False
                 )
    
    print('... plot insets')
    # sort by celerity
    data_last = data_last[data_last['celerity'] <= max_celerity[assumption]]
    data_dx =   data_dx[data_dx['celerity'] <= max_celerity[assumption]]
    data_vd =   data_vd[data_vd['celerity'] <= max_celerity[assumption]]

    plot_bary_by_c(data_last, 
                   ax=mapping_ins[(assumption, 'b')], 
                   palette=wl_palette, 
                   plot_label=False, 
                   plot_legend=False,
                   )
    
    plot_central(data_dx, 'dx', 
             ax = mapping_ins[(assumption, 'dx')], 
             palette=wl_palette,
             plot_label=False, 
             plot_legend=False,
             remove_ylabel=True
             )
    
    plot_central(data_vd, 'vd', 
                 ax = mapping_ins[(assumption, 'vd')], 
                 palette=wl_palette,
                 plot_label=False, 
                 plot_legend=False,
                 plot_sdyna=plot_or_not[assumption],
                 x_sdyna=max_celerity[assumption]
                 )
    
    plot_central(data_last, 'r_d', 
                 ax = mapping_ins[(assumption, 'r_d')], 
                 palette=wl_palette,
                 plot_label=False, 
                 plot_legend=False
                 )
#


if save:
    print('------------------------------------ export figure:')
    plt.savefig(root + '\\' +  path_to_save)
#