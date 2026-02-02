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

def plot_bary_by_c(data, y, ax, palette, plot_label=True, plot_legend=True, plot_sdyna=False):

    sns.scatterplot(data=data,
                    x='celerity',
                    y=y,
                    hue='wavelength',
                    ax=ax,
                    palette=palette,
                    s=150,
                    legend=plot_legend
    )
    
    if ax == ax1:
        y0, y1 = ax.get_ylim()
        ax.text(-.25, .5, 
                r"NDSC cells $\bf f_{dyna}=0$", 
                rotation = 'vertical', fontsize=30, fontweight='bold', va='center', ha='center', transform=ax.transAxes
                )
        #
    #

    if ax == ax5:
        y0, y1 = ax.get_ylim()
        ax.text(-.25, .5, 
                r"DSC cells $\bf f_{dyna}\neq0$", 
                rotation = 'vertical', fontsize=30, fontweight='bold', va='center', ha='center', transform=ax.transAxes
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

    ax.set_ylim(top=.9)

    if plot_legend:
        handles, labels = plt.gca().get_legend_handles_labels()
        mpatch1 = mpatches.Patch(facecolor=palette[0], label=r'20', edgecolor='black')
        mpatch2 = mpatches.Patch(facecolor=palette[1], label=r'40', edgecolor='black')
        mpatch3 = mpatches.Patch(facecolor=palette[2], label=r'80', edgecolor='black')
        handles.extend([mpatch1, mpatch2, mpatch3])
        ax.legend(handles=handles, fontsize=20, title=r"Wavelength (µm)", title_fontsize=22)
        ax.tick_params(axis='both', which='major', labelsize=20)
#

def plot_central(data, y, quantity, ax, palette, plot_label=True, plot_legend=True, plot_sdyna=False, x_sdyna=4., n=3., remove_ylabel=False):
    
    # specify options
    if quantity=='dx':
        infos={'y_label':r'$\overline{D}_{\perp}$'}

        ## dx=0.5 ==> no bias
        ax.axhline(y=.5, linewidth = 2, color='#8C0335', linestyle='--')
        
        ## dx belongs to [0, 1]
        ax.set_ylim(bottom=0., top=1.)

    if quantity=='vd':
        infos={'y_label':r'$\overline{v}$ (µm/s)'}
        ax.set_ylim(bottom=0.)
    #

    if quantity=='rd':
        infos={'y_label':r'$\overline{r}_d$'}
        ax.axhline(y=0., linewidth = 2, color='#8C0335', linestyle='--')
        ax.set_ylim(top=1., bottom=-0.1)
    #

    if plot_sdyna:
        Vm={80.:0.5, 40.:.56,20.:0.6}
        x=np.linspace(0., x_sdyna, 70)
        count=0
        for L in data["wavelength"].unique():
            ax.axvline(x=(n-1.)**(-1./n), color=palette[count])
            #ax.plot(x, x/((x/Vm[L])**n+1.), linewidth=3., color=palette[count])
            count+=1
    #

    if remove_ylabel:
        ax.set_yticklabels([])
    #

    sns.scatterplot(data=data,
                    x='celerity',
                    y=y,
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
    #ax.set_xlim(left=0., right=3.)

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

# study
substrate_regime="dyna"  # dyna or static
NDSC_study="NDSC_RightSinus_20min_VariousCelerity"
DSC_study="DSC_RightSinus_20min_VariousCelerity"

# assets
NDSC_assets=os.path.join(os.getcwd(), substrate_regime, "assets", NDSC_study)
DSC_assets=os.path.join(os.getcwd(), substrate_regime, "assets", DSC_study)

mapping_path_assumption={"NDSC":NDSC_assets, "DSC":DSC_assets}

## metadada
assumptions = ['NDSC', 'DSC']
wavelengths = [20., 40., 80.]
max_celerity ={'NDSC':2.5, 'DSC':2.5}
dt = 15.

# mapping between assumption and quantity and axes
mapping={('NDSC', 'b'):ax1,
         ('NDSC', 'dx'):ax2,
         ('NDSC', 'vd'):ax3,
         ('NDSC', 'rd'):ax4,
         ('DSC', 'b'):ax5,
         ('DSC', 'dx'):ax6,
         ('DSC', 'vd'):ax7,
         ('DSC', 'rd'):ax8
         }

mapping_ins={('NDSC', 'b'):ax_ins4,
             ('NDSC', 'dx'):ax_ins2,
             ('NDSC', 'vd'):ax_ins3,
             ('NDSC', 'rd'):ax_ins1,
             ('DSC', 'b'):ax_ins8,
             ('DSC', 'dx'):ax_ins6,
             ('DSC', 'vd'):ax_ins7,
             ('DSC', 'rd'):ax_ins5
             }
#

# plot sdyna theorical expression
plot_or_not={'NDSC':False, 'DSC':True}

# export figure or not
save = True
path_to_save=os.path.join(os.getcwd(), substrate_regime, "wavelength_effect.svg")


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
    
    path_to_last = mapping_path_assumption[assumption] + '\\' + 'central_last.csv'
    path_to_dx =   mapping_path_assumption[assumption] + '\\' + 'central_dx.csv'
    path_to_vd =   mapping_path_assumption[assumption] + '\\' + 'central_vd.csv'

    data_last = pd.read_csv(path_to_last)
    data_dx =   pd.read_csv(path_to_dx)
    data_vd =   pd.read_csv(path_to_vd)

    # format columns
    data_last.rename(columns={"last_y":"last_y.mean", "last_y.1":"last_y.med", "last_y.2":"last_y.std", "r_d":"rd.mean", "r_d.1":"rd.med", "r_d.2":"rd.std"}, inplace=True)
    data_last.dropna(inplace=True)
    data_last[['last_y.mean', 'last_y.med', 'last_y.std', 'rd.mean', 'rd.med', 'rd.std']]=data_last[['last_y.mean', 'last_y.med', 'last_y.std', 'rd.mean', 'rd.med', 'rd.std']].astype("float")
    
    data_dx.rename(columns={"dx":"dx.mean", "dx.1":"dx.med", "dx.2":"dx.std"}, inplace=True)
    data_dx.dropna(inplace=True)
    data_dx[['dx.mean', 'dx.med', 'dx.std']]=data_dx[['dx.mean', 'dx.med', 'dx.std']].astype("float")

    data_vd.rename(columns={"vd":"vd.mean", "vd.1":"vd.med", "vd.2":"vd.std"}, inplace=True)
    data_vd.dropna(inplace=True)
    data_vd[['vd.mean', 'vd.med', 'vd.std']]=data_vd[['vd.mean', 'vd.med', 'vd.std']].astype("float")

    print('... sort by wavelength and celerity')
    data_last = data_last[data_last['wavelength'].isin(wavelengths)]
    data_dx =   data_dx[data_dx['wavelength'].isin(wavelengths)]
    data_vd =   data_vd[data_vd['wavelength'].isin(wavelengths)]

    data_last = data_last[data_last['celerity'] <= max_celerity[assumption]]
    data_dx =   data_dx[data_dx['celerity'] <=     max_celerity[assumption]]
    data_vd =   data_vd[data_vd['celerity'] <=     max_celerity[assumption]]

    # µm -->  mm
    data_last[['last_y.mean', 'last_y.med', 'last_y.std']] /= 1000.

    #if assumption=='DSC':
    #    for wl in data_last['wavelength'].unique():
    #        f=data_last[data_last['wavelength']==wl]['rd.mean']
    #        dx=data_last[data_last['wavelength']==wl]['celerity']
    #        grad=np.gradient(f, dx)
    #        idx=np.where(grad==grad.min())[0]
    #        print('... min of gradient of rd\n', data_last.iloc[idx])
            

    print('... plot general indicators')
    plot_bary_by_c(data_last, y='last_y.mean',
                   ax=mapping[(assumption, 'b')], 
                   palette=wl_palette,
                   plot_legend=False 
                   )
    
    plot_central(data_dx, y='dx.mean', 
                 quantity='dx',
                 ax=mapping[(assumption, 'dx')], 
                 palette=wl_palette, 
                 plot_legend=True,
                 )
    
    plot_central(data_vd, y='vd.mean', 
                 quantity='vd',
                 ax=mapping[(assumption, 'vd')], 
                 palette=wl_palette, 
                 plot_legend=False,
                 plot_sdyna=plot_or_not[assumption],
                 x_sdyna=max_celerity[assumption],
                 n=3.
                 )
    
    plot_central(data_last, y='rd.mean', 
                 quantity='rd', 
                 ax=mapping[(assumption, 'rd')], 
                 palette=wl_palette, 
                 plot_legend=False
                 )
    
    print('... plot insets')
    ## sort by celerity
    #data_last = data_last[data_last['celerity'] <= max_celerity[assumption]]
    #data_dx =   data_dx[data_dx['celerity'] <= max_celerity[assumption]]
    #data_vd =   data_vd[data_vd['celerity'] <= max_celerity[assumption]]
#
    #plot_bary_by_c(data_last, y='last_y.mean',
    #               ax=mapping_ins[(assumption, 'b')], 
    #               palette=wl_palette, 
    #               plot_label=False, 
    #               plot_legend=False,
    #               )
    #
    #plot_central(data_dx, y='dx.mean',
    #             quantity='dx', 
    #            ax = mapping_ins[(assumption, 'dx')], 
    #            palette=wl_palette,
    #            plot_label=False, 
    #            plot_legend=False,
    #            remove_ylabel=True
    #            )
    #
    #plot_central(data_vd, y='vd.mean',
    #             quantity='vd', 
    #             ax = mapping_ins[(assumption, 'vd')], 
    #             palette=wl_palette,
    #             plot_label=False, 
    #             plot_legend=False,
    #             plot_sdyna=plot_or_not[assumption],
    #             x_sdyna=max_celerity[assumption]
    #             )
    #
    plot_central(data_last, y='rd.mean',
                 quantity='rd', 
                 ax = mapping_ins[(assumption, 'rd')], 
                 palette=wl_palette,
                 plot_label=False, 
                 plot_legend=False
                 )
#


if save:
    print('------------------------------------ export figure:')
    plt.savefig(path_to_save)
#