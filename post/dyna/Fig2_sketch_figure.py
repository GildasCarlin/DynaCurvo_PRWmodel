# ================================================================================================================
# ================================================================================================================
#                                                    packages importation
# ================================================================================================================
# ================================================================================================================

import numpy as np

import os
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
Largeur_ins, Hauteur_ins = 2.5, 2.5
fig1 = plt.figure(num='data_analysis_traj', figsize=(Largeur, Hauteur), constrained_layout=True)
gs1 = GridSpec(NbLignes, NbColonnes, figure=fig1)

## axes pour plot les trajectoires / barycentres en insert
ax_traj1 = fig1.add_subplot(gs1[0, 0])
ax_traj2 = fig1.add_subplot(gs1[0, 1])
ax_traj3 = fig1.add_subplot(gs1[0, 2])
ax_traj4 = fig1.add_subplot(gs1[0, 3])

ax_traj5 = fig1.add_subplot(gs1[1, 0])
ax_traj6 = fig1.add_subplot(gs1[1, 1])
ax_traj7 = fig1.add_subplot(gs1[1, 2])
ax_traj8 = fig1.add_subplot(gs1[1, 3])

## axes d'insert de barycentres
ax_ins1 = inset_axes(ax_traj1, width=Largeur_ins, height=Hauteur_ins)
ax_ins2 = inset_axes(ax_traj2, width=Largeur_ins, height=Hauteur_ins)
ax_ins3 = inset_axes(ax_traj3, width=Largeur_ins, height=Hauteur_ins)
ax_ins4 = inset_axes(ax_traj4, width=Largeur_ins, height=Hauteur_ins)

ax_ins5 = inset_axes(ax_traj5, width=Largeur_ins, height=Hauteur_ins)
ax_ins6 = inset_axes(ax_traj6, width=Largeur_ins, height=Hauteur_ins)
ax_ins7 = inset_axes(ax_traj7, width=Largeur_ins, height=Hauteur_ins)
ax_ins8 = inset_axes(ax_traj8, width=Largeur_ins, height=Hauteur_ins)

## axes pour plot les autres indicateurs (Dx, velocity etc)
ax1 = fig1.add_subplot(gs1[2, :2])
ax2 = fig1.add_subplot(gs1[2, 2:])
ax3 = fig1.add_subplot(gs1[3, :2])
ax4 = fig1.add_subplot(gs1[3, 2:])

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

add_label(ax1, 'C')
add_label(ax2, 'D')
add_label(ax3, 'E')
add_label(ax4, 'F')

## font size
plt.rcParams.update({'font.size': 50})


# ================================================================================================================
# ================================================================================================================
#                                                    intermediate functions
# ================================================================================================================
# ================================================================================================================

def plot_traj(metadata, data, N, mapping, xticks, yticks, center=True):

    ## random cells id to plot
    IDlist = np.random.choice(data['cell'].unique(), size=N, replace=False)

    # select axe to plot
    wl, c = metadata[0], metadata[1]
    ax = mapping[(wl, c)]
    
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

    ax.set_yticks(yticks) 
    ax.set_xticks(xticks)
    ax.set_xlim(xticks[0], xticks[-1])
    ax.set_ylim(yticks[0], yticks[-1])

    if mapping == ax_mapping_traj_NS:
        ax.text(.5, 1.05, 
            r'$\bf \lambda$={}µm, $\bf c$={}µm/s, $\bf T$=20min'.format(int(wl), c), 
            fontsize=20, fontweight='bold', va='center', ha='center', transform=ax.transAxes
            )
        
    if ax == ax_traj1:
        ax.text(-.25, .5, 
                r"NDSC cells $\bf f_{dyna}=0$", 
                rotation = 'vertical', fontsize=30, fontweight='bold', va='center', ha='center', transform=ax.transAxes
                )
        #

        ax.arrow(0.1, .3, 0., 0.4, width=.04, head_length=.1, color='red', transform=ax.transAxes, length_includes_head=True)
        ax.text(0.18, 0.5, 
                r"Wave propagation",
                color='red', rotation = 'vertical', fontsize=20, fontweight='bold', va='center', ha='center', transform=ax.transAxes
                )
        #

        ax.arrow(.2, 0.9, -0.15, 0., width=.04, head_length=.08, color='red', length_includes_head=True, transform=ax.transAxes)
        ax.arrow(.2, 0.9, 0.15, 0., width=.04, head_length=.08, color='red', length_includes_head=True, transform=ax.transAxes)
        ax.text(0.2, 0.8, 
                r"Grooves axis",
                color='red', rotation = 'horizontal', fontsize=20, fontweight='bold', va='center', ha='center', transform=ax.transAxes
                )
        #
    #

    if center:
        ax.scatter([0.0], [0.0], s=200, color='deeppink')
    else:
        rect=mpatches.Rectangle([-0.2, -0.2], width=0.4, height=0.4, ec='deeppink', ls="--", alpha=1., fill=0, lw=4, zorder=2)
        ax.add_patch(rect)

    if ax == ax_traj5:
        ax.text(-.25, .5, 
                r"DSC cells $\bf f_{dyna}\neq 0$", 
                rotation = 'vertical', fontsize=30, fontweight='bold', va='center', ha='center', transform=ax.transAxes
                )
        #
        
        ax.arrow(0.1, .3, 0., 0.4, width=.04, head_length=.1, color='red', transform=ax.transAxes, length_includes_head=True)
        ax.text(0.18, 0.5, 
                r"Wave propagation",
                color='red', rotation = 'vertical', fontsize=20, fontweight='bold', va='center', ha='center', transform=ax.transAxes
                )
        #

        ax.arrow(.2, 0.9, -0.15, 0., width=.04, head_length=.08, color='red', length_includes_head=True, transform=ax.transAxes)
        ax.arrow(.2, 0.9, 0.15, 0., width=.04, head_length=.08, color='red', length_includes_head=True, transform=ax.transAxes)
        ax.text(0.2, 0.8, 
                r"Grooves axis",
                color='red', rotation = 'horizontal', fontsize=20, fontweight='bold', va='center', ha='center', transform=ax.transAxes
                )
        #
#

def plot_baryc(metadata, data, list_frame, mapping, xticks, yticks, palette):

    # select axe to plot
    wl, c, dt = metadata[0], metadata[1], metadata[2]
    axins, plot_legend = mapping[(wl, c)]

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

    axins.set_xlabel('x (mm)', fontsize=20)
    axins.set_ylabel('y (mm)', fontsize=20) 

    axins.set_yticks(yticks) 
    axins.set_xticks(xticks)
    axins.set_xlim(xticks[0], xticks[-1])
    axins.set_ylim(yticks[0], yticks[-1])
    axins.tick_params(axis='both', which='major', labelsize=13)
#

def plot_dist(data, quantity, ax1, palette, whis = 1.5, plot_legend=True):

    # specify options
    if quantity=='dx':
        infos={'y_label1':r'$D_{\perp}$ ',
        }

        ## dx=0.5 ==> no bias
        ax1.axhline(y=.5, linewidth = 2, color='#8C0335', linestyle='--')
        ## dx belongs to [0, 1]
        ax1.set_ylim(bottom=0., top=1.)
        

    elif quantity=='vd':
        infos={'y_label1':'Cell speed (µm/s)',
        }
        ax1.set_ylim(bottom=0., top=4.)

    elif quantity=='r_d':
        infos={'y_label1':r'$r_d$ ',
        }
        ax1.axhline(y=0., linewidth = 2, color='#8C0335', linestyle='--')
        ax1.axhline(y=1., linewidth = 2, color='#8C0335', linestyle='--')
        #ax1.margins(y=.2)
        ax1.set_ylim(bottom=-1.6, top=1.5)
    
    elif quantity=='last_y':
        infos={'y_label1':r'Cell position at $T$=20min (mm)',
        }
        ax1.axhline(y=0., linewidth = 2, color='#8C0335', linestyle='--')
        ax1.set_ylim(bottom=-.25, top=1.5)
    #
        
    # plot distributions
    sns.boxplot(data=data, 
                   x='celerity', 
                   y=quantity, 
                   hue='assumption', 
                   linewidth=1, 
                   palette=palette, 
                   ax=ax1,
                   whis = whis,
                   legend=plot_legend
                   )
    #

    # labels fig 1
    ax1.set_ylabel(infos['y_label1'], fontsize=25)
    ax1.set_xlabel(r'Wave celerity (µm/s)', fontsize=25)
    if plot_legend:
        ax1.legend(fontsize=30, title='Assumption', title_fontsize=32, loc = 'upper left')
    #

    ax1.tick_params(axis='both', which='major', labelsize=20)    
#

def add_stat(ax, eff_size, sign=4*['**'], x=[.12, .37, .62, .87], y=4*[.955], size=0.1):
    for i in range(len(sign)):
        ax.text(x=x[i], y=y[i], s=sign[i] + "\n" + eff_size['ES'] + r"=" + str(eff_size['values'][i]), transform=ax.transAxes, fontsize=15, fontweight='bold', va='center', ha='center')
        
        hlen=.5*size
        ymin, ymax = ax.get_ylim()
        y_convert = ymin + (y[i]-.04) * (ymax - ymin)
        ax.axhline(y=y_convert, xmin=x[i]-hlen, xmax=x[i]+hlen, color='black', linewidth=2.)
    #
#

def format_traj(data, center=False, convert=True):
    # trajectory are centered at origin.
    if center:
        data[['x', 'y']] -= data.groupby('cell')[['x', 'y']].transform('first')
        #for i in data['cell'].unique():
        #    cell_i = data.loc[data['cell']==i,['x','y']].copy()
        #    data.loc[data['cell']==i,['x', 'y']] -= cell_i.iloc[0]
    #

    # js frame vs python frame: the y-axis is reversed
    data['y'] *= -1.0

    # convert µm to mm
    if convert:
        data[['x', 'y']] /= 1000.0
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

## nb oc cell trak to plot
NtoPlot=60

## metadada
wavelength = 20.
dt = 15.

## list of wave celerity to plot distributions
celerities = [.1, .7, .9, 5.]

## list of frames to plot barycenters
times = [10, 40, 80]

# mapping between simulation parameters and ax to plot traj
ax_mapping_traj_NS = {
    (wavelength, celerities[0]): ax_traj1,
    (wavelength, celerities[1]): ax_traj2,
    (wavelength, celerities[2]): ax_traj3,
    (wavelength, celerities[3]):ax_traj4
}

ax_mapping_traj_S = {
    (wavelength, celerities[0]): ax_traj5,
    (wavelength, celerities[1]): ax_traj6,
    (wavelength, celerities[2]): ax_traj7,
    (wavelength, celerities[3]):ax_traj8
}

ax_mapping_baryc_NS = {
    (wavelength, celerities[0]): (ax_ins1, True),
    (wavelength, celerities[1]): (ax_ins2, False),
    (wavelength, celerities[2]): (ax_ins3, False),
    (wavelength, celerities[3]):(ax_ins4, False)
}

ax_mapping_baryc_S = {
    (wavelength, celerities[0]): (ax_ins5, True),
    (wavelength, celerities[1]): (ax_ins6, False),
    (wavelength, celerities[2]): (ax_ins7, False),
    (wavelength, celerities[3]):(ax_ins8, False)
}

# export figure or not
save = True
path_to_save=os.path.join(os.getcwd(), substrate_regime, "traj_analysis.svg")

### ticks for traj plot
xticks_NS = [-1., 0., 1.]
yticks_NS = [-.25, 0., 1.75]

xticks_S = [-1., 0., 1.]
yticks_S = [-.25, 0., 1.75]
#

### palettes de couleur
palette_for_time = ["#3F2773", "#1F9E89", "#F2D027"]
palette_for_assumption = ['#E1BE6A', '#40B0A6']



# ================================================================================================================
# ================================================================================================================
#                                                    MAIN function
# ================================================================================================================
# ================================================================================================================

print('\n------------------------------------ get general indicators')
print('---------------------- read data')

#### NS
assumption = 'NDSC'
last_NS = pd.read_csv(NDSC_assets + '\\' + "last.csv")
last_NS = last_NS[last_NS['celerity'].isin(celerities)]
last_NS = last_NS[last_NS['wavelength']==wavelength]
last_NS.insert(0, 'assumption', len(last_NS) * [assumption])

dx_NS = pd.read_csv(NDSC_assets + '\\' + "dx.csv")
dx_NS = dx_NS[dx_NS['celerity'].isin(celerities)]
dx_NS = dx_NS[dx_NS['wavelength']==wavelength]
dx_NS.insert(0, 'assumption', len(dx_NS) * [assumption])

vd_NS = pd.read_csv(NDSC_assets + '\\' + "vd.csv")
vd_NS = vd_NS[vd_NS['celerity'].isin(celerities)]
vd_NS = vd_NS[vd_NS['wavelength']==wavelength]
vd_NS.insert(0, 'assumption', len(vd_NS) * [assumption])

#### S
assumption = 'DSC'
last_S = pd.read_csv(DSC_assets + '\\' + "last.csv")
last_S = last_S[last_S['celerity'].isin(celerities)]
last_S = last_S[last_S['wavelength']==wavelength]
last_S.insert(0, 'assumption', len(last_S) * [assumption])

dx_S = pd.read_csv(DSC_assets + '\\' + "dx.csv")
dx_S = dx_S[dx_S['celerity'].isin(celerities)]
dx_S = dx_S[dx_S['wavelength']==wavelength]
dx_S.insert(0, 'assumption', len(dx_S) * [assumption])

vd_S = pd.read_csv(DSC_assets + '\\' + "vd.csv")
vd_S = vd_S[vd_S['celerity'].isin(celerities)]
vd_S = vd_S[vd_S['wavelength']==wavelength]
vd_S.insert(0, 'assumption', len(vd_S) * [assumption])

######## concat
all_last = pd.concat([last_NS, last_S], ignore_index=True)
all_dx =   pd.concat([dx_NS, dx_S], ignore_index=True)
all_vd =   pd.concat([vd_NS, vd_S], ignore_index=True)

### convert µm to mm (barycenter position)
all_last['last_y'] /= 1000.

print('---------------------- plot last')
plot_dist(data=       all_last, 
          quantity=   'last_y',  
          ax1=        ax1, 
          palette=    palette_for_assumption,
          plot_legend=True
          )

print('---------------------- plot dx')
plot_dist(data=       all_dx, 
          quantity=   'dx', 
          ax1=        ax2, 
          palette=    palette_for_assumption,
          plot_legend=False
          )

print('---------------------- plot vd')
plot_dist(data=       all_vd, 
          quantity=   'vd', 
          ax1=        ax3, 
          palette=    palette_for_assumption,
          whis=       5.,
          plot_legend=True
          )

print('---------------------- plot r_d')
plot_dist(data=       all_last, 
          quantity=   'r_d',  
          ax1=        ax4, 
          palette=    palette_for_assumption,
          plot_legend=False
          )



print('\n------------------------------------ loop over celerities')
for celerity in celerities:
    metadata = [wavelength, celerity, dt]

    print('\n---------------------- get traj for c =',str(celerity), 'µm/s')

    ######## récupérér les trajectoires ...
    #### NS
    traj_NS=pd.read_csv(os.path.join(NDSC_assets, "L="+str(wavelength), "c="+str(celerity), "data.txt"), comment='#')

    #### S
    traj_S= pd.read_csv(os.path.join(DSC_assets, "L="+str(wavelength), "c="+str(celerity), "data.txt"), comment='#')

       
    print('---------------------- treat traj')
    format_traj(traj_NS, center=True)
    format_traj(traj_S, center=True)

    ######## scatter les traj
    print("---------------------- plot traj")
    #### NS
    plot_traj(metadata=metadata, 
              data=    traj_NS, 
              N=       NtoPlot, 
              mapping= ax_mapping_traj_NS, 
              xticks=  xticks_NS, 
              yticks=  yticks_NS
              )
    
    #### S
    plot_traj(metadata=metadata, 
              data=    traj_S, 
              N=       NtoPlot, 
              mapping= ax_mapping_traj_S, 
              xticks=  xticks_S, 
              yticks=  yticks_S
              )
    
    ######## insert les barycenters
    print("---------------------- insert baryc")
    #### NS
    plot_baryc(metadata=  metadata, 
               data=      traj_NS, 
               list_frame=times, 
               mapping=   ax_mapping_baryc_NS, 
               xticks=    xticks_NS, 
               yticks=    yticks_NS,
               palette=   palette_for_time
               )
    
    #### S
    plot_baryc(metadata=  metadata, 
           data=      traj_S, 
           list_frame=times, 
           mapping=   ax_mapping_baryc_S, 
           xticks=    xticks_S, 
           yticks=    yticks_S,
           palette=   palette_for_time
           )
#

add_stat(ax1, eff_size={'ES':r'$A_w$', 'values':[0.23, 0.07, 0.08, 0.43]}, y=[.59]+3*[.955], sign=3*["**"]+["*"])                   # Cohen's d=[1.4, 2., 47., 10.]
add_stat(ax2, eff_size={'ES':r'$A_w$', 'values':[0.40, 0.83, 0.89, 0.51]}, sign=3*["**"]+["NS"])                                     # Cohen's d=[0.6, 1.2, 14., 5.2]
add_stat(ax3, eff_size={'ES':r'$A_w$', 'values':[0.51, 0.36, 0.31, 0.5]}, y=[.59]+3*[.955], sign=3*["**"]+["NS"])                   # Cohen's d=[0.04, 0.26, 3.1, 0.8]
add_stat(ax4, eff_size={'ES':r'$A_w$', 'values':[0.23, 0.07, 0.08, 0.43]}, sign=3*["**"]+["*"])                                     # Cohen's d=[1.4, 2., 47., 10.]


if save:
    print('\n------------------------------------ export figure:')
    plt.savefig(path_to_save)
#
