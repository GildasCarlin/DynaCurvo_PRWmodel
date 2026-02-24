# ================================================================================================================
# ================================================================================================================
#                                                    packages importation
# ================================================================================================================
# ================================================================================================================

import numpy as np
from sklearn.metrics import r2_score

import os
import regex as re
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
from matplotlib.gridspec import GridSpec
import seaborn as sns

# ================================================================================================================
# ================================================================================================================
#                                                    open figure
# ================================================================================================================
# ================================================================================================================

######## open figure

NbLignes, NbColonnes=1, 3
Largeur, Hauteur = 24, 8
fig1 = plt.figure(num='role_of_Vmax', figsize=(Largeur, Hauteur), constrained_layout=True)
gs1 = GridSpec(NbLignes, NbColonnes, figure=fig1)

# axes pour plot les autres indicateurs (Dx, velocity etc)
ax1 = fig1.add_subplot(gs1[0, 0])
ax2 = fig1.add_subplot(gs1[0, 1])
ax3 = fig1.add_subplot(gs1[0, 2])
#ax4 = fig1.add_subplot(gs1[1, 1])

#ax5 = fig1.add_subplot(gs1[1, 0])
#ax6 = fig1.add_subplot(gs1[1, 1])
#ax7 = fig1.add_subplot(gs1[1, 2])
#ax8 = fig1.add_subplot(gs1[1, 3])

## add labels
def add_label(ax, txt, x=-0.05, y=1.1):
    ax.text(x, y, txt, transform=ax.transAxes, fontsize=20, fontweight='bold', va='top', ha='right')
#


add_label(ax1, 'A')
add_label(ax2, 'B')
add_label(ax3, 'C')
#add_label(ax4, 'D')
#add_label(ax5, 'E')
#add_label(ax6, 'F')
#add_label(ax7, 'G')
#add_label(ax8, 'H')

## font size
plt.rcParams.update({'font.size': 50})


# ================================================================================================================
# ================================================================================================================
#                                                    intermediate functions
# ================================================================================================================
# ================================================================================================================


def plot_bary_by_c(data, ax, palette, plot_linear=False):

    sns.scatterplot(data=data,
                    x='celerity',
                    y='last_y',
                    hue='Vmax',
                    ax=ax,
                    palette=palette,
                    s=150
    )

    if plot_linear:
        x=np.linspace(0., 2.5, 20)
        ax.plot(x, x+1.5, color='grey')
        ax.text(x=.5, y=3., s=r'$y=x$', fontsize=13., rotation=71.)
    #
    
    ax.tick_params(axis='both', which='major', labelsize=20)
    ax.set_xlabel('Wave celerity (µm/s)', fontsize=25)
    ax.set_ylabel(r'Barycenter position at $T$=20min (mm)', fontsize=25)

    handles, labels = plt.gca().get_legend_handles_labels()
    mpatch1 = mpatches.Patch(facecolor=palette[0], label=r'2', edgecolor='black')
    mpatch2 = mpatches.Patch(facecolor=palette[1], label=r'4', edgecolor='black')
    mpatch3 = mpatches.Patch(facecolor=palette[2], label=r'6', edgecolor='black')
    handles.extend([mpatch1, mpatch2, mpatch3])
    ax.legend(handles=handles, fontsize=20, title=r"$V_{max}$ (µm/s)", title_fontsize=22)
    ax.tick_params(axis='both', which='major', labelsize=20)
#

def s_dyna(x, V):
    return x/((x/V)**3.+1.)

def plot_central(data, quantity, ax, palette, plot_legend=True, plot_sdyna=True):
    # specify options
    if quantity=='dx':
        infos={'y_label':r'$\overline{D}_{\perp}$'}

        ## dx=0.5 ==> no bias
        ax.axhline(y=.5, linewidth = 2, color='#8C0335', linestyle='--')
        
        ## dx belongs to [0, 1]
        ax.set_ylim(bottom=0., top=1.)

    if quantity=='vd':
        infos={'y_label':r'$\overline{v}$ (µm/s)'}
        if plot_sdyna:
            x=np.linspace(0., 20., 100)
            ax.plot(x, s_dyna(x, 2.), color=palette[0], linewidth=3., linestyle='--')
            ax.plot(x, s_dyna(x, 4.), color=palette[1], linewidth=3., linestyle='--')
            ax.plot(x, s_dyna(x, 6.), color=palette[2], linewidth=3., linestyle='--')
        #
    #

    if quantity=='r_d':
        infos={'y_label':r'$\overline{r}_d$'}
        ax.axhline(y=0., linewidth = 2, color='#8C0335', linestyle='--')
    #

    sns.scatterplot(data=data,
                    x='celerity',
                    y=quantity,
                    hue='Vmax',
                    ax=ax,
                    palette=palette,
                    legend=plot_legend,
                    s=150
                    )
    #
    
    ax.tick_params(axis='both', which='major', labelsize=20)
    ax.set_xlabel('Wave celerity (µm/s)', fontsize=25)
    ax.set_ylabel(infos['y_label'], fontsize=25)
 
    if plot_legend:
        handles, labels = plt.gca().get_legend_handles_labels()
        mpatch1 = mpatches.Patch(facecolor=palette[0], label=r'2', edgecolor='black')
        mpatch2 = mpatches.Patch(facecolor=palette[1], label=r'4', edgecolor='black')
        mpatch3 = mpatches.Patch(facecolor=palette[2], label=r'6', edgecolor='black')
        point1 = Line2D([0], [0], label=r'$s_{dyna}$', markersize=10, linestyle='--', color='black', linewidth=3.)
        handles.extend([mpatch1, mpatch2, mpatch3, point1])
        ax.legend(handles=handles, fontsize=23, title=r"$V_{max}$ (µm/s)", title_fontsize=25)
        ax.tick_params(axis='both', which='major', labelsize=20)
    #
#

def plot_critical(data, y, r2, ax, palette, plot_legend=False, xticks = [2.0, 4.0, 6.0], alpha=.5, th=(False, .79, 'th')):
    marker={'last_y':'^', 'dx':'o', 'vd':'s'}
    ylabel={'celerity':"Critical wave celerity (µm/s)", 'vd':r"Maximum $\overline{v}$ (µm/s)"}
    
    ## plot theoretical expression
    if th[0]:
        x=np.linspace(2., 6., 3)
        ax.plot(x, th[1]*x, color='black', zorder=1, linestyle='--', linewidth=5.)
    #
    
    ## plot linear interpolation
    sns.lineplot(data=data,
                    x='Vmax',
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
                    x='Vmax',
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


    ## write R^2 values
    for i, val in enumerate(r2):
        txt = f"$R^2$ = {val:.2f}"  
        ax.text(.7, .2 - i * 0.05, txt, transform=ax.transAxes, fontsize=21, va='top', ha='left', color=palette[i])
    
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

## q or Vmax
parameter = 'Vmax'

values = [4.0, 5.0, 6.0]
assumptions = ['cell_is_surfer']

## metadada
wavelengths_explicit = [20.]
wavelengths_critic = [20., 40., 80.]
max_celerity = 20.
dt = 15.

## compute r2 score
perform_stat = False
dependant_variable='vd'

# export figure or not
save = True
path_to_save = 'role_of_' + parameter + '\\role_' + parameter + '.svg'


## platte de couleur 
wl_palette = ['#601A4A',  # wavelength=20 µm
              '#EE442F',  # 40 µm
              '#63ACBE'   # 80 µm
              ] 

Vmax_palette = ["#e5d4e8",
                "#601A4A", 
                "#986eac"]


# ================================================================================================================
# ================================================================================================================
#                                                    MAIN function
# ================================================================================================================
# ================================================================================================================

sources=[]
for value in values:
    for assumption in assumptions:
        ## particular cases
        if (parameter=='Vmax' and value==4.0) or (parameter=='q' and value==1.0):
            sources.append(assumption)
        else:
            path='role_of_' + parameter + '\\' + parameter + '=' + str(value) + '\\' + assumption
            sources.append(path)
        #
    #
#

root = os.getcwd()
list_last, list_dx, list_vd = [], [], []
list_c_last, list_c_dx, list_c_vd = [], [], []

for source in sources:

    print('\n------------------------------------ extract Vmax value')
    test = re.findall(r"Vmax=(\d*.\d*)", source)
    if len(test) == 0:
        Vmax = 4.0
    else:
        Vmax = float(test[0])
    #
    print("... Vmax =", Vmax, "µm/s")
    
    print('------------------------------------ open data')

    path_to_last = root + '\\' + source + '\\analysis_results\\global\\central_last.csv'
    path_to_dx =   root + '\\' + source + '\\analysis_results\\global\\central_dx.csv'
    path_to_vd =   root + '\\' + source + '\\analysis_results\\global\\central_vd.csv'

    data_last = pd.read_csv(path_to_last)
    data_dx =   pd.read_csv(path_to_dx)
    data_vd =   pd.read_csv(path_to_vd)

    print('------------------------------------ extract max/min')
    critic_c_last = data_last.loc[data_last.groupby("wavelength")["last_y"].idxmax(), ["wavelength", "celerity", "last_y"]].reset_index(drop=True)
    critic_c_dx =   data_dx.loc[data_dx.groupby("wavelength")["dx"].idxmin(), ["wavelength", "celerity", "dx"]].reset_index(drop=True)
    critic_c_vd =   data_vd.loc[data_vd.groupby("wavelength")["vd"].idxmax(), ["wavelength", "celerity", "vd"]].reset_index(drop=True)

    critic_c_last.insert(0, 'Vmax', len(critic_c_last)*[Vmax])
    critic_c_dx.insert(0, 'Vmax', len(critic_c_dx)*[Vmax])
    critic_c_vd.insert(0, 'Vmax', len(critic_c_vd)*[Vmax])

    critic_c_last = critic_c_last[critic_c_last['wavelength'].isin(wavelengths_critic)]
    critic_c_dx = critic_c_dx[critic_c_dx['wavelength'].isin(wavelengths_critic)]
    critic_c_vd = critic_c_vd[critic_c_vd['wavelength'].isin(wavelengths_critic)]

    list_c_last.append(critic_c_last)
    list_c_dx.append(critic_c_dx)
    list_c_vd.append(critic_c_vd)


    print('------------------------------------ sort data')
    
    data_last = data_last[data_last['wavelength'].isin(wavelengths_explicit)]
    data_dx =   data_dx[data_dx['wavelength'].isin(wavelengths_explicit)]
    data_vd =   data_vd[data_vd['wavelength'].isin(wavelengths_explicit)]

    data_last = data_last[data_last['celerity'] <= max_celerity]
    data_dx =   data_dx[data_dx['celerity'] <= max_celerity]
    data_vd =   data_vd[data_vd['celerity'] <= max_celerity]

    data_last.insert(0, 'Vmax', len(data_last)*[Vmax])
    data_dx.insert(0, 'Vmax', len(data_dx)*[Vmax])
    data_vd.insert(0, 'Vmax', len(data_vd)*[Vmax])
    
    list_last.append(data_last)
    list_dx.append(data_dx)
    list_vd.append(data_vd)
#

print('\n------------------------------------ concat data')
df_last = pd.concat(list_last, ignore_index=True)
df_last['last_y'] /= 1000.  ## µm -->  mm
df_dx = pd.concat(list_dx, ignore_index=True)
df_vd = pd.concat(list_vd, ignore_index=True)

df_c_last = pd.concat(list_c_last, ignore_index=True)
df_c_dx =   pd.concat(list_c_dx, ignore_index=True)
df_c_vd =   pd.concat(list_c_vd, ignore_index=True)

df_c_last = df_c_last.sort_values(by='Vmax', ascending=True).reset_index(drop=True)
df_c_dx = df_c_dx.sort_values(by='Vmax', ascending=True).reset_index(drop=True)
df_c_vd = df_c_vd.sort_values(by='Vmax', ascending=True).reset_index(drop=True)

if perform_stat:
    print('\n------------------------------------ compute r2 and rmse')

    if dependant_variable in ['last_y', 'r_d']:
        data=df_c_last.copy()
    elif dependant_variable == 'dx':
        data=df_c_dx.copy()
    elif dependant_variable == 'vd':
        data=df_c_vd.copy()
    #

    ## theoretical tables 
    x=np.linspace(2., 6., 3)
    y_true_c, y_true_a=2.**(-1./3.)*x, 2.**(2./3.)/3.*x
    
    r2_c, r2_a = [], []
    for wl in wavelengths_critic:
        print('\n... wavelength =', wl)
        df = data[data['wavelength']==wl]
        r2_c.append(r2_score(y_true_c, df['celerity']))
        print('...r2 for critic celerity', r2_c[-1])
        
        if dependant_variable=='vd':
            r2_a.append(r2_score(y_true_a, df[dependant_variable]))
            print('...r2 for average velocity at critic celerity', r2_a[-1])
        #
    #
#

if save:
    print('\n------------------------------------ plot')
    print("... central curves")
    #plot_bary_by_c(df_last, ax=ax1, palette=Vmax_palette, plot_linear=True)
    #plot_central(df_dx, 'dx', ax = ax2, palette=Vmax_palette, plot_legend=False)
    plot_central(df_vd, 'vd', ax = ax1, palette=Vmax_palette, plot_legend=True)

    print("... critical celerity")
    plot_critical(df_c_vd, y='celerity', r2=r2_c, ax=ax2, palette=wl_palette, plot_legend=True, alpha=1., th=(True, 2.**(-1./3.), r'$w_{sat}$'))
    #plot_critical_c(df_c_dx, ax=ax2, palette=wl_palette)
    #plot_critical_c(df_c_vd, ax=ax2, palette=wl_palette)

    if dependant_variable=='vd':
        print("... Avelocity average at critical celerity")
        plot_critical(df_c_vd, y='vd', r2=r2_a, ax=ax3, palette=wl_palette, plot_legend=True, alpha=1., th=(True, 2.**(2./3.)/3., r'$\max_c (s_{dyna})$'))

    print('------------------------------------ export figure:')
    plt.savefig(root + '\\' +  path_to_save)
