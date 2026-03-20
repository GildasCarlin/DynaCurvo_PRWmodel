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
from matplotlib.lines import Line2D
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

### axes d'insert 
#ax_ins1 = inset_axes(ax1, width=Largeur_ins, height=Hauteur_ins)
#ax_ins2 = ax2.inset_axes([.017, 0.06, .4, .4])
#ax_ins3 = inset_axes(ax3, width=Largeur_ins, height=Hauteur_ins)
#ax_ins4 = ax4.inset_axes([.59, 0.057, .4, .4])
#
#ax_ins5 = inset_axes(ax5, width=Largeur_ins, height=Hauteur_ins)
#ax_ins6 = ax6.inset_axes([0.08, 0.59, .4, .4])
#ax_ins7 = inset_axes(ax7, width=Largeur_ins, height=Hauteur_ins)
#ax_ins8 = ax8.inset_axes([.59, 0.57, .4, .4])

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
#                                                    User interface
# ================================================================================================================
# ================================================================================================================

# study
substrate_regime="dyna"  # dyna or static
NDSC_study="NDSC_CurvatureCorrection_RightSinus_20min_VariousCelerity"
DSC_study="DSC_CurvatureCorrection_RightSinus_20min_VariousCelerity"
csat_study="DSC_lookingForCsat_fstat=0_q=0"
csat_nonoise_study="DSC_lookingForCsat_fstat=0_q=0_noise=0"
csat_taudyna_study="DSC_looking_for_csat_tau_dyna_div_par10"


## metadada
assumptions = ['NDSC', 'DSC']
dt = 15.

# mapping between assumption and quantity and axes
mapping={('NDSC', 'b'):ax1,
         ('NDSC', 'dx'):ax2,
         ('NDSC', 'vd'):ax3,
         ('NDSC', 'rd'):ax3,
         ('DSC', 'b'):ax5,
         ('DSC', 'dx'):ax6,
         ('DSC', 'vd'):ax7,
         ('DSC', 'rd'):ax7
         }

# plot sdyna theorical expression
#plot_or_not={'NDSC':False, 'DSC':True}

# export figure or not
save = True




# ================================================================================================================
# ================================================================================================================
#                                                    intermediate functions
# ================================================================================================================
# ================================================================================================================

def formatData(df_last, df_dx, df_vd, wavelengths=[20., 40., 80.], max_celerity={"NDSC":2.5, "DSC":2.5}):
    
    # format columns names
    df_last.rename(columns={"last_y":"last_y.mean",
                              "last_y.1":"last_y.med",
                              "last_y.2":"last_y.std",
                              "r_d":"rd.mean",
                              "r_d.1":"rd.med",
                              "r_d.2":"rd.std"}, inplace=True)
    df_last.dropna(inplace=True)
    df_last[['last_y.mean','last_y.med', 'last_y.std', 'rd.mean', 'rd.med', 'rd.std']]=df_last[['last_y.mean', 'last_y.med', 'last_y.std', 'rd.mean', 'rd.med', 'rd.std']].astype("float")
    
    df_dx.rename(columns={"dx":"dx.mean",
                            "dx.1":"dx.med",
                            "dx.2":"dx.std"}, inplace=True)
    df_dx.dropna(inplace=True)
    df_dx[['dx.mean', 'dx.med', 'dx.std']]=df_dx[['dx.mean', 'dx.med', 'dx.std']].astype("float")

    df_vd.rename(columns={"speed":"vd.mean",
                            "speed.1":"vd.med",
                            "speed.2":"vd.std"}, inplace=True)
    df_vd.dropna(inplace=True)
    df_vd[['vd.mean', 'vd.med', 'vd.std']]=df_vd[['vd.mean', 'vd.med', 'vd.std']].astype("float")

    ## isolate only the specified wavelength
    print('... sort by wavelength and celerity')
    df_last = df_last[df_last['wavelength'].isin(wavelengths)]
    df_dx =   df_dx[df_dx['wavelength'].isin(wavelengths)]
    df_vd =   df_vd[df_vd['wavelength'].isin(wavelengths)]

    ## isolate only the specified wave celerity
    df_last = df_last[df_last['celerity'] <= max_celerity[assumption]]
    df_dx =   df_dx[df_dx['celerity'] <=     max_celerity[assumption]]
    df_vd =   df_vd[df_vd['celerity'] <=     max_celerity[assumption]]

    # µm -->  mm
    df_last[['last_y.mean', 'last_y.med', 'last_y.std']] /= 1000.

    return df_last, df_dx, df_vd
#

def searchForCsat(df_last, df_dx, df_vd, wavelengths=[20., 40., 80.], Vmax=[.61, .56, .5]):
    print('------------------------------------ extract max/min')
    critic_c_last = df_last.loc[df_last.groupby("wavelength")["last_y.mean"].idxmax(), ["wavelength", "celerity", "last_y.mean"]].reset_index(drop=True)
    critic_c_dx =   df_dx.loc[df_dx.groupby("wavelength")["dx.mean"].idxmin(), ["wavelength", "celerity", "dx.mean"]].reset_index(drop=True)
    critic_c_vd =   df_vd.loc[df_vd.groupby("wavelength")["vd.mean"].idxmax(), ["wavelength", "celerity", "vd.mean"]].reset_index(drop=True)

    critic_c_last.insert(0, 'Vmax', Vmax)
    critic_c_dx.insert(0, 'Vmax', Vmax)
    critic_c_vd.insert(0, 'Vmax', Vmax)
    critic_c_last.insert(0, 'assumption', len(critic_c_last)*[assumption])
    critic_c_dx.insert(0, 'assumption',len(critic_c_dx)*[assumption])
    critic_c_vd.insert(0, 'assumption', len(critic_c_vd)*[assumption])

    critic_c_last = critic_c_last[critic_c_last['wavelength'].isin(wavelengths)]
    critic_c_dx = critic_c_dx[critic_c_dx['wavelength'].isin(wavelengths)]
    critic_c_vd = critic_c_vd[critic_c_vd['wavelength'].isin(wavelengths)]

    return critic_c_last, critic_c_dx, critic_c_vd
#

def plotCentralStat(ass,
                    data,
                    y,          # name of y column to plot
                    quantity,   # specifications (ax label, etc)
                    ax,         # axis to plot 
                    palette=['#601A4A', '#EE442F', '#63ACBE'],      # color palette for wavelength
                    plot_label=True,       
                    plot_legend=True, 
                    plot_sdyna=False, 
                    x_sdyna=4., 
                    n=3., 
                    zorder=10,
                    remove_ylabel=False):
    
    ##### specific option to plot depending on the statistic
    if quantity=='dx':
        infos={'y_label':r'$\overline{D}_{\perp}$'}
        ax.axhline(y=.5, linewidth = 2, color='#8C0335', linestyle='--', zorder=1)
        ax.set_ylim(bottom=0., top=1.)
    #
    if quantity=='vd':
        infos={'y_label':r'$\overline{v}$ (µm/s)'}
        ax.set_ylim(bottom=0.)
    #
    if quantity=='rd':
        infos={'y_label':r'$\overline{r}_d$'}
        ax.axhline(y=0., linewidth = 2, color='#8C0335', linestyle='--', zorder=1)
        ax.set_ylim(top=1., bottom=-0.15)
    #
    if quantity=="b":
        infos={'y_label':r"$B_T$ with $T$=20min (mm)"}
        ax.set_ylim(bottom=-.1, top=.9)
        ax.text(-.25, .5, 
                ass+r" cells", 
                rotation = 'vertical', fontsize=30, fontweight='bold', va='center', ha='center', transform=ax.transAxes
                )
    #


    # theoretical expression of the target dynamic component of speed (WiP)
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

    ##### plot of the central statistic VS the wave celerity
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

    if plot_legend:
        handles, labels = plt.gca().get_legend_handles_labels()
        mpatch1 = mpatches.Patch(facecolor=palette[0], label=r'20', edgecolor='black')
        mpatch2 = mpatches.Patch(facecolor=palette[1], label=r'40', edgecolor='black')
        mpatch3 = mpatches.Patch(facecolor=palette[2], label=r'80', edgecolor='black')
        handles.extend([mpatch1, mpatch2, mpatch3])
        ax.legend(handles=handles, fontsize=20, title=r"Wavelength (µm)", title_fontsize=22)
        ax.tick_params(axis='both', which='major', labelsize=20)
    #
#

def plotCritical(data,
                 y,
                 ax,
                 marker="o",
                 palette=['#601A4A', '#EE442F', '#63ACBE'],
                 plot_legend=False,
                 xticks = [0.5, 0.56, 0.61],
                 alpha=.5,
                 th=(False, .79, 'th')):
    
    #marker={'last_y.mean':'^', 'dx.mean':'o', 'vd.mean':'s'}
    ylabel={"celerity":"Critical wave celerity (µm/s)", 'vd.mean':r"Maximum $\overline{v}$ (µm/s)", 'last_y.mean':r"Maximum $\overline{v}$ (µm/s)"}
    
    ## plot theoretical expression
    if th[0]:
        x=np.linspace(.45, .65, 50)
        ax.plot(x, th[1]*x, color='black', zorder=1, linestyle='--', linewidth=5.)
    #
    
    ## plot linear interpolation
    #sns.lineplot(data=data,
    #                x='Vmax',
    #                y=y,
    #                hue='wavelength',
    #                ax=ax,
    #                palette=palette,
    #                linewidth=3.,
    #                legend=plot_legend,
    #                zorder=2
    #                )
    ##

    ## plot points
    sns.scatterplot(data=data,
                    x='Vmax',
                    y=y,
                    hue='wavelength',
                    ax=ax,
                    palette=palette,
                    edgecolor="black",
                    linewidth=1,
                    marker=marker,
                    s=300,
                    legend=plot_legend,
                    zorder=3
                    )
    #

    ### reduce alpha of points
    #for coll in ax.collections:
    #    facecolors = coll.get_facecolors()
    #    facecolors[:, -1] = alpha  # % de transparence
    #    coll.set_facecolors(facecolors)


    ### write R^2 values
    #for i, val in enumerate(r2):
    #    txt = f"$R^2$ = {val:.2f}"  
    #    ax.text(.7, .2 - i * 0.05, txt, transform=ax.transAxes, fontsize=21, va='top', ha='left', color=palette[i])
    #
    ## axis
    ax.tick_params(axis='both', which='major', labelsize=20)
    ax.set_xlabel("$V_{max}$ (µm/s)", fontsize=25)
    ax.set_ylabel(ylabel[y], fontsize=25)
    ax.set_xticks(xticks)

    ## legend
    if plot_legend:
        handles, labels = plt.gca().get_legend_handles_labels()
        point1 = Line2D([0], [0], label=th[2], markersize=10, linestyle='--', color='black', linewidth=3.)
        point2 = Line2D([0], [0], label="DSCcells", markersize=10, marker='s', color='black', ms=15, linestyle="")
        point3 = Line2D([0], [0], label="DSCcells,\n"+r"$\bf f_{stat}=0$", markersize=10, marker='o', color='black', ms=15, linestyle="")
        #point4 = Line2D([0], [0], label="DSCcells,\n"+r"$\bf f_{stat}=0$, $\bf \sigma=0$", markersize=10, marker='d', color='black', ms=15, linestyle="")
        ax.legend(handles=[point1, point2, point3], fontsize=20, loc="upper left")
        ax.tick_params(axis='both', which='major', labelsize=20)
    #
#




# ================================================================================================================
# ================================================================================================================
#                                                    MAIN function
# ================================================================================================================
# ================================================================================================================

# assets
root = os.getcwd()
NDSC_assets=os.path.join(root, substrate_regime, "assets", NDSC_study)
DSC_assets=os.path.join(root, substrate_regime, "assets", DSC_study)
csat_assets=os.path.join(root, substrate_regime, "assets", csat_study)
csat_nonoise_assets=os.path.join(root, substrate_regime, "assets", csat_nonoise_study)
csat_taudyna_assets=os.path.join(root, substrate_regime, "assets", csat_taudyna_study)
mapping_path_assumption={"NDSC":NDSC_assets, "DSC":DSC_assets}

# path to save figure
path_to_save=os.path.join(root, substrate_regime, "wavelength_effect.svg")


list_c_last, list_c_dx, list_c_vd = [], [], []
for assumption in assumptions:
    print('\n----------------------------------- treatment of ' + assumption)
    print('... open data')
    data_last = pd.read_csv(mapping_path_assumption[assumption] + '\\' + 'central_last.csv')
    data_dx =   pd.read_csv(mapping_path_assumption[assumption] + '\\' + 'central_dx.csv')
    data_vd =   pd.read_csv(mapping_path_assumption[assumption] + '\\' + 'central_vd.csv')

    data_last, data_dx, data_vd = formatData(data_last, data_dx, data_vd)

    ### extract the saturation wave celerity (only for DSC cells)
    if assumption=="DSC":
        csat=searchForCsat(data_last, data_dx, data_vd)
        
        list_c_last.append(csat[0])
        list_c_dx.append(csat[1])
        list_c_vd.append(csat[2])
    #

    print('... plot central statistics')
    quantity="b"
    plotCentralStat(assumption,
                    data_last, y="last_y.mean", 
                    quantity=quantity,
                    ax=mapping[(assumption, quantity)],
                    plot_legend=True,
                    )
    
    quantity="dx"
    plotCentralStat(assumption,
                    data_dx, y="dx.mean", 
                 quantity=quantity,
                 ax=mapping[(assumption, quantity)],
                 plot_legend=True,
                 )
    
    quantity="rd"
    plotCentralStat(assumption,data_last, y="rd.mean", 
                 quantity=quantity,
                 ax=mapping[(assumption, quantity)], 
                 plot_legend=True
                 )
#

## plot critical wave celerity, f_stat!=0 and q=1 
df_c_last = pd.concat(list_c_last, ignore_index=True)
df_c_dx =   pd.concat(list_c_dx, ignore_index=True)
df_c_vd =   pd.concat(list_c_vd, ignore_index=True)
df_c_last = df_c_last.sort_values(by='Vmax', ascending=True).reset_index(drop=True)
df_c_dx = df_c_dx.sort_values(by='Vmax', ascending=True).reset_index(drop=True)
df_c_vd = df_c_vd.sort_values(by='Vmax', ascending=True).reset_index(drop=True)
print(df_c_last)
plotCritical(df_c_last, y='celerity', marker="s", ax=ax4, plot_legend=True, alpha=1., th=(True, 2.**(-1./3.), r'$w_{sat}$'))
plotCritical(df_c_last, y='last_y.mean', marker="s", ax=ax8, plot_legend=True, alpha=1., th=(True, 2.**(2./3.)/3., r'$\max_c (s_{dyna})$'))

### plot critical wave celerity, f_stat=0, tau_dyna/10, q=0
data_last = pd.read_csv(csat_assets + '\\' + 'central_last.csv')
data_dx =   pd.read_csv(csat_assets + '\\' + 'central_dx.csv')
data_vd =   pd.read_csv(csat_assets + '\\' + 'central_vd.csv')

data_last, data_dx, data_vd = formatData(data_last, data_dx, data_vd)

csat=searchForCsat(data_last, data_dx, data_vd)
        
list_c_last=[csat[0]]
list_c_dx=[csat[1]]
list_c_vd=[csat[2]]

df_c_last = pd.concat(list_c_last, ignore_index=True)
df_c_dx =   pd.concat(list_c_dx, ignore_index=True)
df_c_vd =   pd.concat(list_c_vd, ignore_index=True)
df_c_last = df_c_last.sort_values(by='Vmax', ascending=True).reset_index(drop=True)
df_c_dx = df_c_dx.sort_values(by='Vmax', ascending=True).reset_index(drop=True)
df_c_vd = df_c_vd.sort_values(by='Vmax', ascending=True).reset_index(drop=True)
plotCritical(df_c_last, y='celerity', ax=ax4, plot_legend=False, alpha=1., th=(False, 0., ""))
plotCritical(df_c_last, y='last_y.mean', ax=ax8, plot_legend=False, alpha=1., th=(False, 0., ""))

"""
# ce n'est pas le noise qui explique le décalage
### plot critical wave celerity, f_stat=0, q=0 and tau_dyna/10 and no noise
data_last = pd.read_csv(csat_nonoise_assets + '\\' + 'central_last.csv')
data_dx =   pd.read_csv(csat_nonoise_assets + '\\' + 'central_dx.csv')
data_vd =   pd.read_csv(csat_nonoise_assets + '\\' + 'central_vd.csv')

data_last, data_dx, data_vd = formatData(data_last, data_dx, data_vd)

csat=searchForCsat(data_last, data_dx, data_vd)
        
list_c_last=[csat[0]]
list_c_dx=[csat[1]]
list_c_vd=[csat[2]]

df_c_last = pd.concat(list_c_last, ignore_index=True)
df_c_dx =   pd.concat(list_c_dx, ignore_index=True)
df_c_vd =   pd.concat(list_c_vd, ignore_index=True)
df_c_last = df_c_last.sort_values(by='Vmax', ascending=True).reset_index(drop=True)
df_c_dx = df_c_dx.sort_values(by='Vmax', ascending=True).reset_index(drop=True)
df_c_vd = df_c_vd.sort_values(by='Vmax', ascending=True).reset_index(drop=True)
plotCritical(df_c_last, y='celerity', marker="d", ax=ax4, plot_legend=False, alpha=1., th=(False, 0., ""))
plotCritical(df_c_last, y='last_y.mean', marker="d", ax=ax8, plot_legend=False, alpha=1., th=(False, 0., ""))
"""

############# plot critical wave celerity when curvature assumption is deactivated and q=0 and tau_dyna / 2
####data_last = pd.read_csv(csat_taudyna_assets + '\\' + 'central_last.csv')
####data_dx =   pd.read_csv(csat_taudyna_assets + '\\' + 'central_dx.csv')
####data_vd =   pd.read_csv(csat_taudyna_assets + '\\' + 'central_vd.csv')
####
####data_last, data_dx, data_vd = formatData(data_last, data_dx, data_vd)
####
####csat=searchForCsat(data_last, data_dx, data_vd)
####        
####list_c_last=[csat[0]]
####list_c_dx=[csat[1]]
####list_c_vd=[csat[2]]
####
####df_c_last = pd.concat(list_c_last, ignore_index=True)
####df_c_dx =   pd.concat(list_c_dx, ignore_index=True)
####df_c_vd =   pd.concat(list_c_vd, ignore_index=True)
####df_c_last = df_c_last.sort_values(by='Vmax', ascending=True).reset_index(drop=True)
####df_c_dx = df_c_dx.sort_values(by='Vmax', ascending=True).reset_index(drop=True)
####df_c_vd = df_c_vd.sort_values(by='Vmax', ascending=True).reset_index(drop=True)
####plotCritical(df_c_vd, y='celerity', marker="d", ax=ax3, plot_legend=True, alpha=1., th=(False, 0., ""))
####plotCritical(df_c_vd, y='vd.mean', marker="d", ax=ax7, plot_legend=True, alpha=1., th=(False, 0., ""))

if save:
    print('------------------------------------ export figure:')
    plt.savefig(path_to_save)
#