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

# c_sat DSC(+) cells
c_sat="DSC_SpeedDistribution"


#csat_study="DSC_lookingForCsat_q=0"
#csat_nonoise_study="DSC_lookingForCsat_fstat=0_q=0_noise=0"
#csat_taudyna_study="DSC_looking_for_csat_tau_dyna_div_par10"


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

def formatData(df_last, df_dx, wavelengths=[20., 40., 80.], max_celerity={"NDSC":2.5, "DSC":2.5}):
    
    # format columns names
    df_last.rename(columns={"last_y":"last_y.mean",
                              "last_y.1":"last_y.med",
                              "last_y.2":"last_y.std",
                              "r_d":"rd.mean",
                              "r_d.1":"rd.med",
                              "r_d.2":"rd.std"}, inplace=True)
    df_last.drop(0, axis=0, inplace=True)
    df_last[['last_y.mean','last_y.med', 'last_y.std', 'rd.mean', 'rd.med', 'rd.std']]=df_last[['last_y.mean', 'last_y.med', 'last_y.std', 'rd.mean', 'rd.med', 'rd.std']].astype("float")
    
    df_dx.rename(columns={"dx":"dx.mean",
                        "dx.1":"dx.med",
                        "dx.2":"dx.std"}, inplace=True)
    df_dx.dropna(inplace=True)
    df_dx[['dx.mean', 'dx.med', 'dx.std']]=df_dx[['dx.mean', 'dx.med', 'dx.std']].astype("float")

    #df_vd.rename(columns={"speedCell":"vdCell.mean",
    #                      "speedCell.1":"vdCell.med",
    #                      "speedCell.2":"vdCell.std"}, inplace=True)
    #df_vd.dropna(inplace=True)
    #df_vd[['vdCell.mean', 'vdCell.med', 'vdCell.std']]=df_vd[['vdCell.mean', 'vdCell.med', 'vdCell.std']].astype("float")

    ## isolate only the specified wavelength
    df_last = df_last[df_last['wavelength'].isin(wavelengths)]
    df_dx =   df_dx[df_dx['wavelength'].isin(wavelengths)]
    #df_vd =   df_vd[df_vd['wavelength'].isin(wavelengths)]

    ## isolate only the specified wave celerity
    df_last = df_last[df_last['celerity'] <= max_celerity[assumption]]
    df_dx =   df_dx[df_dx['celerity'] <=     max_celerity[assumption]]
    #df_vd =   df_vd[df_vd['celerity'] <=     max_celerity[assumption]]

    # µm -->  mm
    df_last[['last_y.mean', 'last_y.med', 'last_y.std']] /= 1000.

    return df_last, df_dx
#

def formatData1(df_last, wavelengths=[20., 40., 80.], max_celerity={"NDSC":2.5, "DSC":2.5}):
    
    # format columns names
    df_last.rename(columns={"last_y":"last_y.mean",
                              "last_y.1":"last_y.med",
                              "last_y.2":"last_y.std",
                              "r_d":"rd.mean",
                              "r_d.1":"rd.med",
                              "r_d.2":"rd.std"}, inplace=True)
    df_last.drop(0, axis=0, inplace=True)
    df_last[['last_y.mean','last_y.med', 'last_y.std', 'rd.mean', 'rd.med', 'rd.std']]=df_last[['last_y.mean', 'last_y.med', 'last_y.std', 'rd.mean', 'rd.med', 'rd.std']].astype("float")

    ## isolate only the specified wavelength
    df_last = df_last[df_last['wavelength'].isin(wavelengths)]

    ## isolate only the specified wave celerity
    df_last = df_last[df_last['celerity'] <= max_celerity[assumption]]

    # µm -->  mm
    df_last[['last_y.mean', 'last_y.med', 'last_y.std']] /= 1000.

    return df_last
#

def formatDataSpeed(df):
    df=df.dropna()
    df=df.drop(["speedObs", "speedObs.1", "speedObs.2", "speedCell.1", "VzSub", "VzSub.1", "VzSub.2"], axis=1)
    df[["speedCell", "speedCell.2"]]=df[["speedCell", "speedCell.2"]].astype("float")
    df=df.reset_index(drop=True)
    return df
#


def searchForCsat(df_last, wavelengths=[20., 40., 80.], Vmax=[.61, .56, .5]):
    critic_c = df_last.loc[df_last.groupby("wavelength")["last_y.mean"].idxmax(), ["wavelength", "celerity", "last_y.mean"]].reset_index(drop=True)

    critic_c.insert(0, 'Vmax', Vmax)
    critic_c.insert(0, 'assumption', len(critic_c)*[assumption])
    critic_c = critic_c[critic_c['wavelength'].isin(wavelengths)]

    critic_c=pd.concat([critic_c], ignore_index=True)
    critic_c = critic_c.sort_values(by='Vmax', ascending=True).reset_index(drop=True)

    return critic_c
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
    
    name_mapping={"NDSC":"DSC(-)", "DSC":"DSC(+)"}
    
    ##### specific option to plot depending on the statistic
    if quantity=='dx':
        infos={'y_label':r'$\overline{D}_{\mathbf{k}}$'}
        ax.axhline(y=.5, linewidth = 2, color='#8C0335', linestyle='--', zorder=1)
        ax.set_ylim(bottom=0., top=1.)
    #
    if quantity=='vd':
        infos={'y_label':r'$\overline{v}$ (µm/s)'}
        ax.set_ylim(bottom=0.)
    #
    if quantity=='rd':
        infos={'y_label':r'$\overline{f}$'}
        ax.axhline(y=0., linewidth = 2, color='#8C0335', linestyle='--', zorder=1)
        ax.set_ylim(top=.9, bottom=-0.15)
    #
    if quantity=="b":
        infos={'y_label':r"$B_T$ with $T$=20min (mm)"}
        ax.set_ylim(bottom=-.05, top=.55)
        ax.text(-.25, .5, 
                name_mapping[ass]+r" cells", 
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
                    edgecolor=None,
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
                 th=(False, .79, 'th'),
                 sd=True):
    
    ## plot theoretical expression
    if th[0]:
        x=np.linspace(.45, .65, 50)
        ax.plot(x, th[1]*x, color='black', zorder=1, linestyle='--', linewidth=5.)
    #
    
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

    if y=="speedCell" and sd:
        #plot SD on cell speed
        ax.plot((data["Vmax"], data["Vmax"]), (data["speedCell"]+data["speedCell.2"],data["speedCell"]-data["speedCell.2"]),'-',color='black')

    ylabel={"celerity":r"$c_{sat}$ (µm/s)", 'speedCell':r"Average cell speed at $c_{sat}$ (µm/s)"}
    ## axis
    ax.tick_params(axis='both', which='major', labelsize=20)
    ax.set_xlabel(r"$\overline{v}_{\lambda}$ (µm/s)", fontsize=25)
    ax.set_ylabel(ylabel[y], fontsize=25)
    ax.set_xticks(xticks)

    ## legend
    if plot_legend:
        handles, labels = plt.gca().get_legend_handles_labels()
        point1 = Line2D([0], [0], label=th[2], markersize=10, linestyle='--', color='black', linewidth=3.)
        point2 = Line2D([0], [0], label="DSC(+) cells", markersize=10, marker='s', color='black', ms=15, linestyle="")
        point3 = Line2D([0], [0], label="DSC(+) cells\n"+r"and $\bf f_{pass}=0$", markersize=10, marker='o', color='black', ms=15, linestyle="")
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
#csat_assets=os.path.join(root, substrate_regime, "assets", csat_study)
#csat_nonoise_assets=os.path.join(root, substrate_regime, "assets", csat_nonoise_study)
#csat_taudyna_assets=os.path.join(root, substrate_regime, "assets", csat_taudyna_study)
mapping_path_assumption={"NDSC":NDSC_assets, "DSC":DSC_assets}
path_to_assets=os.path.join(root, substrate_regime, "assets")

print("Plot Figure 3\n---------------------------------------\n")

# path to save figure
path_to_save=os.path.join(root, substrate_regime, "wavelength_effect.svg")

print("\nSTEP 1:  plot central statistics (BT, f, Dk)")

for assumption in assumptions:
    print('\n----------------------------------- treatment of ' + assumption)
    print('... open data')
    data_last = pd.read_csv(os.path.join(mapping_path_assumption[assumption], 'central_last.csv'))
    data_dx   =   pd.read_csv(os.path.join(mapping_path_assumption[assumption], 'central_dx.csv'))

    data_last, data_dx = formatData(data_last, data_dx)

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




print("\n\nSTEP 2:  plot critical celerity c_sat")

print('\n... DSC(+) cells,  get c_sat')
# on cherche c_sat à partir du max(BT) par rapport à c
data=pd.read_csv(os.path.join(path_to_assets, DSC_study, "central_last.csv"))
data=formatData1(data)
c_sat=searchForCsat(data)

print('... DSC(+) cells, get cell speed at c_sat')
data_Speed=pd.read_csv(os.path.join(path_to_assets, "DSC_SpeedDistribution", "central_vd.csv"))
data_Speed=formatDataSpeed(data_Speed)
c_sat = c_sat.merge(data_Speed, on=['celerity', 'wavelength'], how='left', indicator=True)

print('... DSC(+) cells, plot c_sat and speed at c_sat')
plotCritical(c_sat, y='celerity', marker="s", ax=ax4, plot_legend=True, alpha=1., th=(True, 2.**(-1./3.), r'$w_{sat}$'))
plotCritical(c_sat, y='speedCell', marker="s", ax=ax8, plot_legend=True, alpha=1., th=(True, 2.**(2./3.)/3., r'$\max_c (s_{dyna})$'), sd=False)



###########
print('\n... DSC(+) cells w fpass=0,  get c_sat')
# on cherche c_sat à partir du max(BT) par rapport à c
dataq0=pd.read_csv(os.path.join(path_to_assets, "DSC_lookingForCsat_q=0", "central_last.csv"))
dataq0=formatData1(dataq0)
c_sat_q0=searchForCsat(dataq0)

print('... DSC(+) cells w fpass=0, get cell speed at c_sat')
data_Speed_q0=pd.read_csv(os.path.join(path_to_assets, "DSC_lookingForCsat_q=0_Zspeed", "central_vd.csv"))
data_Speed_q0=formatDataSpeed(data_Speed_q0)
c_sat_q0 = c_sat_q0.merge(data_Speed_q0, on=['celerity', 'wavelength'], how='left', indicator=True)

print('... DSC(+) cells w fpass=0, plot c_sat and speed at c_sat')
plotCritical(c_sat_q0, y='celerity', marker="o", ax=ax4, plot_legend=False, alpha=1., th=(False, 0., r''))
plotCritical(c_sat_q0, y='speedCell', marker="o", ax=ax8, plot_legend=False, alpha=1., th=(False, 0., r''))


"""
###########
print('\n... DSC(+) cells w fpass=0 & vinst=0,  get c_sat')
# on cherche c_sat à partir du max(BT) par rapport à c
data=pd.read_csv(os.path.join(path_to_assets, "DSC_lookingForCsat_vinst=0_q=0_Zspeed", "central_last.csv"))
data=formatData1(data)
c_sat=searchForCsat(data)

print('... DSC(+) cells w fpass=0 & vinst=0, get cell speed at c_sat')
data_Speed=pd.read_csv(os.path.join(path_to_assets, "DSC_lookingForCsat_vinst=0_q=0_Zspeed", "central_vd.csv"))
data_Speed=formatDataSpeed(data_Speed)
c_sat = c_sat.merge(data_Speed, on=['celerity', 'wavelength'], how='left', indicator=True)

print('... DSC(+) cells w fpass=0 & vinst=0, plot c_sat and speed at c_sat')
plotCritical(c_sat, y='celerity', marker="x", ax=ax4, plot_legend=False, alpha=1., th=(False, 0., r''))
plotCritical(c_sat, y='speedCell', marker="x", ax=ax8, plot_legend=False, alpha=1., th=(False, 0., r''))

"""


print("\n\n----------------------------------")
if save:
    print("... exporting figure")
    plt.savefig(path_to_save)
    print("... figure is exported")


print("\nPlot Figure 3 END OF PROCESS5\n---------------------------------------\n")
