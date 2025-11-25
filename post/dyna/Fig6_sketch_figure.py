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
fig1 = plt.figure(num='role_of_f_stat', figsize=(Largeur, Hauteur), constrained_layout=True)
gs1 = GridSpec(NbLignes, NbColonnes, figure=fig1)

## axes pour scatter les trajectoires
ax_traj1 = fig1.add_subplot(gs1[0, 0])
ax_traj2 = fig1.add_subplot(gs1[0, 1])

## axes pour les statistiques
ax1 = fig1.add_subplot(gs1[1, 0])
ax2 = fig1.add_subplot(gs1[1, 1])
ax3 = fig1.add_subplot(gs1[0, 2])
ax4 = fig1.add_subplot(gs1[1, 2])
ax5 = fig1.add_subplot(gs1[2, 0])
ax6 = fig1.add_subplot(gs1[2, 1])
ax7 = fig1.add_subplot(gs1[2, 2])

ax_ins2 = inset_axes(ax2, width=Largeur_ins, height=Hauteur_ins)
ax_ins5 = inset_axes(ax5, width=Largeur_ins, height=Hauteur_ins)

## add labels
def add_label(ax, txt, x=-0.05, y=1.1):
    ax.text(x, y, txt, transform=ax.transAxes, fontsize=20, fontweight='bold', va='top', ha='right')
#

add_label(ax_traj1, 'A')
add_label(ax_traj2, 'B')
add_label(ax1, 'C')
add_label(ax2, 'D')
add_label(ax3, 'E')
add_label(ax4, 'F')
add_label(ax5, 'G')
add_label(ax6, 'H')
add_label(ax7, 'I')

## font size
plt.rcParams.update({'font.size': 50})


# ================================================================================================================
# ================================================================================================================
#                                                    intermediate functions
# ================================================================================================================
# ================================================================================================================

def plot_traj(metadata, data, N, mapping, xticks, yticks):

    ## random cells id to plot
    IDlist = np.random.choice(data['cell'].unique(), size=N, replace=False)

    # select axe to plot
    ass, wl, c = metadata[0], metadata[1], metadata[2]
    ax = mapping[ass]
    
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
    
    ax.get_legend().set_visible(False)
    ax.set_xlabel('x (mm)', fontsize=22)
    ax.set_ylabel('y (mm)', fontsize=22)
    ax.tick_params(axis='both', which='major', labelsize=15)

    ax.set_yticks(yticks) 
    ax.set_xticks(xticks)
    ax.set_xlim(xticks[0], xticks[-1])
    ax.set_ylim(yticks[0], yticks[-1])
        
    ## print text
    ax.text(0.4, 0.05, r"$f_stat$=0", rotation='horizontal', fontsize=15, va='center', ha='center', transform=ax.transAxes)
    ax.text(0.6, 0.05, r"$f_stat$=1", rotation='horizontal', fontsize=15, va='center', ha='center', transform=ax.transAxes)
    ax.text(0.5, 0.95, ass + " cells", rotation='horizontal', fontsize=20, va='center', ha='center', transform=ax.transAxes)
    
    if ass=='NDSC':
        ax.arrow(0.08, 0.4, 0., 0.25, width=0.015, color='red', transform=ax.transAxes, clip_on=False, head_width=0.05, head_length=0.07) 
        ax.arrow(0.15, 0.92, -0.08, 0., width=0.015, color='red', length_includes_head=True, transform=ax.transAxes, clip_on=False)
        ax.arrow(0.15, 0.92, 0.12, 0., width=0.015, color='red', length_includes_head=True, transform=ax.transAxes, clip_on=False)
        
        ax.text(0.12, 0.55, r"Wave", transform=ax.transAxes, color='red', rotation = 'vertical', fontsize=22, fontweight='bold', va='center', ha='center')
        ax.text(0.17, 0.86, r"Grooves", color='red', rotation='horizontal', fontsize=22, fontweight='bold', va='center', ha='center', transform=ax.transAxes)
    #

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
        ax.legend(handles=handles, fontsize=25, title=r"$f_stat$", title_fontsize=28, loc='upper right')
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
    ax.set_xlabel("$f_stat$", fontsize=25)
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

def plot_critical(data, y, r2, ax, palette, plot_legend=False, xticks = [2.0, 4.0, 6.0], alpha=.5, th=(False, .79, 'th')):
    marker={'last_y':'^', 'dx':'o', 'vd':'s'}
    ylabel={'celerity':r"$c_{sat}$ (µm/s)", 'vd':r"$\overline{v}_{sat}$ (µm/s)"}
    
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

################### parameter
## boucle 1  (barycentres)
para1 = 'Py'
para1_values = [20.0]
para2='Vmax'
para2_values = [4.0]

assumptions = ['NDSC', 'DSC']
alternate_id = {"NDSC":'cell_is_not_surfer', "DSC":'cell_is_surfer'}

## boucle 2  (celerité critique)
para3='Py'
para3_values=[20.0]
para4='Vmax'
para4_values=[2.0, 4.0, 6.0]
assumptions_critic = ['DSC']

################### trajectories scatter parameters
mapping_traj={'NDSC': ax_traj1,
              'DSC': ax_traj2
            }
#

# x and y axis
xticks = [-2., -1., 0., 1., 2.]
yticks = [-.5, 0., 1., 2., 3.5]

################### statistics
mapping_stats={('NDSC', 'bary'):ax1,
               ('DSC', 'bary'):ax2,
                }
#

x_axis_celerity = {'NDSC':5., 'DSC':15.}
plot_sdyna_or_not={'NDSC':False, 'DSC':True}
plot_legend_or_not={'NDSC':True, 'DSC':False}

wavelengths_statistics=[20.]
wavelengths_critic = [20., 40., 80.]

################### compute r2 score
perform_stat = True
dependant_variable='vd'

################### export figure or not
save = True
path_to_save = 'role_of_' + para1 + '\\role_' + para1 + '.svg'


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
list_last, list_dx, list_vd = [], [], []
list_c_last, list_c_dx, list_c_vd = [], [], []

################################### barycenters Vmax=4
source=""
for assumption in assumptions:
    for value in para1_values:
        print('\n------------------------------------ open data')
        ## path to data
        if (para1=='Vmax' and value==4.0) or (para1=='q' and value==1.0):
            source=alternate_id[assumption]
        elif assumption=='NDSC':
            source='role_of_' + para1 + '\\' + para1 + '=' + str(value) + '\\' + alternate_id[assumption]
        elif assumption=='DSC':
            source='role_of_' + para1 + '\\' + para1 + '=' + str(value) + '\\' + alternate_id[assumption] + '\\' + para2 + '=' + str(para2_values[0])
        #
        source+="\\analysis_results\\global\\"
        print("... " + assumption + " " + para1 + "=" + str(value))

        path_to_last = root + "\\" + source + "central_last.csv"
        path_to_vd =   root + '\\' + source + "central_vd.csv"

        data_last = pd.read_csv(path_to_last)
        data_vd =   pd.read_csv(path_to_vd)

        print('------------------------------------ sort data')

        data_last = data_last[data_last['celerity'] <= x_axis_celerity[assumption]]
        data_vd =   data_vd[data_vd['celerity'] <= x_axis_celerity[assumption]]

        data_last.insert(0, para1, len(data_last)*[value])
        data_vd.insert(0, para1, len(data_vd)*[value])

        data_last.insert(0, "assumption", len(data_last)*[assumption])
        data_vd.insert(0, "assumption", len(data_vd)*[assumption])

        data_last = data_last[data_last['wavelength'].isin(wavelengths_statistics)]
        data_vd =   data_vd[data_vd['wavelength'].isin(wavelengths_statistics)]

        list_last.append(data_last)
        list_vd.append(data_vd)
    #
#

for value in para4_values:
    print('\n------------------------------------ open data')
    ## path to data
    source='role_of_' + para3 + '\\' + para3 + '=' + str(para3_values[0]) + '\\' + alternate_id[assumptions_critic[0]] + '\\' + para4 + '=' + str(value)
    
    source+="\\analysis_results\\global\\"
    print("... " + assumption + " " + para4 + "=" + str(value))

    path_to_vd =   root + '\\' + source + "central_vd.csv"
    data_vd =   pd.read_csv(path_to_vd)

    print('------------------------------------ extract max/min')
    critic_c_vd =   data_vd.loc[data_vd.groupby("wavelength")["vd"].idxmax(), ["wavelength", "celerity", "vd"]].reset_index(drop=True)
    critic_c_vd.insert(0, para4, len(critic_c_vd)*[value])
    critic_c_vd.insert(0, "assumption", len(critic_c_vd)*[assumptions_critic[0]])

    critic_c_vd = critic_c_vd[critic_c_vd['wavelength'].isin(wavelengths_critic)]
    critic_c_vd = critic_c_vd[critic_c_vd['assumption'].isin(assumptions_critic)]
    list_c_vd.append(critic_c_vd)
#

print('\n------------------------------------ concat data')
df_last = pd.concat(list_last, ignore_index=True)
df_last['last_y'] /= 1000.  ## µm -->  mm
df_vd = pd.concat(list_vd, ignore_index=True)

df_c_vd = pd.concat(list_c_vd, ignore_index=True)
df_c_vd = df_c_vd.sort_values(by=para4, ascending=True).reset_index(drop=True)

if perform_stat:
    print('\n------------------------------------ compute r2 and rmse')

    if dependant_variable == 'vd':
        data=df_c_vd.copy()
    #

    ## theoretical tables 
    x=np.linspace(2., 6., 3)
    y_true_c, y_true_a=2.**(-1./3.)*x, 2.**(2./3.)/3.*x
    
    r2_c, r2_a = [], []
    for wl in wavelengths_critic:
        print('\n... wavelength =', wl)
        df = data[(data['wavelength']==wl)]
        r2_c.append(r2_score(y_true_c, df['celerity']))
        print('...r2 for critic celerity', r2_c[-1])
        
        if dependant_variable=='vd':
            r2_a.append(r2_score(y_true_a, df[dependant_variable]))
            print('...r2 for average velocity at critic celerity', r2_a[-1])
        #
    #
#

if save:
    
    ######## récupérér les trajectoires ...
    #### NDSCS q=0 and q=1
    traj_NS0 = pd.read_csv(root + "\\role_of_Py\\Py=30.0\\cell_is_not_surfer\\lambda=20.0\\c=0.6\\data.txt", comment='#')
    traj_NS1 = pd.read_csv(root + "\\cell_is_not_surfer\\lambda=20.0\\c=0.6\\data.txt", comment='#')

    #### SDSC q=0 and q=1
    traj_S0 =  pd.read_csv(root + "\\role_of_Py\\Py=30.0\\cell_is_surfer\\Vmax=4.0\\lambda=20.0\\c=3.0\\data.txt", comment='#')
    traj_S1 = pd.read_csv(root + "\\cell_is_surfer\\lambda=20.0\\c=3.0\\data.txt", comment='#')

    print('---------------------- treat traj')
    ######## ... les convertir
    #### invert y-axis
    traj_NS0["y"] *= -1.
    traj_S0["y"]  *= -1.
    traj_NS1["y"] *= -1.
    traj_S1["y"]  *= -1.

    #### convert µm to mm
    traj_NS0[['x', 'y']] /= 1000.
    traj_S0[['x', 'y']]  /= 1000.
    traj_NS1[['x', 'y']] /= 1000.
    traj_S1[['x', 'y']]  /= 1000.

    #### shift traj on the LEFT of the panel (q=0)
    traj_NS0['x'] += -1.
    traj_S0['x'] += -1.

    #### shift traj on the RIGHT of the panel (q=1)
    traj_NS1['x'] += 1.
    traj_S1['x'] += 1.

    ######## scatter les traj
    print("---------------------- plot traj")
    #### NDSC q=0
    plot_traj(metadata=['NDSC', 20.0, 0.6], 
              data=    traj_NS0, 
              N=       30, 
              mapping= mapping_traj, 
              xticks=  xticks, 
              yticks=  yticks
              )
    
    #### DSC q=0
    plot_traj(metadata=['DSC', 20.0, 3.0], 
              data=    traj_S0, 
              N=       30, 
              mapping= mapping_traj, 
              xticks=  xticks, 
              yticks=  yticks
              )
    
    #### NDSC q=1
    plot_traj(metadata=['NDSC', 20.0, 0.6], 
              data=    traj_NS1, 
              N=       30, 
              mapping= mapping_traj, 
              xticks=  xticks, 
              yticks=  yticks
              )
    
    #### DSC q=1
    plot_traj(metadata=['DSC', 20.0, 3.0], 
              data=    traj_S1, 
              N=       30, 
              mapping= mapping_traj, 
              xticks=  xticks, 
              yticks=  yticks
              )
    # 

    print("---------------------- plot central statistics")
    assumption='DSC'
    data=df_vd[(df_vd['assumption']==assumption)]
    plot_central(data, 'vd',para1,ax = ax_ins2, palette=q_palette, plot_sdyna=True, x_sdyna=x_axis_celerity[assumption], plot_legend=False)
    #

    for assumption in assumptions:
        data=df_last[(df_last['assumption']==assumption)]
        plot_central(data,
                     'last_y',
                     para1,
                     ax = mapping_stats[(assumption, 'bary')],
                     palette=q_palette,
                     plot_legend=plot_legend_or_not[assumption],
                     )
        #
    #

    print("---------------------- plot critic")
    plot_critical(df_c_vd, y='celerity', r2=r2_c, ax=ax3, palette=wl_palette,  alpha=1., plot_legend=True, th=(True, 2.**(-1./3.), r'$w_{sat}$'))
    plot_critical(df_c_vd, y='vd', r2=r2_a, ax=ax4, palette=wl_palette,  alpha=1., plot_legend=True, th=(True, 2.**(2./3.)/3., r'$\max_c (s_{dyna})$'))

    print('\n------------------------------------ export figure:')
    plt.savefig(root + '\\' +  path_to_save)
#
        






