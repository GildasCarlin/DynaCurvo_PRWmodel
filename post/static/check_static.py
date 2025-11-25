# ================================================================================================================
# ================================================================================================================
#                                                    packages importation
# ================================================================================================================
# ================================================================================================================

import os
import json

import regex as re

import pandas as pd
import numpy as np

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

NbLignes, NbColonnes=5, 4
Largeur, Hauteur = 25, 30
fig = plt.figure(num='curved_analysis', figsize=(Largeur, Hauteur), constrained_layout=True)
gs = GridSpec(NbLignes, NbColonnes, figure=fig)

## axes pour plot les trajectoires exp
ax_traj_1 = fig.add_subplot(gs[0, 0])
ax_traj_2 = fig.add_subplot(gs[0, 1])
ax_traj_3 = fig.add_subplot(gs[0, 2])
ax_traj_4 = fig.add_subplot(gs[0, 3])

## axes pour plot les trajectoires num
ax_traj_5 = fig.add_subplot(gs[1, 0])
ax_traj_6 = fig.add_subplot(gs[1, 1])
ax_traj_7 = fig.add_subplot(gs[1, 2])
ax_traj_8 = fig.add_subplot(gs[1, 3])

## axes pour plot les autres indicateurs (Dx, velocity, MSD etc)
ax1 = fig.add_subplot(gs[2, :-2])
ax2 = fig.add_subplot(gs[2, 2:])
ax3 = fig.add_subplot(gs[3, :-2])
ax4 = fig.add_subplot(gs[3, 2:])

## add labels

def add_label(ax, txt, x=-.05, y=1.1):
    ax.text(x, y, txt, transform=ax.transAxes, fontsize=20, fontweight='bold', va='top', ha='right')
#

ax_traj_1.text(-0.25, .85, 'Experimental traj.', transform=ax_traj_1.transAxes, rotation=90., fontsize=20, fontweight='bold', va='top', ha='right')
ax_traj_5.text(-0.25, .85, 'Numerical traj.', transform=ax_traj_5.transAxes, rotation=90., fontsize=20, fontweight='bold', va='top', ha='right')

add_label(ax_traj_1, 'A.1')
add_label(ax_traj_2, 'A.2')
add_label(ax_traj_3, 'A.3')
add_label(ax_traj_4, 'A.4')
add_label(ax_traj_5, 'B.1')
add_label(ax_traj_6, 'B.2')
add_label(ax_traj_7, 'B.3')
add_label(ax_traj_8, 'B.4')

add_label(ax1, 'C')
add_label(ax2, 'D')
add_label(ax3, 'E', y=1.15)
add_label(ax4, 'F')

## font size
#plt.rcParams.update({'font.size': 50})


# ================================================================================================================
# ================================================================================================================
#                                                  intermediate functions
# ================================================================================================================
# ================================================================================================================

def search_files_to_analyse(root, cond = 'data.txt'):
    paths_to_files = []

    for root_folder, sub_folder, files in os.walk(root):
        for file in files:
            if re.search(cond, file):
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

def format_traj(data, source='num', center=False):

    # trajectory are centered at origin.
    if center:
        for i in data[id[source][0]].unique():
            cell_i = data.loc[data[id[source][0]]==i,[id[source][1],id[source][2]]].copy()
            data.loc[data[id[source][0]]==i,[id[source][1],id[source][2]]] -= cell_i.iloc[0]
        #
    #

    if source=='num':
        # js frame vs python frame: the y-axis is reversed
        data[id[source][2]] *= -1.
    #
#

def plot_traj(metadata, data, IDlist, mapping, source='num'):
    if metadata['Surface']=='flat':
        return
    #

    # select axe to plot
    wl = metadata['Wavelength']
    ax = mapping[wl]
    
    sns.scatterplot(data=data[data[id[source][0]].isin(IDlist)], 
                    x=id[source][1], 
                    y=id[source][2], 
                    hue=id[source][0], 
                    ax=ax, 
                    palette='deep', 
                    marker='o', 
                    edgecolor='none'
                    )
    #
    
    ax.get_legend().set_visible(False)
    ax.set_xlabel('x (µm)', fontsize=17)
    ax.set_ylabel('y (µm)', fontsize=17)

    if source=='exp':
        ax.text(.65, 1.07, 
            r'$\lambda$={}µm'.format(wl), 
            transform=ax.transAxes, 
            fontsize=20, fontweight='bold', va='top', ha='right'
            )
        #
    #

    yticks = [-800., -600., -300., 0., 300., 300., 800.]
    xticks = [-800., -600., -300., 0., 300., 600.,800.]
    ax.set_yticks(yticks) 
    ax.set_xticks(xticks)
    ax.set_xlim(xticks[0], xticks[-1])
    ax.set_ylim(yticks[0], yticks[-1])
#

def get_Dx(metadata, data, source='num'):
    Dx_by_cell = []
    for i in data[id[source][0]].unique(): # calculates Dx for each cell
        Sx, Sy = data.loc[data[id[source][0]]==i,[id[source][1],id[source][2]]].diff().abs().sum(axis=0)
        Dx_by_cell.append(Sx/(Sx+Sy))
    #

    df = pd.DataFrame(Dx_by_cell, columns=['dx'])
    df.insert(0, 'celerity', len(df)*[metadata['Wave Celerity']])
    df.insert(0, 'wavelength', len(df)*[int(metadata['Wavelength'])])
    df.insert(0, 'source', len(df)*[source])
    return df
# 

def plot_distribution(data, quantity = 'dx', ax=ax1, plot_legend=True):
    spec = {'dx':r'$D_x$', 'speed':'Velocity (µm/s)'}

    if quantity=='dx':
        ax.axhline(y=.5, linewidth = 2, color='#8C0335', linestyle='--')
        ax.set_ylim(bottom=0., top=1.)
    #

    sns.violinplot(
        data=data, 
        x='wavelength', 
        y=quantity, 
        hue ='source', 
        order = ['20', '40', '80', '160'],
        bw_adjust=.5, 
        cut=1, 
        linewidth=1, 
        density_norm='area',
        palette=['#01B085','#2D568C'],
        ax=ax,
        legend=plot_legend
        )
    #

    ax.set_ylabel(spec[quantity], fontsize=25)
    ax.set_xlabel('Wavelength (µm)', fontsize=25)
    ax.tick_params(axis='both', which='major', labelsize=20)
    if plot_legend:
        mpatch1 = mpatches.Patch(facecolor='#01B085', label='T-cells', edgecolor='black')
        mpatch2 = mpatches.Patch(facecolor='#2D568C', label='Model', edgecolor='black')
        ax.legend(handles=[mpatch1, mpatch2], fontsize=20)
    #
#

def get_vd_num(metadata, data):
    dt = metadata['Time Interval']
    list_of_cells = data[id['num'][0]].unique()

    for i in list_of_cells:
        # -------------------- velocity coordinates
        df = data.loc[data['cell']==i, ['x', 'y']].copy().diff().dropna().div(dt).reset_index(drop=True)
        df = pd.DataFrame(df)
        df.columns = ['vx', 'vy']
        df.insert(0, 'celerity', len(df)*[metadata['Wave Celerity']])
        df.insert(0, 'wavelength', len(df)*[int(metadata['Wavelength'])])
        df.insert(0, 'cell', len(df)*[i])
        df.insert(0, 'source', len(df)*['num'])

        # -------------------- speed
        df['speed']=(df['vx'].pow(2) + df['vy'].pow(2)).pow(.5)

        # -------------------- velocity normalized coordinates
        df['nvx'] = df['vx'].div(df['speed'])
        df['nvy'] = df['vy'].div(df['speed'])
        ###df['check_unit_length'] = (df['nvx'].pow(2) + df['nvy'].pow(2)).pow(.5)

        if i==list_of_cells[0]:
            v=df.copy()
        else:
            v=pd.concat([v, df], ignore_index=True)
        #
    return v
#

def get_vd_exp(metadata, data):
    df = data.loc[data['Time Interval']!=0.0, ['Object', 'Time Interval', 'Delta X', 'Delta Y']].copy()

    # -------------------- velocity coordinates
    df['Delta X'] = df['Delta X'].div(df['Time Interval'])
    df['Delta Y'] = df['Delta Y'].div(df['Time Interval'])
    df.rename(columns={'Delta X': 'vx', 'Delta Y': 'vy', 'Object': 'cell'}, inplace=True)

    # -------------------- velocity coordinates
    df['speed']=(df['vx'].pow(2) + df['vy'].pow(2)).pow(.5)

    # -------------------- velocity normalized coordinates
    df['nvx'] = df['vx'].div(df['speed'])
    df['nvy'] = df['vy'].div(df['speed'])
    ###df['check_unit_length'] = (df['nvx'].pow(2) + df['nvy'].pow(2)).pow(.5)

    df.drop(columns=['Time Interval'], inplace=True)

    df.insert(1, 'celerity', len(df)*[metadata['Wave Celerity']])
    df.insert(1, 'wavelength', len(df)*[int(metadata['Wavelength'])])
    df.insert(0, 'source', len(df)*['exp'])
    df.fillna(0.0, inplace=True)
    return df
#

def get_vd(metadata, data, source='num'):
    if source=='num':
        return get_vd_num(metadata, data)
    #

    if source=='exp':
        return get_vd_exp(metadata, data)
    #
#

def get_MSD(data, dt=15., source='num'):
    msd_by_cell = []

    for i in data[id[source][0]].unique():
        cell_i = data.loc[data[id[source][0]]==i, [id[source][1], id[source][2]]].copy()

        traj_len = len(cell_i)
        last_point = 1+int(np.ceil(traj_len/3.))

        msd_by_tau = []
        for m in range(1, last_point):
            msd_by_tau.append(cell_i.diff(periods=m).dropna().pow(2).sum(axis=1).mean(axis=0))
        #
        msd_by_cell.append(msd_by_tau)
    #

    msd_cell = pd.DataFrame(msd_by_cell).T
    msd_cell['moy_pop'] = msd_cell.mean(axis=1)
    msd_cell.insert(0, 'tau', dt * np.arange(1., 1.+msd_cell.shape[0]))
    return msd_cell
#

def plot_time_function(wl, data_num, data_exp, quantity='msd', lim=200., ax=ax3):
    spec = {'msd':r'MSD (µm$^2$)', 'vacf':'Normalized VACF'}
    color = {'20':'#FFE529', '40':'#01B085', '80':'#2D568C', '160':'#4F1167'}

    data_exp_to_plot = data_exp.loc[data_exp['tau']<lim, ['tau', 'moy_pop']]
    data_num_to_plot = data_num.loc[data_num['tau']<lim, ['tau', 'moy_pop']]

    ax.scatter(
        data_exp_to_plot['tau'], 
        data_exp_to_plot['moy_pop'], 
        color=color[str(int(wl))], 
        edgecolor='black', 
        marker='o', 
        label='T-cells', 
        s = 100, 
        zorder=2
        )
    
    ax.scatter(
        data_num_to_plot['tau'], 
        data_num_to_plot['moy_pop'],
        color=color[str(int(wl))], 
        edgecolor='black', 
        marker='^', 
        label = 'Model', 
        s = 100, 
        zorder=3,
        )
    #

    if path==paths_to_simulated_data[0]:
        if quantity=='msd':
            handles, labels = plt.gca().get_legend_handles_labels()
            point1 = Line2D([0], [0], label='T-cells', marker='o', markersize=10, markerfacecolor='black', linestyle='', markeredgecolor='white')
            point2 = Line2D([0], [0], label='Model', marker='^', markersize=10, markerfacecolor='black', linestyle='', markeredgecolor='white')
            mpatch1 = mpatches.Patch(facecolor='#FFE529', label='20', edgecolor='black')
            mpatch2 = mpatches.Patch(facecolor='#01B085', label='40', edgecolor='black')
            mpatch3 = mpatches.Patch(facecolor='#2D568C', label='80', edgecolor='black')
            mpatch4 = mpatches.Patch(facecolor='#4F1167', label='160', edgecolor='black')
            handles.extend([point1, point2, mpatch1, mpatch2, mpatch3, mpatch4])
            ax.legend(handles=handles, fontsize=20)
        #

        ax.set_xlabel(r'$\Delta t$ (s)', fontsize=20)
        ax.set_ylabel(spec[quantity], fontsize=25)
        
        ax.tick_params(axis='both', which='major', labelsize=20)
    #
#

def get_VACF(vd, dt=15., source='num'):
    data = vd.loc[vd['source']==source, ['cell', 'nvx', 'nvy']].copy()
    
    vacf_by_cell = []
    for i in data['cell'].unique():
        cell_i = data.loc[data['cell']==i, ['nvx', 'nvy']].copy()

        N_seg = len(cell_i)
        last_point = 1+int(np.ceil(N_seg/3.))

        vacf_by_tau = []
        for m in range(1, last_point):
            vxy = cell_i[:-m].copy().reset_index(drop=True)
            vxy.columns = ['V0X', 'V0Y']
            vxy[['VmX','VmY']] = cell_i[m:].copy().reset_index(drop=True)
            
            vacf_by_tau.append((vxy['V0X']*vxy['VmX'] + vxy['V0Y']*vxy['VmY']).mean(axis=0))
        #
        vacf_by_cell.append(vacf_by_tau)

    vacf_cell = pd.DataFrame(vacf_by_cell).T
    vacf_cell['moy_pop'] = vacf_cell.mean(axis=1)
    vacf_cell.insert(0, 'tau', dt * np.arange(1., 1.+vacf_cell.shape[0]))

    return vacf_cell
#   


# ================================================================================================================
# ================================================================================================================
#                                                    User interface
# ================================================================================================================
# ================================================================================================================

## nb of cell trak to plot
NtoPlot=30

# mapping between simulation parameters and ax to plot traj
mapping_sim_traj = {
    20.: ax_traj_5,
    40.: ax_traj_6,
    80.: ax_traj_7,
    160.:ax_traj_8
}

# mapping between experimental parameters and ax to plot traj
mapping_exp_traj = {
    20.: ax_traj_1,
    40.: ax_traj_2,
    80.: ax_traj_3,
    160.:ax_traj_4
}

# columns labels
id = {
    'num':('cell', 'x', 'y'), 
    'exp':('Object', 'X', 'Y')
}

# export figure or not
save = True

# parameters
Jacobian='J=gradH'
f='med-flat'
param = 'q'
param_value=0.

# ================================================================================================================
# ================================================================================================================
#                                                    MAIN function
# ================================================================================================================
# ================================================================================================================


principal_path = os.getcwd()+'\\assets\\simulation_data\\'+Jacobian+'\\'+f+'\\'+param+'='+str(param_value)
paths_to_simulated_data = search_files_to_analyse(principal_path)

for path in paths_to_simulated_data:
    # get data and metadata
    print('\n------------------------------------ open data: ...', path[-25:])
    ## load numerical data
    metadata, data = load_trajectory_data(path, ['Surface Index', 'Number of Cells'])
    wl, c = metadata['Wavelength'], metadata['Wave Celerity']
    
    ## load experimental data
    path_to_corresponding_exp_data = search_files_to_analyse(os.getcwd()+'\\assets\\experimental_data', cond=str(int(wl))+'.csv')[0]
    metadata_exp, data_exp = load_trajectory_data(path_to_corresponding_exp_data, ['Surface Index', 'Number of Cells'])

    print('------------------------------------ data analysis: wavelength={}µm, celerity={}µm/s'.format(wl, c))
    print('... trajectories')
    ## plot numerical data
    format_traj(data, center=True)
    plot_traj(metadata, data, np.random.choice(data['cell'].unique(), size=NtoPlot, replace=False), mapping_sim_traj)

    # plot experimental data
    format_traj(data_exp, source='exp', center=True)
    plot_traj(metadata_exp, data_exp, np.random.choice(data_exp['Object'].unique(), size=NtoPlot, replace=False), mapping_exp_traj, source='exp')

    print('... dx, vd')
    # compute dx
    dx_num = get_Dx(metadata, data)
    dx_exp= get_Dx(metadata_exp, data_exp, source='exp')
    dx_by_topo = pd.concat([dx_exp, dx_num], ignore_index=True) 
    
    # compute vd
    vd_num = get_vd(metadata, data)
    vd_exp= get_vd(metadata_exp, data_exp, source='exp')
    vd_by_topo = pd.concat([vd_exp, vd_num], ignore_index=True)

    ## rearrange data
    if path==paths_to_simulated_data[0]:
        dx=dx_by_topo.copy()
        vd=vd_by_topo.copy()
    else:
        dx = pd.concat([dx, dx_by_topo], ignore_index=True)  
        vd = pd.concat([vd, vd_by_topo], ignore_index=True)
    #

    print('... msd')
    plot_time_function(wl, get_MSD(data), get_MSD(data_exp, source='exp'))

    print('... vacf')
    plot_time_function(wl, get_VACF(vd_num), get_VACF(vd_exp, source='exp'), quantity='vacf', ax=ax4)
#

print('\n')
print('------------------------------------ final plots:')
plot_distribution(dx, quantity = 'dx')
plot_distribution(vd, quantity = 'speed', ax=ax2, plot_legend=False)




if save:
    print('------------------------------------ export figure:')
    plt.savefig(principal_path+'\\check_static.svg')
#