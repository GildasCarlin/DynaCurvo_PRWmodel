# ================================================================================================================
# ================================================================================================================
#                                                    packages importation
# ================================================================================================================
# ================================================================================================================

import json
import os
import pandas as pd

import scipy.optimize as opt
from scipy.stats import t
import numpy as np

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.gridspec import GridSpec


# ================================================================================================================
# ================================================================================================================
#                                                    open figure
# ================================================================================================================
# ================================================================================================================

NbLignes, NbColonnes=1, 2
Largeur, Hauteur = 14, 7
fig = plt.figure(num='S_stat', figsize=(Largeur, Hauteur), constrained_layout=True)
gs = GridSpec(NbLignes, NbColonnes, figure=fig)

plt.rcParams.update({'font.size': 18})
## axes pour plot les autres indicateurs (Dx, velocity, MSD etc)
ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[0, 1])


# ================================================================================================================
# ================================================================================================================
#                                                  intermediate functions
# ================================================================================================================
# ================================================================================================================

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

def get_vd(metadata, data):
    # collect data
    wl, c = metadata['Wavelength'], metadata['Wave Celerity']

    df = data.loc[data['Time Interval']!=0.0, ['Object', 'Time Interval', 'Delta X', 'Delta Y']].copy()
    # -------------------- velocity coordinates
    df['Delta X'] = df['Delta X'].div(df['Time Interval'])
    df['Delta Y'] = df['Delta Y'].div(df['Time Interval'])
    df.rename(columns={'Delta X': 'vx', 'Delta Y': 'vy', 'Object': 'cell'}, inplace=True)

    # -------------------- velocity norm
    df['speed']=(df['vx'].pow(2) + df['vy'].pow(2)).pow(.5)

    # -------------------- velocity normalized coordinates
    df['nvx'] = df['vx'].div(df['speed'])
    df['nvy'] = df['vy'].div(df['speed'])
    df['check_unit_length'] = (df['nvx'].pow(2) + df['nvy'].pow(2)).pow(.5)

    df.drop(columns=['Time Interval'], inplace=True)

    n = len(df)
    df.insert(1, 'celerity', n*[c])
    df.insert(1, 'wavelength', n*[wl])
    df.insert(0, 'source', n*['exp'])
    df.fillna(0.0, inplace=True)
    return df
#

def get_CI(metadata, data, conf=.95, central='med'):
    # collect data
    wl = metadata['Wavelength']
    df = data['speed'].copy()
    n = len(df) 
    mea, med, s = df.mean(), df.median(), df.var()

    print('....... topo sur la distribution des vitesses (exp)')
    print("... longueur d'onde={}µm".format(wl))
    print("... taille de l'échantillon={}".format(n))
    print('... médiane={}µm/s, moyenne={}µm/s'.format(med, mea))

    proba = (1.-conf)*100. ; proba = (100.-proba/2.)/100. 
    ddl = n - 1
    intervalle = np.sqrt(s/n) * t.ppf(proba, ddl)
    
    if central=='med':
        return (med, intervalle)
    elif central=='mea':
        return (mea, intervalle)
#

def assumed(theta, x=np.array([20., 40., 80., 160.])):
    A, tau, b = theta[0], theta[1], theta[2]
    return A*np.exp(-1.*x/tau)+b
# 

def fun(theta):  
    return list(reversed(s_stat)) - assumed(theta)
#

# ================================================================================================================
# ================================================================================================================
#                                                    User interface
# ================================================================================================================
# ================================================================================================================

## compute the mean (mea) or the median (med)
central='mea'

## metadata that must be string type
to_string=['Surface Index', 'Number of Cells']

## experimental condition to run to
exps=['flat', 'sin160', 'sin80', 'sin40', 'sin20']
x=[160., 80., 40., 20.]

## show figure or not
sh = True

# ================================================================================================================
# ================================================================================================================
#                                                    MAIN function
# ================================================================================================================
# ================================================================================================================

# source folder
root=os.getcwd()

CI_in_dict = {}
cent_in_dict = {}
for exp in exps:

    ## open experimental data
    print("\n------------------------------------ open experimental data")
    metadata, data = load_trajectory_data(
        root+"\\assets\\experimental_data\\" + exp + ".csv", 
        to_string
        )
    #
    
    print("------------------------------------ compute the velocities")
    vd=get_vd(metadata, data)

    if exp==exps[0]:
        all_v = vd.copy()
    else:
        all_v=pd.concat([all_v, vd])
    #

    print("------------------------------------ compute the confidence interval")
    cent_in_dict[exp], CI_in_dict[exp] = get_CI(metadata, vd, central=central)
#

###### preparation to plot
cent_flat = cent_in_dict.pop('flat', None)
CI_flat = CI_in_dict.pop('flat', None)

list_central=np.array(list(cent_in_dict.values()))
list_CI=np.array(list(CI_in_dict.values()))

xbornes=[10., 170.]

##### plot velocity distribution
ax2.hist(all_v['speed'], bins=50, density=True)
ax2.set_xlim(left=0.)


###### plot median
ax1.scatter(x, list_central, color='black', marker='s')
ax1.axhline(y=cent_flat, color='#C03052', linestyle='--')

###### plot CI
m_plus = list_central + list_CI
m_minus = list_central - list_CI
ax1.plot((x,x), (m_plus,m_minus),'-',color='black')
ax1.fill_between(xbornes, 
                 len(xbornes)*[cent_flat+CI_flat], 
                 len(xbornes)*[cent_flat-CI_flat], 
                 color='#C03052',
                 alpha=.2)
#


###### plot S_stat
s_stat = np.abs(list_central - cent_flat)
print("------------------------------------ experimental values of s_stat")
print('s_stat =',s_stat)
ax1.scatter(x, s_stat, color='blue', s=150)


### perform an exponential fit
theta0 = [3., 3., .5]
res = opt.least_squares(fun, theta0, bounds=([0., 0., 0.], [np.inf, np.inf, np.inf])) 
RMSE=np.sqrt(np.mean((np.square(res.fun))))

print("\n------------------------------------ fit of s_tat")
print("... A={}µm/s, tau={}µm, b={}µm/s".format(res.x[0], res.x[1], res.x[2]))
print('... residuals =', res.fun)
print("... RMSE={}µm/s".format(RMSE))

# plot fit
X=np.linspace(xbornes[0], xbornes[1], 50)
ax1.plot(X, assumed(res.x, x=X), color='blue', linestyle='--')
ax1.text(100., .05, r"RMSE={}µm/s".format(round(RMSE, 3)), fontsize=14)

###### beautiful figure
ax1.set_xlabel("Wavelength (µm)", fontsize=20)
ax1.set_ylabel("Cell velocity (µm/s)", fontsize=20)

y0, y1 = ax1.get_ylim()
ax1.set_ylim(0., y1+.05)
ax1.set_xlim(left=xbornes[0], right=xbornes[1])

###### beautiful legend
handles, labels = plt.gca().get_legend_handles_labels()
point1 = Line2D([0], [0], label=r"mean velocity $\overline{v}_{\lambda}$", marker='s', markersize=10, markerfacecolor='black', linestyle='', markeredgecolor='white')
point2 = Line2D([0], [0], label=r"$|\overline{v}_{\lambda} - \overline{v}_{flat}|$", marker='o', markersize=10, markerfacecolor='blue', linestyle='', markeredgecolor='white')
point3 = Line2D([0], [0], label=r"$s_{stat}$ exponential fit", color='blue', markeredgecolor='white', linestyle='--')
handles.extend([point1, point2, point3])
ax1.legend(handles=handles, fontsize=15)

if sh:
    plt.show()
#


