# ================================================================================================================
# ================================================================================================================
#                                                    packages importation
# ================================================================================================================
# ================================================================================================================

import os
import time
import regex as re

import numpy as np


# ================================================================================================================
# ================================================================================================================
#                                                    User interface
# ================================================================================================================
# ================================================================================================================


## ------------------------------------ folder and paths

# study
substrate_regime="dyna"  # dyna or static
study="SpeedDistribution"

# conditions in the study
Assumptions=["DSC", "NDSC"]
Lambdas=[20.0]
celerities=[0.0, 0.1, 0.4, 0.6, 2.0]
max_c, min_c = 0., 4.0
c_step=.1

# exported file (at the end of the simulation)
home=os.path.expanduser("~")
path_to_export=os.path.join(home, "Downloads", "data.txt")

# model version + path to template param_user.js
path_to_model = os.path.join("model")
path_to_UI = os.path.join("..", "model", "param_user.js")

# assets
path_to_assets=os.path.join("..", "post", substrate_regime, "assets", study)


## ------------------------------------ simulations parameters

# assumptions:    cells are sensitive to curvature (K_assumption) ? cells are sensitive to the wave (W_assumption) ?
K_assumption = r"true"

# substrate pattern:   right_sinus  /  convr_sinus  /   wells_sinus  /  snake_sinus
pattern = "right_sinus"

amplitude=10.0
initial_phase=0.0
angular_velocity=0.0
initial_angle="-Math.PI/2.0"

# condition
simu_length=600.0  # seconds
sampling_rate=0.1

# others
nbCell=150
export_data=r"true"
draw_cells=r"false"


# ================================================================================================================
# ================================================================================================================
#                                                    intermediate functions
# ================================================================================================================
# ================================================================================================================

def createFolders(root, Wavelengths, Assumptions, Celerities=None, min_c=None, max_c=None, step=None):
    # method range: celerities are given by a range and a step, the list is constructed
    # method list: the list of celerities is explicitly given
    
    if Celerities is not None:
        Celerities=Celerities
    elif all(v is not None for v in [min_c, max_c, step]):
        Celerities=list(np.arange(min_c, max_c + step, step))
    Celerities=[round(c, 2) for c in Celerities]

    for assumption in Assumptions:
        for wavelength in Wavelengths:
            for celerity in Celerities:
                path=os.path.join(root, assumption, "L="+str(round(wavelength, 1)), "c="+str(round(celerity, 2)))
                if not os.path.exists(path):
                    os.makedirs(path)
            #
        #
    #
#

def searchFilesToAnalyse(root):
    paths_to_files = []
    for root_folder, sub_folder, files in os.walk(root):
        if re.search('c=', root_folder) and not ('data.csv' in files or 'data.txt' in files):
            paths_to_files.append(os.path.join(root_folder))
    return paths_to_files
#

def prepareUI():
    # this function edits param_user.js but only parameters which are commun for all simulations

    # open UI before editing parameters
    temp_file = open(path_to_UI)
    temp_content = temp_file.read()
    temp_file.close()

    # edit parameters in UI
    temp_content = re.sub(r'let\scurvature_assumption\s*=\s*(true|false);', r'let curvature_assumption='+K_assumption+r';', temp_content)

    temp_content = re.sub(r"let\ssurface_title\s*=\s*(\'|\")\w*(\'|\");", "let surface_title=\'"+pattern+"\';", temp_content)
    temp_content = re.sub(r"let\sa\s*=\s*(\d*\.?\d*);", r'let a='+str(amplitude)+r';', temp_content)
    temp_content = re.sub(r"let\sphi0\s*=\s*(\d*\.?\d*);", r'let phi0='+str(initial_phase)+r';', temp_content)
    temp_content = re.sub(r"let\sw\s*=\s*(\d*\.?\d*);", r'let w='+str(angular_velocity)+r';', temp_content)
    temp_content = re.sub(r"let\spsi0\s*=\s*(\d*\.?\d*);", r'let psi0='+initial_angle+r';', temp_content)

    temp_content = re.sub(r"let\snbCell\s*=\s*(\d*\.?\d*);", r'let nbCell='+str(nbCell)+r';', temp_content)
    temp_content = re.sub(r"let\sTmax\s*=\s*(\d*\.?\d*);", r'let Tmax='+str(simu_length)+r';', temp_content)
    temp_content = re.sub(r"let\sdt_mesure\s*=\s*(\d*\.?\d*);", r'let dt_mesure='+str(sampling_rate)+r';', temp_content)
    temp_content = re.sub(r'let\sdata\s*=\s*(true|false);', r'let data='+export_data+r';', temp_content) 
    temp_content = re.sub(r'let\sdraw_cells\s*=\s*(true|false);', r'let draw_cells='+draw_cells+r';', temp_content)

    # save new UI
    UI_to_write = open(path_to_UI, 'w')
    UI_to_write.write(temp_content)
    UI_to_write.close()
#

def editUI(ass, wl, c):
    # this function edits param_user.js but only parameters which change at each simulation simulations

    def getAssumption(assumption):
        mapping={"NDSC":r"false", "DSC":r"true"}
        return mapping[assumption]
    #
    
    # open UI before editing parameters
    temp_file = open(path_to_UI)
    temp_content = temp_file.read()
    temp_file.close()

    # edit parameters in UI
    temp_content = re.sub(r'let\ssurf_assumption\s*=\s*(true|false);', r'let surf_assumption='+getAssumption(ass)+r';', temp_content) 

    temp_content = re.sub(r'let\sc\s*=\s*(\d*\.?\d*);', r'let c='+str(c)+r';', temp_content)
    temp_content = re.sub(r'let\slambda\s*=\s*(\d*\.?\d*);', r'let lambda='+str(wl)+r';', temp_content)

    # save new UI
    UI_to_write = open(path_to_UI, 'w')
    UI_to_write.write(temp_content)
    UI_to_write.close()
#


# ================================================================================================================
# ================================================================================================================
#                                                    MAIN function
# ================================================================================================================
# ================================================================================================================

print("\nLaunch Multiple Simulations Process")
print("==========================================================\n")
print("substrate regime: " + substrate_regime + "\nstudy: " + study)

createFolders(path_to_assets, Lambdas, Assumptions, Celerities=celerities)

path_to_simu = searchFilesToAnalyse(path_to_assets)
total_number_of_simu=len(path_to_simu)
print("Number of simulations to run:", total_number_of_simu)

# parameters which are commun to all simulations
prepareUI()

count=0
for path in path_to_simu:
    count+=1
    print("\nSimulation number " + str(count) + " out of " + str(total_number_of_simu) + " (" + str(round(100.*count/total_number_of_simu, 2)) + "%)")
    print("------------------------------------------")

    print("...... preparation")
    # parameters to edit
    (ass, wl, c) = re.findall(r'(DSC|NDSC)\\L=(\d+\.\d+|\d+)\\c=(\d+\.\d+|\d+)', path)[0]
    
    print('...  assumption {}, wavelength {}µm, wave celerity {}µm/s'.format(ass, wl, c))
    # parameters which are specific to this simulation
    editUI(ass, wl, c)

    # run simulation
    os.system(r"start http://localhost:8000//" + path_to_model)
    print('...... simulation is running')

    # check when the simulation is done
    simu_is_done = os.path.exists(path_to_export)
    while (not simu_is_done):
        time.sleep(0.1)
        simu_is_done = os.path.exists(path_to_export)
    #

    ## when the simulation is done, move the data file to corresponding data analysis folder ...
    os.system(r'move "'+path_to_export + r'" "' + path)
    print('...... simulation is done, files are exported')

    ##... and close simulation
    os.system(r'taskkill /f /im msedge.exe')
#

print("\nLaunch Multiple Simulations Process DONE")
print("==========================================================\n")
