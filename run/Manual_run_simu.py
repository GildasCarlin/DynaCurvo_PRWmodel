# ================================================================================================================
# ================================================================================================================
#                                                    packages importation
# ================================================================================================================
# ================================================================================================================

import os
import time
import regex as re

import numpy as np

from map import *


# ================================================================================================================
# ================================================================================================================
#                                                    User interface
# ================================================================================================================
# ================================================================================================================


## ------------------------------------ folder and paths

# study
substrate_regime="dyna"  # dyna or static
study="NDSC_speed"

# exported file (at the end of the simulation)
home=os.path.expanduser("~")
path_to_export=os.path.join(home, "Downloads", "data.txt")

# model version + path to template param_user.js
path_to_model = os.path.join("model")
path_to_UI = os.path.join("..", "model", "param_user.js")

# assets
path_to_assets=os.path.join("..", "post", substrate_regime, "assets")

# TEST code (False) before actual run (True).
TEST=False

# ================================================================================================================
# ================================================================================================================
#                                                    intermediate functions
# ================================================================================================================
# ================================================================================================================

def editUI(simu_params):
    # this function edits param_user.js but only parameters which change at each simulation simulations

    # open UI before editing parameters
    with open(path_to_UI, 'r') as temp_file:
        temp_content = temp_file.read()

    for key, value  in simu_params.items():
        pattern = rf"(let\s+{re.escape(key)}\s*=\s*)([^;]+)(;)"

        replacement = rf"\g<1>{value}\g<3>"
        temp_content = re.sub(pattern, replacement, temp_content)
    
    # save new UI
    with open(path_to_UI, 'w') as UI_to_write:
        UI_to_write.write(temp_content)
#

def getSign(metadf):
    if metadf["surface_title"]==r"'convr_sinus'":
        return np.sign(metadf["w"])
    else:
        return None
#

def getPath(values):
    ass={r"true":"DSC", r"false":"NDSC"}
    mapp={(r"'convr_sinus'", 1.0):"converging_sinus",
          (r"'convr_sinus'", -1.0):"diverging_sinus",
          (r"'snake_sinus'", None):"Snake_sinus",
          (r"'wells_sinus'", None):"wells",
        }
    

    subvalues=(values["surface_title"], getSign(values))
    if subvalues[0]==r"'right_sinus'":
        return os.path.join(path_to_assets, study, f"L={values["lambda"]}", f"c={values["c"]}")
    
    else:
        return os.path.join(path_to_assets, "Other_surfaces", mapp[subvalues], ass[values["surf_assumption"]], f"c={values["c"]}_Zspeed")

#



    
#


# ================================================================================================================
# ================================================================================================================
#                                                    MAIN function
# ================================================================================================================
# ================================================================================================================

print("\nLaunch Multiple Simulations Process")
print("==========================================================\n")
print("substrate regime: " + substrate_regime + "\nstudy: " + study)

total_number_of_simu=len(sim)
print("Number of simulations to run:", total_number_of_simu)

count=0
for key, values in sim.items():
    count+=1
    print("\nSimulation number " + str(count) + " out of " + str(total_number_of_simu) + " (" + str(round(100.*count/total_number_of_simu, 2)) + "%)")
    print("------------------------------------------")

    
    ## if a data.txt exists in the target folder, the simulation is considered done and ignored
    path=getPath(values)
    if TEST:
        print(".. work folder:", path)

    if os.path.isfile(os.path.join(path, "data.txt")):
        print("...... simulation already done !")
        continue
    #
    
    print("...... preparation")
    editUI(values)

    ## run simulation
    if not TEST:
        os.system(r"start http://localhost:8000//" + path_to_model)
    print(f'...... simulation {key} is running')

    ## check when the simulation is done
    if not TEST:
        simu_is_done = os.path.exists(path_to_export)
        while (not simu_is_done):
            time.sleep(0.1)
            simu_is_done = os.path.exists(path_to_export)
    ##

    ## when the simulation is done, move the data file to target data analysis folder ...
    if not TEST:
        os.system(r'move "'+path_to_export + r'" "' + path)
    print('...... simulation is done, files are exported')

    ###... and close simulation
    if not TEST:
        os.system(r'taskkill /f /im firefox.exe')   # msedge.exe   firefox.exe
#


print("\nLaunch Multiple Simulations Process DONE")
print("==========================================================\n")
