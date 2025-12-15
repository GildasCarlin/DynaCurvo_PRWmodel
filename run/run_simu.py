# ================================================================================================================
# ================================================================================================================
#                                                    packages importation
# ================================================================================================================
# ================================================================================================================

import os
import time
import regex as re


# ================================================================================================================
# ================================================================================================================
#                                                    intermediate functions
# ================================================================================================================
# ================================================================================================================

def search_files_to_analyse(root):
    paths_to_files = []

    for root_folder, sub_folder, files in os.walk(root):
        if re.search('c=', root_folder) and not ('data.csv' in files or 'data.txt' in files) :
            paths_to_files.append(os.path.join(root_folder))
    return paths_to_files
#

# ================================================================================================================
# ================================================================================================================
#                                                    User interface
# ================================================================================================================
# ================================================================================================================


## ------------------------------------ folder and paths

# study
substrate_regime="dyna"  # dyna or static
study="checkCache"

# exported file (at the end of the simulation)
home=os.path.expanduser("~")
path_to_export=os.path.join(home, "Downloads", "data.txt")

# model version + path to template param_user.js
path_to_model = os.path.join("model") 
path_to_UI = os.path.join("..", "model", "param_user.js")

# assets
path_to_assets=os.path.join("..", "post", substrate_regime, "assets", study)


## ------------------------------------ simulations parameters

# assumption    cell_is_surfer (true) /  cell_is_not_surfer (false)
W_assumption = r"true"
K_assumption = r"true"


# ================================================================================================================
# ================================================================================================================
#                                                    MAIN function
# ================================================================================================================
# ================================================================================================================

print("\nLaunch Multiple Simulations Process")
print("==========================================================\n")
print("substrate regime: " + substrate_regime + "\nstudy: " + study)

path_to_simu = search_files_to_analyse(path_to_assets)
total_number_of_simu=len(path_to_simu)
print("Number of simulations to run:", total_number_of_simu)

count=0
for path in path_to_simu:
    count+=1
    print("\nSimulation number " + str(count) + " out of " + str(total_number_of_simu))
    print("------------------------------------------")

    print("...... preparation")
    # parameters to edit
    #(wl, c) = re.findall(r'lambda=(\d+\.\d+|\d+)\\c=(\d+\.\d+|\d+)', path)[0]
    c=re.findall(r"c=(\d+\.\d+|\d+)", path)[0]
    wl=20.0
    
    print('... parameters  {}µm  {}µm/s'.format(wl, c))
    # open UI before editing parameters
    temp_file = open(path_to_UI)
    temp_content = temp_file.read()
    temp_file.close()

    # edit parameters in UI
    temp_content = re.sub(r'let\sc\s*=\s*(\d*\.?\d*);', r'let c = '+str(c)+r';', temp_content)
    temp_content = re.sub(r'let\slambda\s*=\s*(\d*\.?\d*);', r'let lambda = '+str(wl)+r';', temp_content)
    temp_content = re.sub(r'let\ssurf_assumption\s=\s(true|false);', r'let surf_assumption = '+W_assumption+r';', temp_content) 
    temp_content = re.sub(r'let\scurvature_assumption\s=\s(true|false);', r'let curvature_assumption = '+K_assumption+r';', temp_content)

    # save new UI
    UI_to_write = open(path_to_UI, 'w')
    UI_to_write.write(temp_content)
    UI_to_write.close()

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
