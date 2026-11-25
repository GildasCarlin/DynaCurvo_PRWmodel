"""
 * File              : param_user.js
 * Original article  : doi
 * Authors           : Gildas Carlin, Ian Manifacier and Jean-Louis Milan
 * Date              : 05/05/2025
 * Description       : User interface.
"""


## I. Surface 
##-------------------------------------------------------------------------------------------------------
##
## 'flat' : flat homogeneous surface (static)
##                               -- index: 0
## 'right_sinus': sinusoidal wavy surface  
##                               -- index: 1
##                               -- amplitude a (µm)
##                               -- wavelength lambda (µm)
##                               -- wave celerity c>=0 (µm/s)
##                               -- initial phase phi0 (unitless)
##                               -- angular velocity w (rad/s)
##                               -- initial angle psi0 (rad)
## 'wells_sinus': sum of two perpendicular sinus (such as Pieuchot et al., 10.1038/s41467-018-06494-6)
##                               -- index: 2
##                               -- amplitude a (µm) (commun for the two directions)
##                               -- wavelength lambda (µm) (commun for the two directions)
##                               -- wave celerity c>=0 (µm/s) (only y-axis supports temporal shift)
##                               -- initial phase phi0 (unitless) (only y-axis supports initial shift)
##                               -- initial angle psi0!=0 (rad)
## 'convr_sinus': sum of two converging sinus (WiP)
##                               -- index: 3
##                               -- amplitude a (µm) (commun for the two directions)
##                               -- wavelength lambda (µm) (commun for the two directions)
##                               -- wave celerity c>=0 (µm/s)
##                               -- initial phase phi0 (unitless) (commun for the two directions)
##                               -- rate of angle change w (1/µm^2) 
##                               -- window of the curve (µm)
##                               
## 'snake_sinus': 
##                               -- index: 4
##                               -- amplitude1 a (µm) 
##                               -- amplitude2 as (µm)
##                               -- wavelength1 lambda (µm) 
##                               -- wavelength2 lambdas (µm)
##                               -- wave celerity c>=0 (µm/s)
##                               -- initial phase phi0 (unitless) (commun for the two directions)


sim={1:{"surface_title":r"'right_sinus'",
        "a":10.0,
        "lambda":20.0,
        "c":0.1,
        "phi0":0.0,
        "w":0.0,
        "psi0":r"-Math.PI/2.0",
        "nbCell":150,
        "surf_assumption":r"false",
        "curv_assumption":r"true",
        "dt_mesure":0.1,
        "Tmax":600.0,
        "data":r"true",
        "draw_cells":r"false"},
     2:{"surface_title":r"'right_sinus'",
        "a":10.0,
        "lambda":20.0,
        "c":0.4,
        "phi0":0.0,
        "w":0.0,
        "psi0":r"-Math.PI/2.0",
        "nbCell":150,
        "surf_assumption":r"false",
        "curv_assumption":r"true",
        "dt_mesure":0.1,
        "Tmax":600.0,
        "data":r"true",
        "draw_cells":r"false"},
     3:{"surface_title":r"'right_sinus'",
        "a":10.0,
        "lambda":20.0,
        "c":0.6,
        "phi0":0.0,
        "w":0.0,
        "psi0":r"-Math.PI/2.0",
        "nbCell":150,
        "surf_assumption":r"false",
        "curv_assumption":r"true",
        "dt_mesure":0.1,
        "Tmax":600.0,
        "data":r"true",
        "draw_cells":r"false"},
        }




