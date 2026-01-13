
/**
 * File              : param_user.js
 * Original article  : doi
 * Authors           : Gildas Carlin, Ian Manifacier and Jean-Louis Milan
 * Date              : 05/05/2025
 * Description       : User interface.
*/
//


// I. Surface 
//-------------------------------------------------------------------------------------------------------
//
// 'flat' : flat homogeneous surface (static)
//                               -- index: 0
// 'right_sinus': sinusoidal wavy surface  
//                               -- index: 1
//                               -- amplitude a (µm)
//                               -- wavelength lambda (µm)
//                               -- wave celerity c>=0 (µm/s)
//                               -- initial phase phi0 (unitless)
//                               -- angular velocity w (rad/s)
//                               -- initial angle psi0 (rad)
// 'wells_sinus': sum of two perpendicular sinus (such as Pieuchot et al., 10.1038/s41467-018-06494-6)
//                               -- index: 2
//                               -- amplitude a (µm) (commun for the two directions)
//                               -- wavelength lambda (µm) (commun for the two directions)
//                               -- wave celerity c>=0 (µm/s) (only y-axis supports temporal shift)
//                               -- initial phase phi0 (unitless) (only y-axis supports initial shift)
//                               -- initial angle psi0!=0 (rad)
// 'convr_sinus': sum of two converging sinus (WiP)
//                               -- index: 3
//                               -- amplitude a (µm) (commun for the two directions)
//                               -- wavelength lambda (µm) (commun for the two directions)
//                               -- wave celerity c>=0 (µm/s)
//                               -- initial phase phi0 (unitless) (commun for the two directions)
//                               -- opening angle psi0!=0 (rad)
//                               -- rate of angle change w (rad/µm) 
//                               -- Dpsi>0 (rad)
// 'snake_sinus': 
//                               -- index: 4
//                               -- amplitude a (µm) (commun for the two directions)
//                               -- wavelength1 lambda (µm) 
//                               -- wavelength2 lambdas (µm)
//                               -- wave celerity c>=0 (µm/s)
//                               -- initial phase phi0 (unitless) (commun for the two directions)

let surface_title = 'convr_sinus';

let a = 10.0;
let lambda = 80.0;
let c = 1.0;
let phi0 = 0.0;

let w = -0.02;
let psi0 = 0.0-Math.PI/4.0;
let angleParams=[]; //[-Math.PI/4.0, w, 60.0], [Math.PI/4.0, w, 120.0]
let Dpsi=Math.PI/2.0;
let lambdas=80.0;




// II. Cell 
//-------------------------------------------------------------------------------------------------------

// number of cells
let nbCell=0;

// central metric     median  /  mean
let central_metric = 'mean';

// assumption on cell dynamic    true  /  false
let surf_assumption = true;  
let curv_assumption = true;

// III. Model integration
//-------------------------------------------------------------------------------------------------------

// trajectory sample (seconds)
let dt_mesure = 15.0;

// simulation is stopped when the window is closed. Tmax is the duration in seconds of the trajectories recording.
let Tmax = 1200.0;


// IV. Others
//-------------------------------------------------------------------------------------------------------

// export a .txt file containing trajectories data
let data = false;
let file_name_data = 'data';

// draw 
let draw_cells = true;
let draw_traj = false;   // deprecated if nbCell>5

// record screenshot video
let screen = false;
let captureLength = 2000;  // number of frames recorded






