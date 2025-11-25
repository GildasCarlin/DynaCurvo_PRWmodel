
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
// 'right_sinus': sinusoidal wavy surface  -- amplitude a (µm)
//                                         -- wavelength lambda (µm)
//                                         -- wave celerity c>=0 (µm/s)
//                                         -- initial phase phi0 (unitless)
//                                         -- angular velocity w (rad/s)
//                                         -- initial angle psi0 (rad)

let surface = 'right_sinus';

let a = 10.0;
let lambda = 20.0;
let c = 3.0;
let phi0 = 0.0;
let w = 1.0;
let psi0 = 0.0;
let angleParams=[[Math.PI/4.0, w, 600.0], 
                 [-Math.PI/4.0, w, 1200.0],
                ];




// II. Cell 
//-------------------------------------------------------------------------------------------------------

// number of cells
let nbCell = 300;

// central metric     median  /  mean
let central_metric = 'mean';

// assumption on cell dynamic    true  /  false
let surf_assumption = true;  
let curvature_assumption = true;

// III. Model integration
//-------------------------------------------------------------------------------------------------------

// trajectory sample (seconds)
let dt_mesure = 15.0;

// simulation is stopped when the window is closed. Tmax is the duration in seconds of the trajectories recording.
let Tmax = 1800.0;


// IV. Others
//-------------------------------------------------------------------------------------------------------

// export a .csv file containing trajectories data
let data = true;
let file_name_data = 'data';

// draw 
let draw_cells = true;
let draw_traj = false;   // deprecated if nbCell>5

// record screenshot video
let screen = false;
let captureLength = 2000;  // number of frames recorded






