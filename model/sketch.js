/**
 * File              : sketch.js
 * Original article  : doi
 * Authors           : Gildas Carlin, Ian Manifacier and Jean-Louis Milan
 * Date              : 05/05/2025
 * Description       : Main file performing simulation and simulation sketch.
 */
//

// Variables initialization
//-------------------------------------------------------------------------------------------------------

// I. Surface 
//-------------------------------------------------------------------------------------------------------

// a surface is defined by an explicit string 'flat' or 'right_sinus' or by an index 0 or 1 respectively
let surface_index = 0;
if (surface=='right_sinus'){
  surface_index = 1;
}  

// wave number
let k = 2.0*Math.PI/lambda;
let kc = k * c;

// minimal and maximal curvature
Kmax = a*k**2.0;
Kmin = -1.0*Kmax;


// II. Cell 
//-------------------------------------------------------------------------------------------------------

// trajectories are stored in a table
let cellID = 0;
let cells = [];

// persistence and speed parameters are stored in a table
let all_parameters;
let parameters;


// III. Model integration 
//-------------------------------------------------------------------------------------------------------

// integrative step time
let dt = dt_mesure/100.0; 
let it_mesure = dt_mesure/dt;

// total number of iterations
let iTmax = Tmax/dt;

// chronometer init
let t = 0;

// trajectories sampling parameters
let it_sauv=0;
let t_last_mesure = 0;


// IV. Others
//-------------------------------------------------------------------------------------------------------

// drawing trajectories
let current;
let previous;
let painting = false;
let paths = [];
let next = 0.0;

// background sketch
let theShader;

// font
let inconsolata;

// table to export init
let row = []; 
let table_data;

// screenshot video
let capturer = new CCapture({
  frameRate: 60,
  format: "webm",
});


// Functions definition
//-------------------------------------------------------------------------------------------------------

function preload(){
  // persistence and speed parameters init
  all_parameters = loadTable('assets/parameters.csv', 'csv', 'header');
  
  inconsolata = loadFont('assets/Inconsolata.otf');
  
  theShader = loadShader('shader.vert', 'shader.frag'); 
  current = createVector(0,0);
  previous = createVector(0,0); 
}

// sketch window and cells initialization
function setup() {
  createCanvas(800, 800, WEBGL);
  noStroke();
  shaderGraphics = createGraphics(width, height, WEBGL);
  shaderGraphics.noStroke();

  // select the correct set of parameters
  parameters = all_parameters.getColumn('flat');

  let V_flat = parseFloat(parameters[6]);  // the mean by default
  if (central_metric=='median'){
    V_flat = parseFloat(parameters[5]);
  }


  if (surface=='right_sinus'){
    parameters = all_parameters.getColumn(Math.floor(lambda).toString());
  }
  parameters = parameters.map(parseFloat);


  table_data = initTable();
  for (let i = 0; i < nbCell; i++) {
    cells[i] = new Cell(parameters, V_flat, i+1);
    cells[i].sauvegarde(t);
    paths.push([]);
  }
  it_sauv += 1;
}

function drawScaleBar(font) {
  push();
  stroke(0)
  fill('black');
  rect(-0.47*width, -300.0, 115.0, 50.0);
  textFont(font, 22);
  fill('white');
  text('80µm', -0.428*width, -280.0, 10.0);
  fill('white');
  rect(-0.445*width,-270.0,80,10.0);
  pop();
}

function drawTimer(font, t) {
  let minutes = Math.floor(t/60.0).toString();
  let seconds = (t % 60.0).toFixed(0).toString();

  push();
  fill('black');
  textFont(font, 22);
  rect(-0.47*width, -0.425*height,115.0,27);
  fill('white');
  text(minutes+'min '+seconds+'s', -0.45*width, -0.4*height, 800);
  pop();
}

function initBackground(t){
  background(240);
  
  push(); 
  shaderGraphics.shader(theShader);

  // Set the resolution as a uniform
  theShader.setUniform('u_resolution', [width, height]);
  theShader.setUniform('surface', surface_index);
  theShader.setUniform('k', k);
  theShader.setUniform('kc', kc);
  theShader.setUniform('phi0', phi0);
  theShader.setUniform('psi0', psi0);
  theShader.setUniform('t', t);
  theShader.setUniform('s', height);
  theShader.setUniform('num_angles', angleParams.length);
  theShader.setUniform('temporal_angle', angleParams.flat());

  shaderGraphics.rect(0.0, 0.0, height, width);
  image(shaderGraphics, -width/2.0, -height/2.0, width, height);

  // center of the window
  fill('grey');
  ellipse(0.0, 0.0, 10.0, 10.0);
  pop();
}

function initTable(){
  let table_data = [];

  table_data.push('# Metadata');
  table_data.push(`# Surface=${surface}`);
  table_data.push(`# Surface Index=${surface_index}`);
  table_data.push(`# Amplitude=${a}`);
  table_data.push(`# Wavelength=${lambda}`);
  table_data.push(`# Wave Celerity=${c}`);
  table_data.push(`# Angular Velocity=${w}`);
  table_data.push(`# Initial Phase=${phi0}`);
  table_data.push(`# Initial Angle=${psi0}`);
  table_data.push(`# Number of Cells=${nbCell}`);
  table_data.push(`# Simulation Duration=${Tmax}`);
  table_data.push(`# Time Interval=${dt_mesure}`);
  table_data.push(`# Px=${parameters[0]}`);
  table_data.push(`# Py=${parameters[1]}`);
  table_data.push(`# Sx=${parameters[2]}`);
  table_data.push(`# Sy=${parameters[3]}`);
  table_data.push(`# alpha=${parameters[4]}`);
  table_data.push(`# tau_stat=${parameters[7]}`);
  table_data.push(`# tau_dyna=${parameters[8]}`);
  table_data.push(`# q=${parameters[10]}`);
  table_data.push(`# Vmed=${parameters[5]}`);
  table_data.push(`# Vmean=${parameters[6]}`);
  table_data.push(`# Vmax=${parameters[9]}`);
  table_data.push(`# surf_assumption=${surf_assumption}`);
  table_data.push(`# curvature_assumption=${curvature_assumption}`);
  table_data.push('frame,time,cell,x,y, theta');
  
  return table_data;
}


// MAIN
//-------------------------------------------------------------------------------------------------------

function draw(){

  // update time in seconds
  let timer_sec = t*dt;

  //start record video
  if (screen){
    if (frameCount == 1) {
      capturer.start();
    }
  }

  // background sketch 
  initBackground(timer_sec);

  // sketch cells trajectories. Deprecated if nbCells>5
  if (draw_traj){
    trajectories();
  }
  
  // cell update
  for (let i = 0; i<cells.length; i++){
    cells[i].run(timer_sec);
  }

  t += 1;
  t_last_mesure += 1;

  // sketch of the scale bar and the chronometer
  drawScaleBar(inconsolata);
  drawTimer(inconsolata, timer_sec);

  //stop recording video
  if (screen){ 
    if (frameCount < captureLength) {
      capturer.capture(canvas);
    } else if (frameCount === captureLength) {
      capturer.save();
      capturer.stop();
    }
  }

  // data 
  export_the_data();
}



















