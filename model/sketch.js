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
let surface_repertory={"flat":0,
                       "right_sinus":1,
                       "wells_sinus":2,
                       "convr_sinus":3,
                       "snake_sinus":4
                      };

let surface_representations=[surface_title, surface_repertory[surface_title]];

// wave number
let k = 2.0*Math.PI/lambda;
let kc = k*c;

let ks = 2.0*Math.PI/lambdas;
let ksc = ks*c;


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

  if ([1, 2].includes(surface_representations[1])){
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

function drawSpeedometer(font, val) {
  if (nbCell==1){
    push();
    fill('black');
    textFont(font, 22);
    rect(-0.47*width, -0.3*height,115.0,60);
    fill('white');
    text("Cell Speed", -0.466*width, -0.275*height, 800);
    text(val.toFixed(2).toString(), -0.4257*width, -0.24*height, 800); // µm/s
    pop();
  }
}

function initBackground(t){
  background(240);
  
  push(); 
  shaderGraphics.shader(theShader);

  // Set the resolution as a uniform
  theShader.setUniform('u_resolution', [width, height]);
  theShader.setUniform('surface', surface_representations[1]);
  theShader.setUniform('as', as);  
  theShader.setUniform('k', k);
  theShader.setUniform('kc', kc);
  theShader.setUniform('ks', ks);
  theShader.setUniform('ksc', ksc);
  theShader.setUniform('phi0', phi0);
  theShader.setUniform('psi0', psi0);
  theShader.setUniform('w', w);
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
  table_data.push(`# Surface=${surface_representations[0]}`);
  table_data.push(`# Surface Index=${surface_representations[1]}`);
  table_data.push(`# Amplitude=${a}`);
  table_data.push(`# Amplitude (snake)=${as}`);
  table_data.push(`# Wavelength=${lambda}`);
  table_data.push(`# Wavelength (snake)=${lambdas}`);
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
  table_data.push(`# curvature_assumption=${curv_assumption}`);
  table_data.push('frame,time,cell,x,y,z,theta,v_dot_e1,meanCurvature');
  
  return table_data;
}


// MAIN
//-------------------------------------------------------------------------------------------------------
let max=0.0;
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
  if (nbCell==1){
    drawSpeedometer(inconsolata, sqrt(MetricDotProduct(cells[0].position, timer_sec, angleParams, cells[0].velocity, cells[0].velocity)));
  }
  

  //stop recording video
  if (screen){ 
    if (frameCount < captureLength) {
      capturer.capture(canvas);
    } else if (frameCount === captureLength) {
      capturer.save();
      capturer.stop();
    }
  }

  let X=[];// [0.0, -50.0, -110.0, 110.0, 50.0, -300.0, 300.0];
  let Y=[0.0];
  console.log("-------------------");
  for (const x of X){
    for (const y of Y){
      let posiF=createVector(x, y);
      let H=meanCurvature(posiF, timer_sec, angleParams);
      let K=gaussianCurvature(posiF, timer_sec, angleParams);

      push();
      fill('black');
      stroke(0);
      translate(x, y);
      ellipse(0, 0, 2.0, 2.0);
      pop();

      let [gradf, gradpf]=gradCurvature(posiF, timer_sec, angleParams);
      let gradGf=gradGaussianCurvature(posiF, timer_sec, angleParams);
      drawVector(x, y, gradf.copy().mult(15000), 'red');
      //drawVector(x, y, gradGf.copy().mult(700000), 'green');

      let [k1, k2, v1, v2]=PrincipalCurvature(posiF, timer_sec, angleParams);
      let [kxx, kyy, kxy] = curvatureTensor(posiF, timer_sec, angleParams)
      
      console.log("height=", substrate(posiF, timer_sec, angleParams).toFixed(5));
      //console.log("H, K=", (H**2.0).toFixed(5), K.toFixed(5));
      //console.log("kxx, kyy, kxy", kxx.toFixed(5), kyy.toFixed(5), kxy.toFixed(5))
      //console.log("k1=", k1.toFixed(5));
      //console.log("k2=", k2.toFixed(5));
      //drawVector(x, y, v1.copy().mult(25), 'red');
      //drawVector(x, y, v2.copy().mult(20), 'green');
    }
  }

  //let xm=-200.0;//+Math.cos(psi0)*c*timer_sec;  //x*Math.sin(psi) + y*Math.cos(psi)
  //let ym=0.0-c*timer_sec;
  //let posiM=createVector(xm, ym);
  //
  //let xf=60.0, yf=0.0;
  //let posiF=createVector(xf, yf);



  //push();
  //fill('red');
  //stroke(0);
  //translate(xm, ym);
  //ellipse(0, 0, 10.0, 10.0);
  //pop();



  //let [gradf, gradpf]=gradCurvature(posiF, timer_sec, angleParams);
  //let [gradm, gradpm]=gradCurvature(posiM, timer_sec, angleParams);
  //let gradGf=gradGaussianCurvature(posiF, timer_sec, angleParams);
  //drawVector(xf, yf, gradf.copy().mult(15000), 'red');
  //drawVector(xf, yf, gradGf.copy().mult(1000000), 'green');
  //drawVector(xm, ym, gradm.copy().mult(10000000), 'red');
  //let speedM=phaseVelocity(posiM, timer_sec, angleParams);
  //let speedF=phaseVelocity(posiF, timer_sec, angleParams);
  //drawVector(xf, yf, speedF.copy().mult(100), 'red');
  //drawVector(xm, ym, speedM.copy().mult(100), 'red');


  //console.log("H:", H.toFixed(5));
  //console.log("K:", K.toFixed(5));
  //console.log("H max:", (0.5*a*k**2.0).toFixed(5));
  //console.log("psi :", psi0*180.0/Math.PI);
  //console.log("cos, sin :", Math.cos(psi0).toFixed(3), Math.sin(psi0).toFixed(3));
  //console.log("heigth fixed point :", substrate(posiF, timer_sec, angleParams).toFixed(3));
  //console.log("heigth moving point :", substrate(posiM, timer_sec, angleParams).toFixed(4));
  //console.log("gradient selon x :", gradf.x.toFixed(5));
  //console.log("gradient selon y :", gradf.y.toFixed(5));
  //console.log("speed moving point :", speedM.mag().toFixed(4));
  //console.log("speed fixed point :", speedF.mag().toFixed(4));

  // data 
  export_the_data();
}



















