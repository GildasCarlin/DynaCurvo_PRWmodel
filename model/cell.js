
/**
 * File              : cell.js
 * Original article  : doi
 * Authors           : Gildas Carlin, Ian Manifacier and Jean-Louis Milan
 * Date              : 05/05/2025
 * Description       : Cell class
 */
//


class Cell {
  constructor(set_of_parameters, central_vel, i) {
    this.cellID = i; 

    this.P2 = set_of_parameters[0];
    this.P1 = set_of_parameters[1];
    this.S2 = set_of_parameters[2];
    this.S1 = set_of_parameters[3];
    this.alpha = set_of_parameters[4];
    this.q = set_of_parameters[10];
    this.tau_stat = set_of_parameters[7];
    this.tau_dyna = set_of_parameters[8];

    this.Vmed = set_of_parameters[5];
    this.Vmean = set_of_parameters[6];
    this.Vmax = set_of_parameters[9];
    this.central_V_flat = central_vel;
    
    this.inv_polarityMatrix = [
      [1.0/this.P1, 0.0],
      [0.0, 1.0/this.P2], 
    ];

    this.diffusionMatrix = [ 
      [this.S1 / Math.sqrt(this.P1), 0.0],
      [0.0, this.S2 / Math.sqrt(this.P2)]
    ];

    this.beta = 0.0;
    //WiP
    this.f_mag=0.0;
    //this.f_ang=0.0;
    //this.f = createVector(Math.cos(this.f_ang), Math.sin(this.f_ang)).mult(this.f_mag);

    this.f_ext = createVector(0.0, 0.0);
    this.velocity = createVector(0.01*Math.random()-0.005, 0.01*Math.random()-0.005);
    //this.velocity = createVector(0.0, 0.0);
    this.position = createVector(400.0 * Math.random() - 200.0, 400.0 * Math.random() - 200.0);
    //this.position = createVector(0.0, 0.0);

    // cell shape (visual only)
    this.cell_radius_min = 10.0;  
    this.cell_radius_max = 15.0;

    let [gradH, gradH_perp] = gradCurvature(this.position, 0.0, angleParams);    // local basis
    this.k_angle=createVector(1.0, 0.0).angleBetween(gradH);

    //test
    this.V1=createVector(1.0, 1.0);
    this.V2=this.V1.copy();
  }
  
  // Run simulation
  // ------------------------------------------------------------------------------------
  run(t){

    // 1.trajectory recording
    this.sauvegarde_interval(t);

    // 1.update cell velocity and position 
    this.update(t);
 
    // 3.sketch
    this.render(this.position.copy(), this.velocity.copy());
  }

  // Update cell velocity and position
  // ------------------------------------------------------------------------------------
  update(t){

    // Polarity matrix and Diffusion matrix are updated
    this.Inferred_parameters_mapping(this.position, t);
    
    //-------------------------------------------- external forces caused by curvature and dynamic wave
    this.f_ext = createVector(0.0, 0.0);
    if (curvature_assumption) this.f_ext.add(this.get_Psi_curvature(this.position.copy(), this.velocity.copy(), t));
    if (surf_assumption) this.f_ext.add(this.get_Psi_wave(this.position.copy(), t));

    //let test=this.f_ext.copy();
    //test=this.From_local_to_global(test, t);
    //drawVector(this.position.x, this.position.y, test.copy().mult(10000), 'red');
    //-------------------------------------------- velocity computation
    
    let DfDr_dot_v=this.MultiplyMatrixWVector(this.get_Jacobian_Matrix(this.position.copy(), t), this.From_global_to_local(this.velocity.copy(), t));
    DfDr_dot_v=this.From_local_to_global(DfDr_dot_v,t);
    //drawVector(this.position.x, this.position.y, DfDr_dot_v.copy().mult(10000.0), 'red');

    let Pers_dot_pol=this.MultiplyMatrixWVector(this.inv_polarityMatrix, this.From_global_to_local(this.velocity.copy(), t).sub(this.f_ext.copy()));
    Pers_dot_pol=this.From_local_to_global(Pers_dot_pol,t);

    this.velocity = this.velocity.add(
      DfDr_dot_v.mult(0.0*dt*this.q)
    ).sub(
      Pers_dot_pol.mult(dt)
    ).add(
      this.MultiplyMatrixWVector(this.diffusionMatrix, this.get_random_vector()).mult(Math.sqrt(dt))
    );
    //this.velocity=this.From_local_to_global(this.velocity, t);
    //this.velocity = this.velocity.add(
    //  this.MultiplyMatrixWVector(this.get_Jacobian_Matrix(this.position.copy(), t), this.velocity.copy()).mult(dt*this.q)
    //).sub(
    //  this.MultiplyMatrixWVector(this.inv_polarityMatrix, this.velocity.copy().sub(this.f_ext.copy())).mult(dt)
    //).add(
    //  this.MultiplyMatrixWVector(this.diffusionMatrix, this.get_random_vector()).mult(Math.sqrt(dt))
    //);

    //-------------------------------------------- velocity integration to get cell position
    this.position.add(this.velocity.copy().mult(dt));
  }
  
  // Intermediate functions
  // ------------------------------------------------------------------------------------
  
  get_random_vector() {   
    // Box-Muller transform 
    let u = 0.0, v = 0.0; 
    while(u === 0.0) u = Math.random();  // Converting [0,1) to (0,1)
    while(v === 0.0) v = Math.random();
    let noise = createVector(
      Math.sqrt(-2.0*Math.log(u)) * Math.cos(2.0*Math.PI*v),
      Math.sqrt(-2.0*Math.log(u)) * Math.sin(2.0*Math.PI*v)
    );
    return noise;
  }

  get_Psi_curvature(pos, vel, t){
    // this function returns the contribution of f_ext caused by the static curvature. The formulation is discontinuous for the angle

    vel=this.From_global_to_local(vel.copy(),t);
    // 1. mesurements of the angle between e2 and the velocity
    let theta=0.0;
    if (meanCurvature(pos, t, angleParams)>=0.0){
      // normalized local basis (e1, e2) 
      let [e1, e2]=Normalize_Vec(this.From_global_to_local(gradCurvature(pos, t, angleParams), t));
      //let [e1, e2]=Normalize_Vec(gradCurvature(pos, t, angleParams));

      // is the cell polarized along +e2 (True) or -e2 (False) ?
      if (vel.dot(e2)<0.0){
        e2.mult(-1.0);
      } 
      theta=e2.angleBetween(vel);
    }
    
    // 2. define the orientation of f_stat
    let w=Normalize_Vec(vel.copy().rotate(-this.alpha*theta));

    // 3. norm of psi
    let f_d = 0.24*Math.exp(-lambda/19.8) + 0.02;               
    this.f_mag -= dt/this.tau_stat * (this.f_mag-f_d);     
    
    return w.mult(this.f_mag);
  }

  getPsiStatic(pos, pol, t){
    // not finished
    // this function returns the contribution of f_ext caused by the static curvature. The formulation is continuous
    if (meanCurvature(pos, t, angleParams)<0.0){
      this.f=createVector(0.0, 0.0);
    } else {                                
      let [e1, e2] = gradCurvature(pos, t, angleParams);
      let f_d = Math.abs(this.Vmean - this.central_V_flat);

      // is the cell polarized along +e2 (True) or -e2 (False) ?
      if (pol.dot(e2)<0.0){
        e2.mult(-1.0);
      }
      let equilibrium=Normalize_Vec(e2).mult(f_d);

      this.f_mag -= dt/this.tau_stat * (this.f_mag-f_d);          
      this.f_ang -= dt/this.tau_stat * (this.f_ang-createVector(1.0, 0.0).angleBetween(equilibrium)); 
      this.f=createVector(Math.cos(this.f_ang), Math.sin(this.f_ang)).mult(this.f_mag);
   }
   return this.f;
  }

  get_Psi_wave(pos, t){
    let W = this.From_global_to_local(phase_velocity(pos, t, angleParams), t);
    //let W = phase_velocity(pos, t, angleParams);
    let Wspeed = W.mag();

    // target velocity
    let s = Wspeed / (1.0+(Wspeed/this.Vmax)**3.0);

    this.beta -= dt/this.tau_dyna * (this.beta-s);       // update beta
    return Normalize_Vec(W).mult(this.beta);
  }

  get_Jacobian_Matrix(pos, t){
    // local basis
    let [gradH, gradH_perp]=gradCurvature(pos, t, angleParams);
    // global basis
    //let [X, Y]=[createVector(1.0, 0.0), createVector(0.0, 1.0)];

    // base locale
    let sgn=Math.sign(gradH.dot(this.velocity));
    let Mxx=gradH.mag() * sgn;
    let Myy=0.0;
    let Mxy=0.0;

    //let theta=gradH.angleBetween(X);
    //let sgn=Math.sign(gradH.dot(this.velocity));
    ////let sgn=Math.sign(gradH.y);
    //let curv=gradH.mag();
    //let Mxx=sgn*curv*Math.cos(theta)**2.0;
    //let Myy=sgn*curv*Math.sin(theta)**2.0;
    //let Mxy=sgn*Math.abs(curv)*Math.cos(theta)*Math.sin(theta);
  
    //let Mxx=gradH.x;
    //let Myy=gradH.y;
    //let Mxy=0.0;

    //let e1=Normalize_Vec(gradH.copy());
    //let Mxx=-gradH.x * e1.y;
    //let Myy=gradH.y * e1.x;
    //let Mxy=-gradH.y * e1.y;
    //let Myx=gradH.x * e1.x;

    //let Mxx=gradH.x*Math.sin(theta)**2.0 + gradH.y*Math.cos(theta)**2.0;
    //let Myy=gradH.x*Math.cos(theta)**2.0 + gradH.y*Math.sin(theta)**2.0;
    //let Mxy=Math.abs(gradH.x-gradH.y)*Math.cos(theta)*gradH.y*Math.sin(theta);

    return [[Mxx, Mxy], [Mxy, Myy]];
    //return [[Mxy, Myy], [-Mxx, -Mxy]];
  }

  Inferred_parameters_mapping(pos, t){
    // local basis expressed in global basis
    let [e1, e2]=Normalize_Vec(gradCurvature(pos, t, angleParams));

    // U : matrice de changement de base de local à global
    // Soit D une matrice 2x2 diagonale exprimée en base locale, alors D en base globale est donné par,
    // D(globale) = U.D(locale).U^T
    let U=[[e1.x, e2.x],
           [e1.y, e2.y]]

    // persistence matrix (from local basis P1, P2 to global basis Pxx, Pxy, Pyy)
    let Pxx=1.0/this.P1*(U[0][0])**2.0 + 1.0/this.P2*(U[0][1])**2.0;
    let Pyy=1.0/this.P1*(U[1][0])**2.0 + 1.0/this.P2*(U[1][1])**2.0;
    let Pxy=1.0/this.P1*U[0][0]*U[1][0] + 1.0/this.P2*U[0][1]*U[1][1];
    this.inv_polarityMatrix = [
      [Pxx, Pxy],
      [Pxy, Pyy], 
    ];

    // diffusion matrix
    let Sxx=this.S1/Math.sqrt(this.P1)*U[0][0]**2.0 + this.S2/Math.sqrt(this.P2)*U[0][1]**2.0;
    let Syy=this.S1/Math.sqrt(this.P1)*U[1][0]**2.0 + this.S2/Math.sqrt(this.P2)*U[1][1]**2.0;
    let Sxy=this.S1/Math.sqrt(this.P1)*U[0][0]*U[1][0] + this.S2/Math.sqrt(this.P2)*U[0][1]*U[1][1];
    this.inv_polarityMatrix = [
      [Sxx, Sxy],
      [Sxy, Syy], 
    ];

  }

  MultiplyMatrixWVector(A, v){
    let v1=v.x, v2=v.y;
    return createVector(A[0][0]*v1 + A[0][1]*v2, A[1][0]*v1 + A[1][1]*v2);
  }

  DoubleDotProduct(matrixA, matrixB) {
  let sum = 0.;
  for (let i = 0; i < 2; i++) {
    for (let j = 0; j < 2; j++) {
      sum += matrixA[i][j] * matrixB[i][j];
    }
  }
  return sum;
}

  Heaviside_step_function(s){
    if (s<=0.0){
      return 0.0;
    } else {
      return 1.0;
    }
  }

  From_global_to_local(vecs, t){
    // take as input: a vector or an array of vectors expressed in the global static base (x, y)
    // returns: the same vector(s) expressed in the local basis (e1, e2) 

    const isArray = Array.isArray(vecs);
    const array_of_vectors = isArray ? vecs : [vecs];  // if vecs is a single vector in argument, it is stored in an array
    
    let [e1, e2]=gradCurvature(this.position, t, angleParams);
    let theta=createVector(1.0, 0.0).angleBetween(e1);
    let c = Math.cos(theta), s = Math.sin(theta);
    
    let rotation_done=[];
    for (let i=0; i<array_of_vectors.length; i++){
      let vx=array_of_vectors[i].x, vy=array_of_vectors[i].y;
      rotation_done.push(createVector(c*vx + s*vy, -s*vx + c*vy));
    }

    return isArray ? rotation_done : rotation_done[0];
  }

  From_local_to_global(vecs, t){
    // take as input: a vector or an array of vectors expressed in the local basis (e1, e2)
    // returns: the same vector(s) expressed in the global static base (x, y) 

    const isArray = Array.isArray(vecs);
    const array_of_vectors = isArray ? vecs : [vecs];  // if vecs is a single vector in argument, it is stored in an array
    
    let [e1, e2]=gradCurvature(this.position, t, angleParams);
    let theta=createVector(1.0, 0.0).angleBetween(e1);
    let c = Math.cos(theta), s = Math.sin(theta);
    
    let rotation_done=[];
    for (let i=0; i<array_of_vectors.length; i++){
      let vx=array_of_vectors[i].x, vy=array_of_vectors[i].y;
      rotation_done.push(createVector(c*vx - s*vy, s*vx + c*vy));
    }

    return isArray ? rotation_done : rotation_done[0];
  }
 

  // 3.sketch
  render(pos, vel){
    if (draw_cells && this.isInBorders(pos)){
      this.draw_cell(pos, vel);
    }
  }

  isInBorders(pos) { 
  //is the cell inside the window?
    if (pos.x < -width/2 || pos.y < -height/2 || pos.x > width/2 || pos.y > height/2){
      return false;
    } else {
      return true;
    }
  }

  draw_cell(pos, vel){
    let x = pos.x, y = pos.y;
    
    // cytoplasm
    push();
    fill(215, 95, 76);
    stroke(0);
    translate(x, y);
    rotate(createVector(1.0, 0.0).angleBetween(vel));

    let a = this.cell_radius_min + this.cell_radius_max*(1.0-exp(-vel.mag()/(0.5*(this.S1+this.S2))));
    let b = 0.5*a;
    if (b<this.cell_radius_min){
      b = this.cell_radius_min;
    }
    ellipse(0, 0, 2.0*a, 2.0*b);
    pop();
    
    // nucleus is drawn only if there are not too much cells
    if (nbCell<50){
      push();
      fill(242, 142, 19);
      stroke(0);
      translate(x, y);
      ellipse(0, 0, this.cell_radius_min, this.cell_radius_min);
      pop();
    }
    
    // direction of migration is drawn only if there are not too much cells
    if (nbCell<50){
      drawVector(x, y, vel.copy().mult(25.0), color(1, 3, 38));
    }
  }
  
  // 1.trajectory recording
  sauvegarde_interval(t){
    if (t_last_mesure>=it_mesure){
      this.sauvegarde(t);
      if (this.cellID == nbCell){ 
        t_last_mesure = 0;
        it_sauv += 1;
      }
    }
  }
  
  sauvegarde(t){
    table_data.push(`${it_sauv}, ${t}, ${this.cellID}, ${this.position.x}, ${this.position.y}, ${this.k_angle}`);
  }
  
}
