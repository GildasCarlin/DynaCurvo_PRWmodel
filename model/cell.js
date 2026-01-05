
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
    
    this.Local_Persistence_Matrix = [
      [1.0/this.P1, 0.0],
      [0.0, 1.0/this.P2], 
    ];

    this.Local_Diffusion_Matrix = [ 
      [this.S1 / Math.sqrt(this.P1), 0.0],
      [0.0, this.S2 / Math.sqrt(this.P2)]
    ];

    this.Global_Persistence_Matrix = [
      [0.0, 0.0],
      [0.0, 0.0], 
    ];

    this.Global_Diffusion_Matrix = [ 
      [0.0, 0.0],
      [0.0, 0.0]
    ];

    this.beta = 0.0;
    //WiP
    this.f_mag=0.0;
    //this.f_ang=0.0;
    //this.f = createVector(Math.cos(this.f_ang), Math.sin(this.f_ang)).mult(this.f_mag);

    this.f_ext = createVector(0.0, 0.0);
    //this.velocity = createVector(0.01*Math.random()-0.005, 0.01*Math.random()-0.005);
    this.velocity = createVector(0.0, 0.0);
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

    //-------------------------------------------- Polarity matrix and Diffusion matrix are updated in the global basis
    this.InferredParametersGlobalMapping(this.position, t);
    
    //-------------------------------------------- external forces caused by curvature and dynamic wave
    this.f_ext = createVector(0.0, 0.0);
    if (curv_assumption) this.f_ext.add(this.getPsiCurv(this.position, this.velocity.copy(), t));
    if (surf_assumption) this.f_ext.add(this.getPsiWave(this.position, t));

    //-------------------------------------------- velocity computation
    this.velocity = this.velocity.add(
      MatrixVectorProduct(this.get_Jacobian_Matrix(this.position, t), this.velocity).mult(dt*this.q)
    ).sub(
      MatrixVectorProduct(this.Global_Persistence_Matrix, this.velocity.copy().sub(this.f_ext)).mult(dt)
    ).add(
      MatrixVectorProduct(this.Global_Diffusion_Matrix, getRandomVector()).mult(Math.sqrt(dt))
    );

    //-------------------------------------------- velocity integration to get cell position
    this.position.add(this.velocity.copy().mult(dt));
  }
  
  // Intermediate functions
  // ------------------------------------------------------------------------------------
  
  getPsiCurv(pos, vel, t){
    // this function returns the contribution of f_ext caused by the static curvature. The formulation is discontinuous for the angle

    vel=NormalizeVector(this.From_global_to_local(vel, t));
    drawVector(0.0, 0.0, this.From_local_to_global(vel.copy(), t).mult(15.0), 'red');
    
    // 1. mesurements of the angle between e2 and the velocity
    let theta=0.0;
    if (meanCurvature(pos, t, angleParams)>=0.0){
      // local basis (e1, e2) 
      let [e1, e2]=NormalizeVector(this.From_global_to_local(gradCurvature(pos, t, angleParams), t));
      drawVector(0.0, 0.0, this.From_local_to_global(e1, t).copy().mult(15.0), 'black');
      drawVector(0.0, 0.0, this.From_local_to_global(e2, t).copy().mult(15.0), 'blue');
      //let [e1, e2]=Normalize_Vec(gradCurvature(pos, t, angleParams));

      // is the cell polarized along +e2 (True) or -e2 (False) ?
      if (vel.dot(e2)<0.0){
        e2.mult(-1.0);
      } 
      theta=e2.angleBetween(vel);
    }
    console.log("--------- debut du print ---------")
    console.log("theta :", (theta*180.0/Math.PI).toFixed(5));
    console.log("sin, cos :", Math.sin(theta).toFixed(5), Math.cos(theta).toFixed(5));
    console.log("vx, vy :", vel.x.toFixed(5), vel.y.toFixed(5));
    console.log("atan2", (atan2(vel.x, vel.y)*180.0/Math.PI).toFixed(5))
    
    // 2. define the orientation of f_stat
    let w=NormalizeVector(vel.copy().rotate(-this.alpha*theta));

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
      let equilibrium=NormalizeVector(e2).mult(f_d);

      this.f_mag -= dt/this.tau_stat * (this.f_mag-f_d);          
      this.f_ang -= dt/this.tau_stat * (this.f_ang-createVector(1.0, 0.0).angleBetween(equilibrium)); 
      this.f=createVector(Math.cos(this.f_ang), Math.sin(this.f_ang)).mult(this.f_mag);
   }
   return this.f;
  }

  getPsiWave(pos, t){
    let W = this.From_global_to_local(phase_velocity(pos, t, angleParams), t);
    //let W = phase_velocity(pos, t, angleParams);
    let Wspeed = W.mag();

    // target velocity
    let s = Wspeed / (1.0+(Wspeed/this.Vmax)**3.0);

    this.beta -= dt/this.tau_dyna * (this.beta-s);       // update beta
    return NormalizeVector(W).mult(this.beta);
  }

  get_Jacobian_Matrix(pos, t){
    // local basis
    let [gradH, gradH_perp]=gradCurvature(pos, t, angleParams);
    let curv=gradH.mag();
    let sgn=Math.sign(gradH.dot(this.velocity));
    let Local_Jac=[[curv*sgn, 0.0], 
                   [0.0, 0.0]];

    let [e1, e2]=NormalizeVector([gradH, gradH_perp]);
    let U=[[e1.x, e2.x],
           [e1.y, e2.y]];

    return MatrixMatrixProduct(U, MatrixMatrixProduct(Local_Jac, MatrixTranspose(U)));
  }

  InferredParametersGlobalMapping(pos, t){
    // local basis expressed in global basis
    let [e1, e2]=NormalizeVector(gradCurvature(pos, t, angleParams));

    // U: matrice de changement de base de local à global
    // Soit D une matrice 2x2 diagonale exprimée en base locale, alors D en base globale est donné par,
    // D(globale) = U.D(locale).U^T
    let U=[[e1.x, e2.x],
           [e1.y, e2.y]]

    // Matrix rotation
    this.Global_Persistence_Matrix=MatrixMatrixProduct(U, MatrixMatrixProduct(this.Local_Persistence_Matrix, MatrixTranspose(U)));
    this.Global_Diffusion_Matrix=MatrixMatrixProduct(U, MatrixMatrixProduct(this.Local_Diffusion_Matrix, MatrixTranspose(U)));
  }

  From_global_to_local(vecs, t){
    // take as input: a vector or an array of vectors expressed in the global static base (x, y)
    // returns: the same vector(s) expressed in the local basis (e1, e2) 

    const isArray = Array.isArray(vecs);
    const array_of_vectors = isArray ? vecs : [vecs];  // if vecs is a single vector in argument, it is stored in an array
    
    let [e1, e2]=NormalizeVector(gradCurvature(this.position, t, angleParams));
    let U=[[e1.x, e2.x],
           [e1.y, e2.y]];
    
    let rotation_done=[];
    for (let i=0; i<array_of_vectors.length; i++){
      rotation_done.push(MatrixVectorProduct(MatrixTranspose(U), array_of_vectors[i]));
    }
    return isArray ? rotation_done : rotation_done[0];
  }

  From_local_to_global(vecs, t){
    // take as input: a vector or an array of vectors expressed in the local basis (e1, e2)
    // returns: the same vector(s) expressed in the global static base (x, y) 

    const isArray = Array.isArray(vecs);
    const array_of_vectors = isArray ? vecs : [vecs];  // if vecs is a single vector in argument, it is stored in an array
    
    let [e1, e2]=NormalizeVector(gradCurvature(this.position, t, angleParams));
    let U=[[e1.x, e2.x],
           [e1.y, e2.y]];
    
    let rotation_done=[];
    for (let i=0; i<array_of_vectors.length; i++){
      rotation_done.push(MatrixVectorProduct(U, array_of_vectors[i]));
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
