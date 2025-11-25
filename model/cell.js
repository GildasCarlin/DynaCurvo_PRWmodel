
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

    this.Px = set_of_parameters[0];
    this.Py = set_of_parameters[1];
    this.Sx = set_of_parameters[2];
    this.Sy = set_of_parameters[3];
    this.alpha = set_of_parameters[4];
    this.q = set_of_parameters[10];
    this.tau_stat = set_of_parameters[7];
    this.tau_dyna = set_of_parameters[8];

    this.Vmed = set_of_parameters[5];
    this.Vmean = set_of_parameters[6];
    this.Vmax = set_of_parameters[9];
    this.central_V_flat = central_vel;
    
    this.inv_polarityMatrix = [
      [1.0/this.Px, 0.0],
      [0.0, 1.0/this.Py], 
    ];

    this.diffusionMatrix = [ 
      [this.Sx / Math.sqrt(this.Px), 0.0],
      [0.0, this.Sy / Math.sqrt(this.Py)]
    ];

    this.beta = 0.0;
    this.f = 0.0;

    this.f_ext = createVector(0.0, 0.0);
    this.velocity = createVector(0.0, 0.0);
    this.position = createVector(400.0 * Math.random() - 200.0, 400.0 * Math.random() - 200.0);

    // cell shape (visual only)
    this.cell_radius_min = 10.0;  
    this.cell_radius_max = 15.0;

    let [gradH, gradH_perp] = gradCurvature(this.position, 0.0, angleParams);    // local basis
    this.theta=createVector(1.0, 0.0).angleBetween(gradH);
    
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

    // external forces caused by curvature ... 
    if (curvature_assumption){
      this.f_ext = this.get_Psi_curvature(this.position.copy(), this.velocity.copy(), t);
    } else {
      this.f_ext = createVector(0.0, 0.0);
    }

    //let test=this.f_ext.mag();
    //console.log(test);
    //drawVector(this.position.x, this.position.y, this.f_ext.copy().mult(1000));

    
    if (surf_assumption){
      // ... and dynamic wave
      this.f_ext.add(this.get_Psi_wave(this.position.copy(), t));
    }

    // P and S parameters are exressed in the local basis, they are expressed on the global x, y basis
    let [gradH, gradH_perp] = gradCurvature(this.position, t, angleParams);    // local basis
    let [Kxx, Kyy, Kxy]=curvature_tensor(this.position, t);
    //console.log(Kxx.toFixed(2), Kyy.toFixed(2), Kxy.toFixed(2))
    //console.log(Math.sqrt(Kxx**2.0+ Kyy**2.0+ 2.0*Kxy**2.0).toFixed(2), Math.sqrt(gradH.x**2.0+ gradH.y**2.0).toFixed(2));
    //console.log(Kyy, Kyy.toFixed(2));
    

    this.theta=Math.abs(createVector(1.0, 0.0).angleBetween(gradH));
    if (this.theta > Math.PI/2) {
      this.theta = Math.PI - this.theta;
    }

    let c2=Math.cos(this.theta)**2.0; let s2=Math.sin(this.theta)**2.0;

    let a=1.0/this.Px; let b=1.0/this.Py;
    this.inv_polarityMatrix = [
      [a*c2+b*s2, 0.0],
      [0.0, b*c2+a*s2], 
    ];

    a=this.Sx / Math.sqrt(this.Px); b=this.Sy / Math.sqrt(this.Py);

    this.diffusionMatrix = [ 
      [a*c2+b*s2, 0.0],
      [0.0, b*c2+a*s2]
    ];

    // velocity computation
    this.velocity = this.velocity.add(
      this.MultiplyMatrixWVector(this.get_Jacobian_Matrix(this.position.copy(), t), this.velocity.copy()).mult(dt*this.q)
    ).sub(
      this.MultiplyMatrixWVector(this.inv_polarityMatrix, this.velocity.copy().sub(this.f_ext.copy())).mult(dt)
    ).add(
      this.MultiplyMatrixWVector(this.diffusionMatrix, this.get_random_vector()).mult(Math.sqrt(dt))
    );


    // velocity integration to get cell position
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

  get_Psi_curvature(pos, pol, t){
    // init
    let n = Normalize_Vec(pol);  // unit-length direction of migration
    let u = n.x, v = n.y;        // coordinates

    // curvature
    let [Kxx, Kyy, Kxy] = curvature_tensor(pos, t, angleParams);     //
    let kappa = Kxx*u**2.0 + Kyy*v**2.0 + 2.0*Kxy*u*v;  // directional curvature
    let [gradH, gradH_perp] = gradCurvature(pos, t, angleParams);    // local basis
    
    // direction of Psi if the cell is in concave
    let Dtheta = 0.0;
    if (kappa>0.0){
      // is the cell polarized along +gradH_perp (True) or -gradH_perp (False) ?
      let test = pol.dot(gradH_perp)<0.0;
      if (test){
        gradH_perp.mult(-1.0);
      }
      Dtheta = pol.angleBetween(gradH_perp) * this.alpha;
    }

    // norm of psi
    let f_d = Math.abs(this.Vmean - this.central_V_flat);
    this.f -= dt/this.tau_stat * (this.f-f_d);          // update f
    
    let w = this.multiplyRotationMatrixVector(Dtheta, n);
    return w.mult(this.f);
  }

  get_Psi_wave(pos, t){
    let W = phase_velocity(pos, t, angleParams);
    let Wspeed = W.mag();

    // target velocity
    let s = Wspeed / (1.0+(Wspeed/this.Vmax)**3.0);

    this.beta -= dt/this.tau_dyna * (this.beta-s);       // update beta
    return Normalize_Vec(W).mult(this.beta);
  }

  get_Jacobian_Matrix(pos, t){
    //let [Kxx, Kyy, Kxy]=curvature_tensor(pos, t);
    //let jacobian = [
    //  [Kxx, Kxy],
    //  [Kxy, Kyy]
    //];

    let [gradH, gradH_perp] = gradCurvature(pos, t, angleParams);
    let norm=gradH.mag();
    let sgn=Math.sign(gradH.dot(this.velocity));
    let theta=Math.abs(gradH.angleBetween(this.velocity));
    if (this.theta > Math.PI/2) {
      this.theta = Math.PI - this.theta;
    }
    let jacobian = [
      [norm*Math.sin(theta)*sgn, 0.0],
      [0.0, norm*Math.cos(theta)*sgn]
    ];
    //console.log(Math.sin(theta), Math.cos(theta));
    return jacobian;
  }

  MultiplyMatrixWVector(A, v){
    return createVector(A[0][0]*v.x + A[0][1]*v.y, A[1][0]*v.x + A[1][1]*v.y);
  }

  Heaviside_step_function(s){
    if (s<=0.0){
      return 0.0;
    } else {
      return 1.0;
    }
  }

  multiplyRotationMatrixVector(angle, v) {
    let c = Math.cos(angle), s = Math.sin(angle);
    let x = v.x, y = v.y;
    return createVector(c*x - s*y, s*x + c*y);
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

    let a = this.cell_radius_min + this.cell_radius_max*(1.0-exp(-vel.mag()/(0.5*(this.Sx+this.Sy))));
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
    table_data.push(`${it_sauv}, ${t}, ${this.cellID}, ${this.position.x}, ${this.position.y}, ${this.theta}`);
  }
  
}
