/**
 * File              : curvature.js
 * Original article  : doi
 * Authors           : Gildas Carlin, Ian Manifacier and Jean-Louis Milan
 * Date              : 05/05/2025
 * Description       : Surface and curvature description.
 */
//





// shape function r, offset phi, angle psi and their derivatives
//-----------------------------------------------------------------
function logi(t, param){
  // logistic function A/(1+exp(-w(t-tau))) with param=[A, w, tau]
  return -param[0]/(1.0+Math.exp(-param[1]*(t-param[2])));
}

function angle(t, params=angleParams){
  let psi = -psi0;
  for (const param of params) {
    psi += logi(t, param);
  }
  return psi;  
}

function angleDerivatives(t, params=angleParams){
  let psi_dot=0.0;
  for (const param of params){
    // derivative of logistic function f : df/dt=w*f*(1-f/A)
    let f=-logi(t, param);
    psi_dot += param[1]*f*(1.0-f/param[0]);
  }
  return psi_dot;
}

function phase(t){
  return phi0-kc*t;
}

function phaseDerivatives(t){
  return kc;
}

function shape(pos, t, params){
  let x = pos.x;
  let y = pos.y;
  let psi = angle(t, params);

  return k * (-y*cos(psi)+x*sin(psi));
}

function shapeSpatialDerivatives(pos, t, params){
  let [fx, fxx, fxxx, fy, fyy, fyyy, fxy, fxxy, fxyy] = new Array(9).fill(0.0);
  let psi = angle(t, params);

  fx = k * sin(psi);
  fy = -k * cos(psi);

  return [fx, fxx, fxxx, fy, fyy, fyyy, fxy, fxxy, fxyy];
}

function shapeTimeDerivatives(pos, t, params){
  let [ft, ftx, fty] = new Array(3).fill(0.0);
  let x=pos.x, y=pos.y;

  // offset and its derivative
  let psi=angle(t, params), psi_dot=angleDerivatives(t, params);
  let c=cos(psi), s=sin(psi);

  ft = k*psi_dot*(x*c + y*s);

  ftx = k*psi_dot*c;
  fty = k*psi_dot*s;

  return [ft, ftx, fty];
}



// substrate and its derivatives
//-----------------------------------------------------------------
function substrate(pos, t, params){
  let f = 0.0;
  
  if (surface == 'right_sinus'){  
    // phase and shape
    let phi = phase(t);
    let r = shape(pos, t, params);
    
    f = a*sin(r + phi);
  }
  return f;
}

function substrateSpatialDerivatives(pos, t, params){
  let [fx, fxx, fxxx, fy, fyy, fyyy, fxy, fxxy, fxyy] = new Array(9).fill(0.0);

  if (surface == 'right_sinus'){
    // phase and shape
    let phi = phase(t);
    let r = shape(pos, t, params);
    let arg = r + phi;
    let c = cos(arg), s = sin(arg);
    let [rx, rxx, rxxx, ry, ryy, ryyy, rxy, rxxy, rxyy] = shapeSpatialDerivatives(pos, t, params);

    fx = rx*c;
    fy = ry*c;

    fxx = rxx*c - rx**2.0*s;
    fyy = ryy*c - ry**2.0*s;
    fxy = rxy*c - rx*ry*s;

    fxxx = (rxxx-rx**3.0)*c - 3.0*rxx*rx*s;
    fyyy = (ryyy-ry**3.0)*c - 3.0*ryy*ry*s;
    fxxy = (rxxy-rx**2.0*ry)*c - (rxx*ry+2.0*rxy*rx)*s;
    fxyy = (rxyy-rx*ry**2.0)*c - (ryy*rx+2.0*rxy*ry)*s;
  }

  return [fx, fxx, fxxx, fy, fyy, fyyy, fxy, fxxy, fxyy].map(f => a*f);
}

function substrateTimeDerivatives(pos, t, params){
  let [ft, ftx, fty] = new Array(3).fill(0.0);

  // phase and its derivative
  let phi=phase(t), phi_t=phaseDerivatives(t);

  // shape and its derivatives
  let r = shape(pos, t, params);
  let [rx, rxx, rxxx, ry, ryy, ryyy, rxy, rxxy, rxyy] = shapeSpatialDerivatives(pos, t, params);
  let [rt, rtx, rty] = shapeTimeDerivatives(pos, t, params);
  
  // convenient notation
  let c=cos(r+phi), s=sin(r+phi);

  ft = (rt+phi_t) * cos(r+phi);

  ftx = rtx*c - (rt+phi_t)*rx*s;
  fty = rty*c - (rt+phi_t)*ry*s;

  return [ft, ftx, fty].map(f => a*f);
}

function phase_velocity(pos, t, params){
  let psi = angle(t, params);
  let psi_dot=angleDerivatives(t, params);
  return createVector(sin(psi), -cos(psi)).mult(c).add(createVector(-pos.y, -pos.x).mult(psi_dot)); 
}

//  Curvature metrics
//-----------------------------------------------------------------

function metric_tensor(pos, t, params){
  let [hx, hxx, hxxx, hy, hyy, hyyy, hxy, hxxy, hxyy] = substrateSpatialDerivatives(pos, t, params);

  // metric tensor
  let gxx = 1.0+hx**2.0;
  let gxy = hx*hy;
  let gyy = 1.0+hy**2.0;
  let g = gxx*gyy-gxy**2.0;  // determinant 

  return [g, gxx, gyy, gxy];
}

function curvature_tensor(pos, t, params){
  let [hx, hxx, hxxx, hy, hyy, hyyy, hxy, hxxy, hxyy] = substrateSpatialDerivatives(pos, t, params);
  let [g, gxx, gyy, gxy] = metric_tensor(pos, t, params);
 
  let sqrg = sqrt(g);

  // curvature tensor
  let Kxx = hxx/sqrg;
  let Kxy = hxy/sqrg;
  let Kyy = hyy/sqrg;

  return [Kxx, Kyy, Kxy];
}

function gradCurvature(pos, t, params){
  let [hx, hxx, hxxx, hy, hyy, hyyy, hxy, hxxy, hxyy] = substrateSpatialDerivatives(pos, t, params);
  let [g, gxx, gyy, gxy] = metric_tensor(pos, t, params);

  let u = hxx*gyy + hyy*gxx - 2.0*hxy*gxy;
  let v = 2.0 * g**1.5;

  let ux = gyy*hxxx + gxx*hxyy + 2.0*hxx*(hxy*hy+hyy*hx) - 2.0*(gxy*hxxy+hxy*(hxx*hy+hxy*hx));
  let vx = 6.0 * (hxx*hx+hxy*hy) * g**0.5;

  let uy = gyy*hxxy + gxx*hyyy + 2.0*hyy*(hxx*hy+hxy*hx) - 2.0*(gxy*hxyy+hxy*(hxy*hy+hyy*hx));
  let vy = 6.0 * (hxy*hx+hyy*hy) * g**0.5;

  let gradH = createVector(v*ux - u*vx,  v*uy - u*vy).mult(0.25*g**-3.0);
  let gradH_perp = createVector(-1.0*gradH.y, gradH.x);
  return [gradH, gradH_perp];
}


//  Vector normalization
//-----------------------------------------------------------------

function Normalize_Vec(vec){
  if (vec.mag()<pow(10.0, -10)){
    return createVector(0.0, 0.0);
  } else {
    return vec.normalize();
  }
}










