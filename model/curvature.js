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

function angle(pos, t, params=angleParams){
  let psi=psi0;

  if (surface_representations[0]=='right_sinus' && angleParams.length!=0){
    for (const param of params){
      psi += LogisticFunction(t, param);
    }
  }

  if (surface_representations[0]=='convr_sinus'){
    let params=[psi0, w, 0.0];
    psi=LogisticFunction(pos.x, params)-0.5*psi;
  }

  return psi;
}

function angleTimeDerivative(t, params=angleParams){
  let psi_dot=0.0;
  for (const param of params){
    // derivative of logistic function f : f'(t)=w*f(t)*(1-f(t)/A)
    let f=LogisticFunction(t, param);
    psi_dot += param[1]*f*(1.0-f/param[0]);
  }
  return psi_dot;
}

function angleSpatialDerivatives(pos, t){
  let [fx, fxx, fxxx, fy, fyy, fyyy, fxy, fxxy, fxyy] = new Array(9).fill(0.0);
  
  if (surface_representations[0]=="convr_sinus"){
    let params=[psi0, w, 0.0];
    let f=LogisticFunction(pos.x, params);
    fx=f*(1.0-f/psi0);
    fxx=fx*(1.0-2.0*f/psi0);
    fxxx=fxx*(1.0-2.0*f/psi0) - 2.0/psi0*fx**2.0;
  }

  return [fx, fxx, fxxx, fy, fyy, fyyy, fxy, fxxy, fxyy].map(f => w*f);
}

function phase(pos, t){
  if (surface_representations[1]==1 || surface_representations[1]==2 || surface_representations[1]==3){
    return phi0-kc*t;
  }

  if (surface_representations[0]=="snake_sinus"){
    return -sin(ks*pos.y + ksc*t);
  }
}



function shape(pos, t, params){
  let x = pos.x, y = pos.y;
  let psi = angle(pos, t, params);

  if (surface_representations[0] == 'right_sinus'){
    return k*(x*cos(psi) + y*sin(psi));
  }

  if (surface_representations[0] == 'wells_sinus'){
    let shape=[x, (x*cos(psi) + y*sin(psi))];
    return shape.map(f => k*f);
  }

  if (surface_representations[0] == 'convr_sinus'){
    return k*(x*sin(psi) + y*cos(psi));
  }

  if (surface_representations[0] == 'snake_sinus'){
    return k*x;
  }
}

function shapeSpatialDerivatives(pos, t, params){
  let [fx, fxx, fxxx, fy, fyy, fyyy, fxy, fxxy, fxyy] = new Array(9).fill(0.0);
  
  if (surface_representations[1] == 1 || surface_representations[1] == 2){
    let psi = angle(pos, t, params);
    fx = cos(psi);
    fy = sin(psi);
  }

  if (surface_representations[0] == 'convr_sinus'){
    let psi = angle(pos, t, params);
    let c=cos(psi), s=sin(psi);
    let x=pos.x, y=pos.y;
    let [px, pxx, pxxx, py, pyy, pyyy, pxy, pxxy, pxyy]=angleSpatialDerivatives(pos, t);

    fx=s + x*px*c - y*px*s;
    fy=c + x*py*c - y*py*s;

    fxx=2.0*px*c + x*(pxx*c-px**2.0*s) - y*(pxx*s+px**2.0*c);
    fxy=py*c-px*s + x*(pxy*c-px*py*s) - y*(pxy*s+px*py*c);
    fyy=-2.0*py*s - y*(pyy*s+py**2.0*c) + x*(pyy*c-py**2.0*s);

    fxxx=2.0*(pxx*c-px**2.0*s) + pxx*c-px**2.0*s + x*(pxxx*c-3.0*pxx*px*s-px**3.0*c) - y*(pxxx*s+3.0*pxx*px*c-px**3.0*s);
    fxxy=2.0*(pxy*c-px*py*s) - (pxx*s+px**2.0*c) + x*(pxxy*c-pxx*py*s-2.0*pxy*px*s-px**2.0*py*c) - y*(pxxy*s+pxx*py*c-2.0*pxy*px*c-px**2.0*py*s);
    fxyy=-2.0*(pxy*s+px*py*c) + pyy*c-py**2.0*s - y*(pxyy*s+pyy*px*c+2.0*pxy*py*c-py**2.0*px*s) + x*(pxyy*c-pyy*px*s-2.0*pxy*py*s-py**2.0*px*c);
    fyyy=-2.0*(pyy*s+py**2.0*c) - (pyy*s+py**2.0*c) - y*(pyyy*s+3.0*pyy*py*c-py**3.0*s) + x*(pyyy*c-3.0*pyy*py*s-py**3.0*c);
  }

  if (surface_representations[0] == 'snake_sinus'){

  }

  return [fx, fxx, fxxx, fy, fyy, fyyy, fxy, fxxy, fxyy].map(f => k*f);
}

//function shapeTimeDerivatives(pos, t, params){
//  let [ft, ftx, fty] = new Array(3).fill(0.0);
//  let x=pos.x, y=pos.y;
//
//  // offset and its derivative
//  let psi=angle(t, params), psi_dot=angleDerivatives(t, params);
//  let c=cos(psi), s=sin(psi);
//
//  ft = k*psi_dot*(x*c + y*s);
//
//  ftx = k*psi_dot*c;
//  fty = k*psi_dot*s;
//
//  return [ft, ftx, fty];
//}



// substrate and its derivatives
//-----------------------------------------------------------------
function substrate(pos, t, params){
  let f = 0.0;
  
  if (surface_representations[1]==1 || surface_representations[1]==3){
    let phi = phase(pos, t);
    let r = shape(pos, t, params);
    
    f = sin(r + phi);
  }

  if (surface_representations[0] == 'wells_sinus'){
    let phi = phase(pos, t);
    let r = shape(pos, t, params);

    f=0.5*(sin(r[0]) + sin(r[1] + phi));
  }

  if (surface_representations[0] == 'snake_sinus'){
    let phi = phase(pos, t);
    let r = shape(pos, t, params);
    f=sin(r+phi);
  }

  return a*f;
}

function substrateSpatialDerivatives(pos, t, params){
  let [fx, fxx, fxxx, fy, fyy, fyyy, fxy, fxxy, fxyy] = new Array(9).fill(0.0);

  if (surface_representations[1]==1 || surface_representations[1]==3){
    let phi = phase(pos, t);
    let r = shape(pos, t, params);
    let arg = r+phi;
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

  if (surface_representations[0] == 'wells_sinus'){
    let phi = phase(pos, t);
    let r = shape(pos, t, params);
    let cx=cos(r[0]), sx=sin(r[0]);
    let cy=cos(r[1]+phi), sy=sin(r[1]+phi);
    let [rx, rxx, rxxx, ry, ryy, ryyy, rxy, rxxy, rxyy] = shapeSpatialDerivatives(pos, t, params);

    fx = 0.5*(k*cx + rx*cy);
    fy = 0.5*ry*cy;

    fxx = -0.5*(k**2.0*sx + rx**2.0*sy);
    fyy = -0.5*ry**2.0*sy;
    fxy = -0.5*rx*ry*sy;

    fxxx = -0.5*(k**3.0*cx + rx**3.0*cy);
    fyyy = -0.5*ry**3.0*cy;
    fxxy = -0.5*rx**2.0*ry*cy;
    fxyy = -0.5*rx*ry**2.0*cy;
  }

  return [fx, fxx, fxxx, fy, fyy, fyyy, fxy, fxxy, fxyy].map(f => a*f);
}

//function substrateTimeDerivatives(pos, t, params){
//  let [ft, ftx, fty] = new Array(3).fill(0.0);
//
//  if (surface_representations[0] == 'right_sinus'){
//    // phase and its derivative
//    let phi=phase(t), phi_t=phaseDerivatives(t);
//
//    // shape and its derivatives
//    let r = shape(pos, t, params);
//    let [rx, rxx, rxxx, ry, ryy, ryyy, rxy, rxxy, rxyy] = shapeSpatialDerivatives(pos, t, params);
//    let [rt, rtx, rty] = shapeTimeDerivatives(pos, t, params);
//    let c=cos(r+phi), s=sin(r+phi);
//
//    ft = (rt+phi_t) * cos(r+phi);
//
//    ftx = rtx*c - (rt+phi_t)*rx*s;
//    fty = rty*c - (rt+phi_t)*ry*s;
//  }
//
//  if (surface_representations[0] == 'wells_sinus'){
//    let phi=phase(t), phi_t=phaseDerivatives(t);
//    let r = shape(pos, t, params);
//    let cy=cos(r[1]+phi), sy=sin(r[1]+phi);
//
//    ft = 0.5*phi_t*cy;
//    fty = -0.5*k*phi_t*sy;
//  }
//  
//  return [ft, ftx, fty].map(f => a*f);
//}

function phase_velocity(pos, t, params){
  let psi = angle(pos, t, params);
  let psi_dot=angleTimeDerivative(t, params);
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

function meanCurvature(pos, t, params){
  let [hx, hxx, hxxx, hy, hyy, hyyy, hxy, hxxy, hxyy] = substrateSpatialDerivatives(pos, t, params);
  return 0.5 * ((1.0+hx**2.0)*hyy + (1.0+hy**2.0)*hxx - 2.0*hx*hy*hxy) / (1.0+hx**2.0+hy**2.0)**1.5;
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














