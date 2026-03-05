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

  if (surface_representations[0]=='right_sinus' && angleParams.length>0){
    for (const param of params){
      psi += LogisticFunction(t, param);
    }
  }

  if (surface_representations[1]==2 || surface_representations[1]==3 || surface_representations[1]==4){
    psi=0.0-HALF_PI;
  }

  return psi;
}

function angleTimeDerivative(t, params=angleParams){
  let psi_dot=0.0;
  if (surface_representations[0]=='right_sinus' && angleParams.length>0){
    for (const param of params){
      // derivative of logistic function f : f'(t)=w*f(t)*(1-f(t)/A)
      let f=LogisticFunction(t, param);
      psi_dot += param[1]*f*(1.0-f/param[0]);
    }
  }

  return psi_dot;
}

function angleSpatialDerivatives(pos, t){
  let [fx, fxx, fxxx, fy, fyy, fyyy, fxy, fxxy, fxyy] = new Array(9).fill(0.0);
  
  //if (surface_representations[0]=="convr_sinus"){
  //  let params=[1.0, w, 0.0];
  //  let f=LogisticFunction(pos.x, params);
  //  fx=psi0*f*(1.0-f);
  //  fxx=psi0*fx*(1.0-2.0*f);
  //  fxxx=psi0*(fxx-2.0*(fx**2.0+f*fxx));
  //}

  return [fx, fxx, fxxx, fy, fyy, fyyy, fxy, fxxy, fxyy].map(f => w*f);
}

function phase(pos, t){
  if (surface_representations[1]==1 || surface_representations[1]==2){
    return phi0-kc*t;
  }

  if (surface_representations[0]=="convr_sinus"){
    let X=abs(pos.x);
    let U=HeavisideStepFunction(X-psi0);
    return phi0+kc*t + (1.0-U)*w*X**2.0 + U*w*psi0*(2.0*X-psi0);
  }

  if (surface_representations[0]=="snake_sinus"){
    return as*sin(ks*pos.y + ksc*t);
  }
}

function phaseSpatialDerivatives(pos, t){
  let [fx, fxx, fxxx, fy, fyy, fyyy, fxy, fxxy, fxyy] = new Array(9).fill(0.0);

  if (surface_representations[0]=="convr_sinus"){
    let X=abs(pos.x);
    let U=HeavisideStepFunction(X-psi0);

    fx=Math.sign(pos.x)*2.0*w*((1.0-U)*X + U*psi0);
    fxx=2.0*w*(1.0-U);
  }
  
  if (surface_representations[0]=="snake_sinus"){
    let arg=ks*pos.y+ksc*t;
    fy=as*ks*cos(arg);
    fyy=-1.0*as*ks**2.0*sin(arg);
    fyyy=-1.0*as*ks**3.0*cos(arg);
  }

  return [fx, fxx, fxxx, fy, fyy, fyyy, fxy, fxxy, fxyy];
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
    return k*y;
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
    fy=1.0;
    //let psi = angle(pos, t, params);
    //let c=cos(psi), s=sin(psi);
    //let x=pos.x, y=pos.y;
    //let [px, pxx, pxxx, py, pyy, pyyy, pxy, pxxy, pxyy]=angleSpatialDerivatives(pos, t);
//
    //fx=s + px*(x*c - y*s);
    //fy=c;//c + x*py*c - y*py*s;
//
    //fxx=2.0*px*c + x*(pxx*c-px**2.0*s) - y*(pxx*s+px**2.0*c);
    //fxy=-px*s;//py*c-px*s + x*(pxy*c-px*py*s) - y*(pxy*s+px*py*c);
    //fyy=0.0; //-2.0*py*s - y*(pyy*s+py**2.0*c) + x*(pyy*c-py**2.0*s);
//
    //fxxx=2.0*(pxx*c-px**2.0*s) + pxx*c-px**2.0*s + x*(pxxx*c-3.0*pxx*px*s-px**3.0*c) - y*(pxxx*s+3.0*pxx*px*c-px**3.0*s);
    //fxxy=-(pxx*s+px**2.0*c);  //2.0*(pxy*c-px*py*s) - (pxx*s+px**2.0*c) + x*(pxxy*c-pxx*py*s-2.0*pxy*px*s-px**2.0*py*c) - y*(pxxy*s+pxx*py*c-2.0*pxy*px*c-px**2.0*py*s);
    //fxyy=0.0;//-2.0*(pxy*s+px*py*c) + pyy*c-py**2.0*s - y*(pxyy*s+pyy*px*c+2.0*pxy*py*c-py**2.0*px*s) + x*(pxyy*c-pyy*px*s-2.0*pxy*py*s-py**2.0*px*c);
    //fyyy=0.0;//-2.0*(pyy*s+py**2.0*c) - (pyy*s+py**2.0*c) - y*(pyyy*s+3.0*pyy*py*c-py**3.0*s) + x*(pyyy*c-3.0*pyy*py*s-py**3.0*c);
  }

  if (surface_representations[0] == 'snake_sinus'){
    fx=1.0;
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
  
  if (surface_representations[1]==1 || surface_representations[1]==3 || surface_representations[1]==4){
    let phi = phase(pos, t);
    let r = shape(pos, t, params);
    f = sin(r + phi);
  }

  if (surface_representations[0] == 'wells_sinus'){
    let phi = phase(pos, t);
    let r = shape(pos, t, params);
    f=0.5*(sin(r[0]) + sin(r[1] + phi));
  }

  return a*f;
}

function substrateSpatialDerivatives(pos, t, params){
  let [fx, fxx, fxxx, fy, fyy, fyyy, fxy, fxxy, fxyy] = new Array(9).fill(0.0);

  if (surface_representations[1]==1 || surface_representations[1]==3 || surface_representations[1]==4){
    let phi=phase(pos, t), r=shape(pos, t, params);

    let c=cos(r+phi), s=sin(r+phi);
    let [rx, rxx, rxxx, ry, ryy, ryyy, rxy, rxxy, rxyy]=shapeSpatialDerivatives(pos, t, params);
    let [px, pxx, pxxx, py, pyy, pyyy, pxy, pxxy, pxyy]=phaseSpatialDerivatives(pos, t);
    
    fx = (rx+px)*c;
    fy = (ry+py)*c;

    fxx = (rxx+pxx)*c - (rx+px)**2.0*s;
    fyy = (ryy+pyy)*c - (ry+py)**2.0*s;
    fxy = (rxy+pxy)*c - (rx+px)*(ry+py)*s;

    fxxx = (rxxx+pxxx-(rx+px)**3.0)*c - 3.0*(rxx+pxx)*(rx+px)*s;
    fyyy = (ryyy+pyyy-(ry+py)**3.0)*c - 3.0*(ryy+pyy)*(ry+py)*s;
    fxxy = (rxxy+pxxy-(rx+px)**2.0*(ry+py))*c - ((rxx+pxx)*(ry+py)+2.0*(rxy+pxy)*(rx+px))*s;
    fxyy = (rxyy+pxyy-(rx+px)*(ry+py)**2.0)*c - ((ryy+pyy)*(rx+px)+2.0*(rxy+pxy)*(ry+py))*s;
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

function phaseVelocity(pos, t, params){
  let psi = angle(pos, t, params);
  let psi_dot=angleTimeDerivative(t, params);
  return createVector(cos(psi), sin(psi)).mult(c).add(createVector(-pos.y, pos.x).mult(psi_dot));
}

//  Curvature metrics
//-----------------------------------------------------------------

function metricTensor(pos, t, params){
  let [hx, hxx, hxxx, hy, hyy, hyyy, hxy, hxxy, hxyy] = substrateSpatialDerivatives(pos, t, params);

  // metric tensor
  let gxx = 1.0+hx**2.0;
  let gxy = hx*hy;
  let gyy = 1.0+hy**2.0;
  let g = gxx*gyy-gxy**2.0;  // determinant 

  return [g, gxx, gyy, gxy];
}

function curvatureTensor(pos, t, params){
  let [hx, hxx, hxxx, hy, hyy, hyyy, hxy, hxxy, hxyy] = substrateSpatialDerivatives(pos, t, params);
  let [g, gxx, gyy, gxy] = metricTensor(pos, t, params);
 
  return [hxx, hyy, hxy].map(K => K/sqrt(g));
}

function Christoffel(pos, t, params){
  let [hx, hxx, hxxx, hy, hyy, hyyy, hxy, hxxy, hxyy] = substrateSpatialDerivatives(pos, t, params);
  let g=1.0+hx**2.0+hy**2.0;

  return [hxx*hx, hxy*hx, hyy*hx, hxx*hy, hxy*hy, hyy*hy].map(G => G/g);
}

function meanCurvature(pos, t, params){
  let [hx, hxx, hxxx, hy, hyy, hyyy, hxy, hxxy, hxyy] = substrateSpatialDerivatives(pos, t, params);
  return 0.5 * ((1.0+hx**2.0)*hyy + (1.0+hy**2.0)*hxx - 2.0*hx*hy*hxy) / (1.0+hx**2.0+hy**2.0)**1.5;
}

function gaussianCurvature(pos, t, params){
  let [hx, hxx, hxxx, hy, hyy, hyyy, hxy, hxxy, hxyy] = substrateSpatialDerivatives(pos, t, params);
  return (hxx*hyy-hxy**2.0) / (1.0+hx**2.0+hy**2.0)**2.0;

}

function gradGaussianCurvature(pos, t, params){
  let [hx, hxx, hxxx, hy, hyy, hyyy, hxy, hxxy, hxyy] = substrateSpatialDerivatives(pos, t, params);

  let Kx=(hxxx*hyy+hxyy*hxx-2.0*hxxy*hxy)*(1.0+hx**2.0+hy**2.0) - 4.0*(hxx*hyy-hxy**2.0)*(hxx*hx+hxy*hy);
  let Ky=(hyyy*hxx+hxxy*hyy-2.0*hxyy*hxy)*(1.0+hx**2.0+hy**2.0) - 4.0*(hxx*hyy-hxy**2.0)*(hyy*hy+hxy*hx);

  return createVector(Kx, Ky).div((1.0+hx**2.0+hy**2.0)**3.0);
}

function gradCurvature(pos, t, params){
  // returns the gradient of mean curvature in the euclidean (x,y) plane
  // returns also its perpendicular vector which forms a local basis
  let [hx, hxx, hxxx, hy, hyy, hyyy, hxy, hxxy, hxyy] = substrateSpatialDerivatives(pos, t, params);
  let [g, gxx, gyy, gxy] = metricTensor(pos, t, params);

  let u = hxx*gyy + hyy*gxx - 2.0*hxy*gxy;
  let v = 2.0*g**1.5;

  let ux = gyy*hxxx + gxx*hxyy + 2.0*hxx*(hxy*hy+hyy*hx) - 2.0*(gxy*hxxy+hxy*(hxx*hy+hxy*hx));
  let vx = 6.0 * (hxx*hx+hxy*hy) * g**0.5;

  let uy = gyy*hxxy + gxx*hyyy + 2.0*hyy*(hxx*hy+hxy*hx) - 2.0*(gxy*hxyy+hxy*(hxy*hy+hyy*hx));
  let vy = 6.0 * (hxy*hx+hyy*hy) * g**0.5;

  let gradH = createVector(v*ux - u*vx,  v*uy - u*vy).mult(0.25*g**-3.0);
  let gradH_perp = createVector(-1.0*gradH.y, gradH.x);
  return [gradH, gradH_perp];
}

function TangentLocalBasis(pos, t, params){
  let [hx, hxx, hxxx, hy, hyy, hyyy, hxy, hxxy, hxyy] = substrateSpatialDerivatives(pos, t, params);
  return [createVector(1.0, 0.0, hx), createVector(0.0, 1.0, hy), createVector(-hx, -hy, 1.0)];

}

function MetricDotProduct(pos, t, params, V1, V2){
  // returns the dot product of two vectors expressed in the local tangent plane.
  let [g, gxx, gyy, gxy] = metricTensor(pos, t, params);
  return gxx*V1.x*V2.x + gxy*V1.x*V2.y + gxy*V1.y*V2.x + gyy*V1.y*V2.y;
}

function TangentGradCurvature(pos, t, params){
  // returns the gradient of mean curvature in the local tangent plane 
  // returns also its perpendicular vector (included in the plane) which forms a local basis. Both are normalized.
  let [gradH, gradH_perp] = gradCurvature(pos, t, params);
  let [g, gxx, gyy, gxy] = metricTensor(pos, t, params);

  let g_contra=[[gyy, -gxy], [-gxy, gxx]].map(row => row.map(val => val/g));
  let gradSH = MatrixVectorProduct(g_contra, gradH);
  let norm = sqrt(MetricDotProduct(pos, t, params, gradSH, gradSH));
  gradSH=gradSH.div(norm);

  let gradSH_perp = createVector(-1.0*gradSH.y, gradSH.x);
  gradSH_perp=gradSH_perp.div(sqrt(MetricDotProduct(pos, t, params, gradSH_perp, gradSH_perp)));

  return [norm, gradSH, gradSH_perp];
}

function PrincipalCurvature(pos, t, params){
  [kxx, kyy, kxy] = curvatureTensor(pos, t, params);

  let sqr=sqrt(kxx**2.0+kyy**2.0-2.0*kxx*kyy+4.0*kxy**2.0);
  let k1=0.5*(kxx+kyy + sqr);
  let k2=0.5*(kxx+kyy - sqr);
  let kmax=Math.max(k1, k2);
  let kmin=Math.min(k1, k2);

  let v1=NormalizeVector(createVector(1.0, -(kxx-k1)/kxy));
  let v2=createVector(1.0, -(kxx-k2)/kxy).normalize();

  //sqr=sqrt(kxy**2.0 - kyy*(kxx-kmax));
  //let tan1 = (-kxy + sqr)/kyy;
  //let tan2 = (-kxy - sqr)/kyy;
//
  //let theta1=atan(tan1);
  //let theta2=atan(tan2);
  //let v1=createVector(cos(theta1), sin(theta1));
  //let v2=createVector(cos(theta2), sin(theta2));

  return[k1, k2, v1, v2];
}















