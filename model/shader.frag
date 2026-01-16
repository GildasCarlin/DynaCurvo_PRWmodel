/**
 * File              : shader.frag
 * Original article  : doi
 * Authors           : Gildas Carlin, Ian Manifacier and Jean-Louis Milan
 * Date              : 05/05/2025
 * Description       : Background sketch (visual only). 
 */
//

#ifdef GL_ES
precision mediump float;
#endif

#define MAX_ANGLES 10 // Taille maximale

uniform vec2 u_resolution;
uniform int surface;
uniform float k;
uniform float kc;
uniform float ks;
uniform float ksc;
uniform float phi0;
uniform float psi0;
uniform float w;
uniform float t;
uniform float s;
uniform int num_angles;
uniform vec3 temporal_angle[MAX_ANGLES];


float angle(float x, float y){
  float psi = psi0;
  if (surface==1){
    for (int i = 0; i < MAX_ANGLES; i++){
      if (i >= num_angles) break; 

      float A = temporal_angle[i].x; 
      float speed = temporal_angle[i].y; 
      float tau = temporal_angle[i].z; 
      psi += A/(1.0+exp(-speed*(t-tau)));
      }
  }
  
  if (surface==3){
  psi=psi/(1.0+exp(-w*x))-0.5*psi;  //
  }

  return psi;
}

float phase(float x, float y){

  if (surface==1 || surface==2 || surface==3){
    return phi0 - kc*t;
  }

  if (surface==4){  // snake_sinus
    return sin(ks*y + ksc*t);
  }
}

float Heaviside(float y, float m){
  if(y>0.0){
    return 0.0;
  }else if(y<0.0){
    return 1.0;
  }else{
    return m;
  }
}

float surface_equation(float x, float y){
  // x and y are in the shader frame, they are transformed in the js frame 
  x -= 0.5*u_resolution.x;
  y = -1.0*y+0.5*u_resolution.y;

  // In shaders, the RGB color goes from 0 - 1 instead of 0 - 255. Make sure to normalize the height.
  float m = 0.0;     // flat

  if (surface == 1){ // right sinus
    float phi = phase(x, y);
    float psi = angle(x, y);
    m = 0.5 * (1.0 + sin(k*(y*sin(psi) + x*cos(psi)) + phi));
  }

  if (surface==2){  // wells_sinus
    float phi = phase(x, y);
    float psi = angle(x, y);
    m = 0.5 + 0.25*(sin(k*x) + sin(k*(y*sin(psi) + x*cos(psi)) + phi));
  }

  if (surface==3){  // convr_sinus
    float phi = phase(x, y);
    float psi = angle(x, y);
    m = 0.5 * (1.0 + sin(k*(x*sin(psi) + y*cos(psi)) - phi));
  }

  if (surface==4){   // snake_sinus
    float phi = phase(x, y);
    m = 0.5 * (1.0 + sin(k*x + phi));
  }
  return m;
}

void main() {
  vec2 st = gl_FragCoord.xy/u_resolution.xy; 
  // Use the x position of the pixel to define its hue (with saturation and brightness to 1)
  // Then convert the HSB color to RGB
  
  // Color the pixel
  float alpha = 1.0;
  float h = surface_equation(gl_FragCoord.x, gl_FragCoord.y);

  vec3 color_max = vec3(0.976, 0.976, 0.976);  // convex
  vec3 color_min = vec3(0.063, 0.396, 0.671);  // concave
  vec3 color = mix(color_min, color_max, h);

  if (h>0.9995){
    color = vec3(0.95, 0.42, 0.6);  
  }

  if (h<0.0005){
    color = vec3(0.79, 0.95, 0.07);  
  }

  
  //float x0=0.0;
  //float y0=0.0;
  //x0 -= 0.5*u_resolution.x;
  //y0 += -1.0*y0+0.5*u_resolution.y;

  //float Dx_p=10.0;
  //float Dx_m=10.0;
  //float Dy_p=10.0;
  //float Dy_m=10.0;
  //if (gl_FragCoord.x<x0+Dx_p && gl_FragCoord.x>x0-Dx_m){
  //  if (gl_FragCoord.y<y0+Dy_p && gl_FragCoord.y>y0-Dy_m){
  //    color = vec3(0.95, 0.42, 0.6);
  //  }
  //}

  gl_FragColor = vec4(color, alpha);

}
