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
uniform float phi0;
uniform float psi0;
uniform float t;
uniform float s;
uniform int num_angles;
uniform vec3 temporal_angle[MAX_ANGLES];


float angle(){
  float psi = psi0;
  for (int i = 0; i < MAX_ANGLES; i++) {
    if (i >= num_angles) break; 

    float A = temporal_angle[i].x; 
    float speed = temporal_angle[i].y; 
    float tau = temporal_angle[i].z; 
    psi += A/(1.0+exp(-speed*(t-tau)));
  }
  return psi;
}

float phase(){
  return 3.141595 + phi0 + kc*t;
}

float surface_equation(float x, float y){
  // x and y are in the shader frame, they are transformed in the js frame 
  x -= 0.5*u_resolution.x;
  y = -1.0*y + 0.5*u_resolution.y;

  // In shaders, the RGB color goes from 0 - 1 instead of 0 - 255. Make sure to adapt the surface equation consequently.
  float m = 0.0;
  if (surface == 1){
    float phi = phase();
    float psi = angle();
    m = 0.5 * (1.0 + sin(k*(y*cos(psi) + x*sin(psi)) + phi));
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

  //if (h<0.002){
  //  color = vec3(0.95, 0.42, 0.6);  
  //}


  
  gl_FragColor = vec4(color, alpha);

}
