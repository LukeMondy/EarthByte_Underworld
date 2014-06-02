//#version 120 //Required to use gl_PointCoord
uniform float pointScale;   // scale to calculate size in pixels
uniform int pointDist;   // Scale by distance

attribute float aSize;
attribute float aPointType;

varying float vSmooth;
varying float vPointType;
varying vec3 vPosEye;

void main(void)
{
   vSmooth = aSize;
   float pSize = abs(aSize);

   // calculate window-space point size
   vec3 posEye = vec3(gl_ModelViewMatrix * gl_Vertex);
   float dist = 1.0;
   if (pointDist > 0)
      dist = length(posEye);
   //Limit scaling, overly large points are very slow to render
   gl_PointSize = max(1.0, min(40.0, pointScale * pSize / dist));
   gl_Position = gl_ModelViewProjectionMatrix * gl_Vertex;
   gl_FrontColor = gl_Color;
   vPosEye = posEye;
   vPointType = aPointType;
}

