#version 120 //Required to use gl_PointCoord
//uniform sampler2D texture;
varying float vSmooth;
varying vec3 vPosEye;
varying float vPointType;
uniform int uPointType;
uniform float uOpacity;
void main(void)
{
   float alpha = gl_Color.a;
   if (uOpacity > 0.0) alpha *= uOpacity;
   gl_FragColor = gl_Color;
   float pointType = uPointType;
   if (vPointType >= 0) pointType = vPointType;

   //Flat, square points, fastest
   if (pointType == 4 || vSmooth < 0.0) 
      return;

   // calculate normal from point/tex coordinates
   vec3 N;
   N.xy = gl_PointCoord * vec2(2.0, -2.0) + vec2(-1.0, 1.0);
   float mag = dot(N.xy, N.xy);
   if (alpha < 0.01 || mag > 1.0) discard;   // kill pixels outside circle radius and transparent pixels

   if (pointType < 2)
   {
      if (pointType == 0)
         gl_FragColor.a = alpha * 1.0-sqrt(mag);  //Gaussian
      else
         gl_FragColor.a = alpha * 1.0-mag;      //Linear
      return;
   }
   N.z = sqrt(1.0-mag);

   // calculate lighting
   vec3 lightDir = normalize(vec3(1,1,1) - vPosEye);
   float diffuse = max(0.0, dot(lightDir, N));

   // compute the specular term if diffuse is larger than zero 
   vec3 specular = vec3(0.0,0.0,0.0);
   if (pointType == 3 && diffuse > 0.0)
   {
      vec3 lightPos = lightDir*2.0; //vec3(1.0, 1.0, 1.0);        //Fixed light position
      float shininess = 32.0;
      vec3 specolour = vec3(0.5, 0.5, 0.5);   //Color of light
      // normalize the half-vector, and then compute the 
      // cosine (dot product) with the normal
      //vec3 halfVector = normalize(lightPos - vec3(gl_PointCoord, 0));
      vec3 halfVector = normalize(lightPos - gl_TexCoord[0].xyz);
      float NdotHV = max(dot(N, halfVector), 0.0);
      specular = specolour * pow(NdotHV, shininess);
   }

   gl_FragColor = vec4(gl_Color.rgb * diffuse + specular, alpha);
   //gl_FragColor = vec4(gl_Color.rgb * diffuse, alpha);
}
