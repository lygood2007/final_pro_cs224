// IoR of Air is   1.0002926
// IoR of Water is 1.33335  (at 20 deg C)
// n = n2/n1, (n - 1)/(n + 1) = 0.0399

const float Eta = .04; // Ratio of indices of refraction
const float FresnelPower = 5.0;

//const float F  = ((1.0-Eta) * (1.0-Eta)) / ((1.0+Eta) * (1.0+Eta));
const float F  = .0399; //caclulating it out once for water
varying vec3  Reflect;
varying vec3  Refract;
varying float Ratio;

uniform sampler2D NormalMap;

void main()
{

    vec3 vertex  = (gl_ModelViewMatrix * gl_Vertex).xyz;
    vec3 n = normalize(gl_NormalMatrix * gl_Normal);
    vec3 i = normalize(vertex);


    Ratio   = F + (1.0 - F) * pow((1.0 - dot(-i, n)), FresnelPower);

    Refract = refract(i, n, Eta);
    //I should probably have this sample from the terrain and not the skymap
    Refract = vec3(gl_TextureMatrix[0] * vec4(Refract, 1.0));

    Reflect = reflect(i, n);
    Reflect = vec3(gl_TextureMatrix[0] * vec4(Reflect, 1.0));

    gl_Position = ftransform();
}
