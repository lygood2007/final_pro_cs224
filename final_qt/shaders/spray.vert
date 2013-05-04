
//varying vec3 n, lightDir, r;
//const vec3 L = vec3(0.,0.,1.);

//void main()
//{
//        vec3 vertex = vec3(gl_ModelViewMatrix * gl_Vertex).xyz;
//        n = normalize( gl_NormalMatrix * gl_Normal );
//        vec3 i = normalize(vertex);
//        lightDir = vec3(L - vVertex);

//         r = reflect(I,normal);

//    gl_Position = ftransform();
////    gl_PointSize = 64.0/gl_Position.w;
////    gl_PointSize  = 64.0;
//}
// IoR of Air is   1.0002926
// IoR of Water is 1.33335  (at 20 deg C)
// (n1-n2)/(n1+n2) = -0.1427199692

const float Eta = 0.143;         // Ratio of indices of refraction
const float EtaR = 0.75;
const float EtaG = 0.85;
const float EtaB = 0.9;

const float FresnelPower = 5.0;

const float F  = ((1.0-Eta) * (1.0-Eta)) / ((1.0+Eta) * (1.0+Eta));

varying vec3  Reflect;
varying vec3  Refract;
varying vec3  RefractR;
varying vec3  RefractG;
varying vec3  RefractB;
varying float Ratio;

uniform sampler2D NormalMap;

void main()
{

    vec3 vertex  = (gl_ModelViewMatrix * gl_Vertex).xyz;
    // Extract the normal from the normal map
//    vec3 bump = normalize(texture2D(NormalMap, gl_TexCoord[0].st).rgb * 2.0 - 1.0);

    vec3 n = normalize(gl_NormalMatrix * gl_Normal);
//    vec3 normal = (bump + n) / 2.0;
    vec3 i = normalize(vertex);


    Ratio   = F + (1.0 - F) * pow((1.0 - dot(-i, n)), FresnelPower);

    Refract = refract(i, n, Eta);
    RefractR = refract(i, n, EtaR);
    RefractG = refract(i, n, EtaG);
    RefractB = refract(i, n, EtaB);

    //I should probably have this sample from the terrain and not the skymap
    Refract = vec3(gl_TextureMatrix[0] * vec4(Refract, 1.0));
    RefractR = vec3(gl_TextureMatrix[0] * vec4(RefractR, 1.0));
    RefractG = vec3(gl_TextureMatrix[0] * vec4(RefractG, 1.0));
    RefractB = vec3(gl_TextureMatrix[0] * vec4(RefractB, 1.0));

    Reflect = reflect(i, n);
    Reflect = vec3(gl_TextureMatrix[0] * vec4(Reflect, 1.0));

    gl_Position = ftransform();
}
