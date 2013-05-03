//attribute vec4 qt_Vertex;
//attribute vec4 qt_MultiTexCoord0;
//uniform mat4 qt_ModelViewProjectionMatrix;
//varying vec4 qt_TexCoord0;

//void main(void)
//{
//    gl_Position = qt_ModelViewProjectionMatrix * qt_Vertex;
//    qt_TexCoord0 = qt_MultiTexCoord0;
//}


varying vec3 n, lightDir, r;
const vec3 L = vec3(0.,0.,1.);

void main()
{
        vec3 vertex = vec3(gl_ModelViewMatrix * gl_Vertex).xyz;
        n = normalize( gl_NormalMatrix * gl_Normal );
        vec3 i = normalize(vertex);
        lightDir = vec3(L - vVertex);

         r = reflect(I,normal);


    gl_Position = ftransform();
}

/*
 * Vertex shader for chromatic aberration effect.
 *
 * Author: Randi Rost
 *
 * Copyright (c) 2003-2006: 3Dlabs, Inc.
 * See 3Dlabs-License.txt for license information
 *
 * Modified from the original example and combined with phong shading
 * by Karsten Schmidt.
 */

//const float F = 0.05;

//uniform float etaR;
//uniform float etaG;
//uniform float etaB;
//uniform float fresnelPower;

//varying vec3  reflect;
//varying vec3  refractR;
//varying vec3  refractG;
//varying vec3  refractB;
//varying float ratio;

//varying vec3 normal, lightDir, eyeVec;

//void main() {
//    vec4 ecPosition  = gl_ModelViewMatrix * gl_Vertex;
//    vec3 ecPosition3 = ecPosition.xyz / ecPosition.w;

//    vec3 i = normalize(ecPosition3);
//    vec3 n = normalize(gl_NormalMatrix * gl_Normal);

//    ratio   = F + (1.0 - F) * pow(1.0 - dot(-i, n), fresnelPower);

//    refractR = refract(i, n, etaR);
//    refractR = vec3(gl_TextureMatrix[0] * vec4(refractR, 1.0));

//    refractG = refract(i, n, etaG);
//    refractG = vec3(gl_TextureMatrix[0] * vec4(refractG, 1.0));

//    refractB = refract(i, n, etaB);
//    refractB = vec3(gl_TextureMatrix[0] * vec4(refractB, 1.0));

//    reflect  = reflect(i, n);
//    reflect  = vec3(gl_TextureMatrix[0] * vec4(reflect, 1.0));

//    normal = gl_NormalMatrix * gl_Normal;

//    vec3 vVertex = vec3(gl_ModelViewMatrix * gl_Vertex);

//    lightDir = vec3(gl_LightSource[0].position.xyz - vVertex);
//    eyeVec = -vVertex;

//    gl_Position = ftransform();
//}
