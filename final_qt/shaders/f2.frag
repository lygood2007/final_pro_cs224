varying vec3  Reflect;
varying vec3  Refract;
varying float Ratio;

uniform samplerCube Cubemap;
uniform vec4 CurrColor;

void main()
{
    vec3 refractColor = vec3(textureCube(Cubemap, Refract));
    vec3 reflectColor = vec3(textureCube(Cubemap, Reflect));

//    refractColor = vec3(1,0,1);
//    reflectColor = vec3(0,0,1);

    vec3 color   = mix(refractColor, reflectColor, Ratio);

    gl_FragColor = mix(CurrColor, vec4(color, 1.0), .6);
//    gl_FragColor = vec4(color, 1.0);
//    gl_FragColor = CurrColor;

//    gl_FragColor = vec4(reflectColor, 1.0);
}
