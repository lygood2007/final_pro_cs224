uniform samplerCube CubeMap; //the skymap
uniform vec4 CurrColor; //the current color of the water
varying vec3 n, lightDir, r;

void main (void)
{
//	vec4 current_color = CurrColor;
//        vec4 final_color = textureCube( CubeMap, r);
        vec4 final_color = vec4(1.0,1.0,1.0,0.5);
        vec3 N = normalize(n);
        vec3 L = normalize(lightDir);
        float lambertTerm = dot(N,L);
        if(lambertTerm > 0.0)
        {
                // Specular
                final_color += textureCube( CubeMap, r);
//            final_color += vec4(1.0, 1.0, 1.0, 1.0);
        }
//                gl_FragColor = CurrColor;
//        gl_FragColor = CurrColor *(0.5) +  final_color * (0.5);
        gl_FragColor = mix(CurrColor, final_color,0.90);
//        gl_FragColor = final_color;
}


//varying vec3  Reflect;
//varying vec3  Refract;
//varying vec3  RefractR;
//varying vec3  RefractG;
//varying vec3  RefractB;
//varying float Ratio;

//uniform samplerCube Cubemap;
//uniform vec4 CurrColor;

//void main()
//{
////    vec3 refractColor = vec3(textureCube(Cubemap, Refract));

//    vec3 refractColor = vec3(vec3(textureCube(Cubemap, RefractR)).r,
//                        vec3(textureCube(Cubemap, RefractG)).g,
//                        vec3(textureCube(Cubemap, RefractB)).b);


//    vec3 reflectColor = vec3(textureCube(Cubemap, Reflect));

//    vec4 light_color = vec4(1.0,1.0,1.0,0.5); //instead of using watercolor use bright white

//    vec3 color   = mix(refractColor, reflectColor, Ratio);

//    gl_FragColor = mix(light_color, vec4(color, 0.9), .35);
////    gl_FragColor = vec4(color, 1.0);
////    gl_FragColor = CurrColor;

////    gl_FragColor = vec4(reflectColor, 1.0);
//}
