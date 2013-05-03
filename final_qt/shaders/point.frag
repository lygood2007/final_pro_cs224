uniform samplerCube CubeMap; //the skymap
uniform vec4 CurrColor; //the current color of the water
varying vec3 n, lightDir, r;

void main (void)
{
//	vec4 current_color = CurrColor;
//        vec4 final_color = textureCube( CubeMap, r);
        vec4 final_color = vec4(1.0,1.0,1.0,1.0);
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
        gl_FragColor = mix(CurrColor, final_color,0.85);
//        gl_FragColor = final_color;
}
