//uniform sampler2D circle_texture;
uniform sampler2D texture;
uniform vec4 CurrColor; //the current color of the water
varying vec3 n, lightDir, r;

void main (void)
{
        vec4 final_color = vec4(1.0,1.0,1.0,0.2);
//        final_color = texture2D(circle_texture, gl_PointCoord[0].st);
        final_color = texture2D(texture, gl_TexCoord[0].st);
        vec3 N = normalize(n);
        vec3 L = normalize(lightDir);
        float lambertTerm = dot(N,L);
        if(lambertTerm > 0.0)
        {
                // Specular
                final_color += vec4(1.0,1.0,1.0,0.2); //textureCube( CubeMap, r);
        }
        gl_FragColor = mix(CurrColor, final_color,0.85);
}
