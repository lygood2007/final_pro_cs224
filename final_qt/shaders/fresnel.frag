
varying vec3 vertex;		// The position of the vertex, in eye space
varying vec3 light;		// The normalized vector from the vertex to the light
varying vec3 eye;		// The normalized vector from the vertex to the eye
varying vec3 normal;		// The normal vector of the vertex, in eye space

uniform float r0;		// The R0 value to use in Schlick's approximation
uniform float rS;		// The R0 value to use in Schlick's approximation
uniform float eta1D;		// The eta value to use initially
uniform vec3  eta;		// Contains one eta for each channel (ue eta.r, eta.g, eta.b in your code)

// IoR of Air is   1.0002926
// IoR of Water is 1.33335  (at 20 deg C)
// (n1-n2)/(n1+n2) = -0.1427199692

uniform samplerCube CubeMap;
uniform vec4 CurrColor; //the current color of the water

void main()
{
    vec3 n = normalize(normal);
    vec3 e = normalize(eye);
    vec3 i = normalize(vertex - eye);
    
    vec3 incident = reflect(i,n);
    incident = (gl_ModelViewMatrixInverse * vec4(incident, 0.0)).xyz;
    vec4 reflectColor = textureCube(CubeMap, incident);
    
//    vec3 eta = vec3(0.8,0.8,0.8);

    vec3 vRr = refract(i, n, eta.r);
    vec3 vRg = refract(i, n, eta.g);
    vec3 vRb = refract(i, n, eta.b);
    
    vRr = (gl_ModelViewMatrixInverse * vec4(vRr, 0.0)).xyz;
    vRg = (gl_ModelViewMatrixInverse * vec4(vRg, 0.0)).xyz;
    vRb = (gl_ModelViewMatrixInverse * vec4(vRb, 0.0)).xyz;
      
    //fresnel equation
//    float Rs = 0.1427199692;
    float r0 = 0.4;
    float costheta = dot(-i,n);
    float f = rS + (1.0 - rS) * pow(1.0 - costheta,5.0);
    float t = 1.0 - f;
    
//      gl_FragColor = CurrColor;
//        gl_FragColor = reflectColor;
//            gl_FragColor = CurrColor * 0.3 + reflectColor * 0.7;


//        gl_FragColor =vec4((textureCube(CubeMap,vRr).r * t) +(f * reflectColor.r),
//                            (textureCube(CubeMap,vRg).g * t) +(f * reflectColor.g),
//                            (textureCube(CubeMap,vRb).b * t) +(f * reflectColor.b), 1.0);

//    reflectColor = vec4(1,0,0,1);
//            gl_FragColor =vec4((textureCube(CubeMap,vRr).r * t) +(f * reflectColor.r),
//                               (textureCube(CubeMap,vRg).g * t) +(f * reflectColor.g),
//                                (textureCube(CubeMap,vRb).b * t) +(f * reflectColor.b), 1.0);

//    gl_FragColor = vec4(costheta,costheta,costheta,1);
//    gl_FragColor = vec4(f, f, f, 1);
    gl_FragColor = CurrColor * 0.5 +
                    vec4((textureCube(CubeMap,vRr).r * t) +(f * reflectColor.r),
                        (textureCube(CubeMap,vRg).g * t) +(f * reflectColor.g),
                        (textureCube(CubeMap,vRb).b * t) +(f * reflectColor.b), 1.0) * 0.5;

}
