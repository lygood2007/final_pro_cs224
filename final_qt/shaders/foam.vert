
varying vec3 n, lightDir, r;
const vec3 L = vec3(0.,0.,1.);

void main()
{
        vec3 vertex = vec3(gl_ModelViewMatrix * gl_Vertex).xyz;
        n = normalize( gl_NormalMatrix * gl_Normal );
        gl_PointCoord[0] = gl_MultiTexCoord0;
        vec3 i = normalize(vertex);
        lightDir = vec3(L - vVertex);

         r = reflect(I,normal);

    gl_Position = ftransform();
}
