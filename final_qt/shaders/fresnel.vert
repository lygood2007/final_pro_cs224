
varying vec3 vertex;		// The position of the vertex, in eye space
varying vec3 eye;		// The normalized vector from the vertex to the eye
varying vec3 normal;		// The normal vector of the vertex, in eye space

void main() 
{
    vertex = (gl_ModelViewMatrix * gl_Vertex).xyz;
    normal = normalize(gl_NormalMatrix * gl_Normal);
    eye = vec3(gl_ProjectionMatrixInverse*vec4(0,0,-1,0)).xyz;

    gl_Position = ftransform();
}
