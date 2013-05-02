//attribute vec4 qt_Vertex;
//attribute vec4 qt_MultiTexCoord0;
//uniform mat4 qt_ModelViewProjectionMatrix;
//varying vec4 qt_TexCoord0;

//void main(void)
//{
//    gl_Position = qt_ModelViewProjectionMatrix * qt_Vertex;
//    qt_TexCoord0 = qt_MultiTexCoord0;
//}

uniform vec2 windowSize;
varying vec2 screenPos;
varying float radius;

void main() {
  gl_Position = gl_ModelViewProjectionMatrix * gl_Vertex;
  gl_PointSize = 100.0;

  // Convert position to window coordinates
  vec2 halfsize = vec2(windowSize.x, windowSize.y) * 0.5;
  screenPos = halfsize + ((gl_Position.xy / gl_Position.w) * halfsize);

  // Convert radius to window coordinates
  radius = gl_PointSize * 0.5;
}
