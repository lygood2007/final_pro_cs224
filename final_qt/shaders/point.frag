//uniform sampler2D qt_Texture0;
//varying vec4 qt_TexCoord0;

//void main(void)
//{
//    gl_FragColor = texture2D(qt_Texture0, qt_TexCoord0.st);
//}
varying vec2 screenPos;
varying float radius;

void main() {

  // Sphere shaded
  if( distance(gl_FragCoord.xy, screenPos) > radius ) discard;
  gl_FragData[0] = gl_Color;

}
