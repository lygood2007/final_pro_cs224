uniform sampler2D tex;
const vec3 avgVector = vec3(0.299, 0.587, 0.114);
void main(void) {
    vec4 sample = texture2D(tex, gl_TexCoord[0].st);
    float luminance = max(0.0, dot(avgVector, sample.rgb));

    // set fragment color to black if the luminance is less than some number
    //adjust these values until you like the effect you are getting
    if(luminance < 1.0) {
	sample.rgb = vec3(0.0, 0.0, 0.0);
    } else {
        luminance = luminance/(luminance + 1.0 );
	sample.rgb = vec3(luminance, luminance, luminance);
    }
    gl_FragColor = sample;
}
