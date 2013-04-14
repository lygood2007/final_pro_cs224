const int MAX_KERNEL_SIZE = 128;
uniform sampler2D tex;
uniform int arraySize;
uniform vec2 offsets[MAX_KERNEL_SIZE]; 
uniform float kernel[MAX_KERNEL_SIZE];
void main(void) { 
    // TODO: Step 2 - Fill this in!
    vec4 sample = vec4(0.0);
    for(int i = 0; i < arraySize; i++) 
    {
	sample += texture2D(tex, gl_TexCoord[0].st + offsets[i])  * kernel[i];
    }
    gl_FragColor = sample;
}
