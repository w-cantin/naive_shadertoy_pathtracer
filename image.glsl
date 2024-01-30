

void mainImage(out vec4 fragColor,in vec2 fragCoord)
{
    vec2 uv=fragCoord/iResolution.xy;
    
    vec4 colorChannel0=texture(iChannel0,uv).rgba;
    
    //The other shader accumulates the total amount of luminance collected
    // by the multiple pass. We store the amount of samples in the alpha channel
    // and then calculate the average of the luminance / number of samples
    //for each color
    fragColor=vec4(colorChannel0/colorChannel0.a);
    
    //We apply a gamma correction to
    // For more info: https://learnopengl.com/Advanced-Lighting/Gamma-Correction
    const float GAMMA=2.2;
    fragColor.rgb=pow(fragColor.rgb,vec3(1./GAMMA));
    
}
