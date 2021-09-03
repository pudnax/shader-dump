#version 460

// In the beginning, colours never existed. There's nothing that can be done
// before you...

layout(location = 0) in vec2 in_uv;
layout(location = 0) out vec4 out_color;

layout(set = 0, binding = 0) uniform sampler2D previous_frame;
layout(set = 0, binding = 1) uniform sampler2D generic_texture;
layout(set = 0, binding = 2) uniform sampler2D dummy_texture;
#define T(t) (texture(t, vec2(in_uv.x, -in_uv.y)))
#define T_off(t, off) (texture(t, vec2(in_uv.x + off.x, -(in_uv.y + off.y))))

layout(set = 0, binding = 3) uniform sampler2D float_texture1;
layout(set = 0, binding = 4) uniform sampler2D float_texture2;

layout(set = 1, binding = 0) uniform sampler1D fft_texture;

layout(std430, push_constant) uniform PushConstant {
    vec3 pos;
    float time;
    vec2 resolution;
    vec2 mouse;
    bool mouse_pressed;
    uint frame;
    float time_delta;
} pc;

#define blur 1.

#define DITHER
#define QUALITY 2

#define DECAY .974
#define EXPOSURE .24
#if (QUALITY == 2)
#define SAMPLES 64
#define DENSITY .97
#define WEIGHT .25
#else
#if (QUALITY == 1)
#define SAMPLES 32
#define DENSITY .95
#define WEIGHT .25
#else
#define SAMPLES 16
#define DENSITY .93
#define WEIGHT .36
#endif
#endif

#define ViewZoom 3.
#define fsaa 14. / min(pc.resolution.x, pc.resolution.y)
#define fra(u) (u - 0.5 * pc.resolution.xy) * ViewZoom / pc.resolution.y

#define iterBayerMat 1
#define bayer2x2(a) (4 - (a).x - ((a).y << 1)) % 4

float GetBayerFromCoordLevel(vec2 pixelpos) {
    ivec2 p = ivec2(pixelpos);
    int a = 0;
    for (int i = 0; i < iterBayerMat; i++) {
        a += bayer2x2(p >> (iterBayerMat - 1 - i) & 1) << (2 * i);
    }
    return float(a) / float(2 << (iterBayerMat * 2 - 1));
}
// https://www.shadertoy.com/view/XtV3RG

float bayer2  (vec2 a){a=floor(a);return fract(dot(a,vec2(.5, a.y*.75)));}
float bayer4  (vec2 a){return bayer2 (  .5*a)*.25    +bayer2(a);}
float bayer8  (vec2 a){return bayer4 (  .5*a)*.25    +bayer2(a);}
float bayer16 (vec2 a){return bayer4 ( .25*a)*.0625  +bayer4(a);}
float bayer32 (vec2 a){return bayer8 ( .25*a)*.0625  +bayer4(a);}
float bayer64 (vec2 a){return bayer8 (.125*a)*.015625+bayer8(a);}
float bayer128(vec2 a){return bayer16(.125*a)*.015625+bayer8(a);}
#define dither2(p)   (bayer2(  p)-.375      )
#define dither4(p)   (bayer4(  p)-.46875    )
#define dither8(p)   (bayer8(  p)-.4921875  )
#define dither16(p)  (bayer16( p)-.498046875)
#define dither32(p)  (bayer32( p)-.499511719)
#define dither64(p)  (bayer64( p)-.49987793 )
#define dither128(p) (bayer128(p)-.499969482)


float iib(vec2 u) {
    return dither128(u);
    return GetBayerFromCoordLevel(u*999.);
}

vec3 sun(vec2 uv) {
    vec2 p = fra(abs(vec2(pc.mouse_pressed))) * .666;
    if (!pc.mouse_pressed) p = vec2(sin(pc.time), sin(pc.time * 0.5) * 0.5);
    vec3 res;
    float di = distance(uv, p);
    res.x = di <= .3333 ? sqrt(1. - di * 3.0) : 0.0;
    res.yz = p;
    res.y /= (pc.resolution.x / pc.resolution.y);
    res.yz = (res.yz + 1.) * 0.5;
    return res;
}

#define SS blur/min(pc.resolution.x, pc.resolution.y)

float circle(vec2 p, float r) {
    return smoothstep(SS, -SS, length(p) - r);
}

vec4 BA(in vec2 uv) {
    uv = uv * 2. - 1.;
    float aspect = pc.resolution.x / pc.resolution.y;
    uv.x *= aspect;
    vec2 m = pc.mouse;
    if (m.x <= 0.) m = vec2(pc.resolution * 0.1);
    else m = (m / pc.resolution.xy) * 2. - 1.;
    m.x *= aspect;

    float occluders = circle(uv - vec2(-0.66, 0.), 0.366) -
                      circle(uv + vec2(0.75, 0.1), 0.18);
    occluders += circle(uv - vec2(0.6, 0.2), 0.23);
    float mouse =
        smoothstep(SS, -SS, abs(abs(length(uv - m) - 0.2) - 0.05) - 0.2);
    float mouse2 = smoothstep(SS, -SS, abs(uv.x - m.x) - .03);  // vertical bar
    mouse = max(mouse, mouse2);                                 // union
    mouse2 = smoothstep(SS, -SS, length(uv - m - vec2(.5, 0)) - .2);
    mouse = max(mouse, mouse2 * .5);  // union

    occluders += mouse;
    occluders = min(occluders, 1.);
    vec3 light = min(sun(uv), 1.);
    float col = max(light.x - occluders, 0.);
    return vec4(col, occluders, light.yz);
}

void main() {
    vec2 uv = (in_uv + -0.5) * 2.0 / vec2(pc.resolution.y / pc.resolution.x, 1);
    uv = in_uv / vec2(pc.resolution.y / pc.resolution.x, 1);
    vec2 coord = uv;
    vec4 ic = BA(uv);
    vec2 lightpos = ic.zw;
    float occ = ic.x;
    float obj = ic.y;
    float dither = iib(in_uv*pc.resolution);
    vec2 dtc = (coord - lightpos)*(1. / float(SAMPLES)*DENSITY);
    float illumdecay = 1.;

    for (int i = 0; i < SAMPLES; ++i) {
        coord -= dtc;
#ifdef DITHER
        float s = BA(coord + (dtc * dither)).x;
#else
        float s = BA(coord).x;
#endif
        s *= illumdecay * WEIGHT;
        occ += s;
        illumdecay *= DECAY;
    }

    vec3 col = vec3(uv, 1.);
    col = vec3(0.5,0.7,0.1)*obj/3. + occ*EXPOSURE;
    out_color = vec4(col, 1.0);
}
