#version 460

// In the beginning, colours never existed. There's nothing that can be done before you...

layout(location = 0) in vec2 in_uv;
layout(location = 0) out vec4 out_color;

layout(set = 0, binding = 0) uniform sampler2D previous_frame;
layout(set = 0, binding = 1) uniform sampler2D generic_texture;
layout(set = 0, binding = 2) uniform sampler2D dummy_texture;
#define T(t) (texture(t, vec2(in_uv.x, -in_uv.y)))
#define T_off(t,off) (texture(t, vec2(in_uv.x + off.x, -(in_uv.y + off.y))))

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

#define IVORY vec3(0.0, 0.3, 0.3)
#define WHITE vec3(1.0, 1.0, 1.0)
#define GOLDEN vec3(1.0, 0.7, 0.0)

#define PATTERN_RAND

const float PI = acos(-1.);
const float TWO_PI = 2. * PI;

mat2 rot(float a) {
	float c = cos(a), s = sin(a);
	return mat2(c, -s, s, c);
}

vec4 PatternRand( uint seed )
{
    return vec4(
        float((seed*0x73494U)&0xfffffU)/float(0x100000),
    	float((seed*0xAF71FU)&0xfffffU)/float(0x100000),
        float((seed*0x67a42U)&0xfffffU)/float(0x100000), // a bit stripey against x and z, but evens out over time
        float((seed*0x95a8cU)&0xfffffU)/float(0x100000) // good vs x or y, not good against z
    );
}

uvec4 Hash(uint seed) {
    // integer hash from Hugo Elias
    seed = (seed << 13U) ^ seed;
    seed = seed * (seed * seed * 15731U + 789221U) + 1376312589U;
    return seed * uvec4(seed, seed * 16807U, seed * 48271U, seed * 31713U);
}

vec4 Rand(uint seed) {
#if defined(PATTERN_RAND)
    return PatternRand(seed);
#else
    return vec4(Hash(seed) & 0x7fffffffU) / float(0x7fffffffU);
#endif
}

#define MAX_LEVEL 3
#define bayer2x2(a) (4-(a).x-((a).y<<1))%4
float GetBayerFromCoordLevel(vec2 pixelpos)
{
    ivec2 ppos = ivec2(pixelpos);
    int sum = 0;
    for(int i=0; i<MAX_LEVEL; i++)
    {
        sum += bayer2x2(ppos>>(MAX_LEVEL-1-i)&1)<<(2*i);
    }

    return float(sum) / float(2<<(MAX_LEVEL*2-1));
}

float scene(in vec2 uv, float time) {
    float t = mod(time, 3);
    uv /= 0.3;
    uv = rot(t * 30 * exp(-t * t * 2.5)) * uv;

    float pix = fwidth(6 * uv.y);
    /* pix = 1.0 / pc.resolution.y; */
    float res = smoothstep(0.0 - pix, 0 + pix,
                           length(uv * 1) * (abs(uv.y) * abs(uv.x)) + -0.30 +
                               -0.3 * pow(1 - t / 5, 3));
    /* smoothstep(0.0 - pix, 0 + pix,  (abs(uv.y) * abs(uv.x)) + -0.1 ); */
    res = max(res, length(uv) - 2.0);
    return res;
}

#define NUM_SAMPLES 1
void main() {
    vec2 uv = (in_uv + -0.5) * 2.0 / vec2(pc.resolution.y / pc.resolution.x, 1);
    vec2 frag_coord = in_uv*pc.resolution * 3;

    float motionJitter = GetBayerFromCoordLevel(in_uv*pc.resolution);
    /* motionJitter = Hash(uint(uv.x * uv.y)).x; */

    uint pseed = /*uint(pc.frame<<16)*/ + (uint(frag_coord.y)) * uint(frag_coord.x);
    vec4 prand = Rand( pseed );

/*     motionJitter = prand.y; */

    float res = 0;
    uint n = uint(NUM_SAMPLES);
    for (uint i = 0U; i < n; ++i) {
        float f = (float(i) + motionJitter) / float(n) - 0.5;
        float time = pc.time + pc.time_delta * f;
        res += scene(uv, time);
    }
    res /= float(n);

    /* res = scene(uv, pc.time); */

    vec3 col = mix(GOLDEN, IVORY, res);
    col = pow(col, vec3(2.4));
    out_color = vec4(col, 1.0);
}
