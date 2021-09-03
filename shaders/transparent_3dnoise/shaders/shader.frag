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
} pc;

float smax( float a, float b, float k ) {
    k *= 1.4;
    float h = max(k-abs(a-b),0.0);
    return max(a, b) + h*h*h/(6.0*k*k);
}

#define MAX_LEVEL 1
#define bayer2x2(a) (4-(a).x-((a).y<<1))%4
float GetBayerFromCoordLevel(vec2 pixelpos) {
    ivec2 ppos = ivec2(pixelpos);
    int sum = 0;
    for(int i=0; i<MAX_LEVEL; i++)
    {
        sum += bayer2x2(ppos>>(MAX_LEVEL-1-i)&1)<<(2*i);
    }

    return float(sum) / float(2<<(MAX_LEVEL*2-1));
}

#define PATTERN_RAND
vec4 PatternRand( uint seed ) {
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

vec3 hash33(vec3 p) {
	float n = sin(dot(p, vec3(7, 157, 113)));
	return fract(vec3(2097152, 262144, 32768) * n);
}

float noise3D(vec3 p) {
	p *= 1.1;
	const vec3 s = vec3(15, 157, 113);

	vec3 ip = floor(p);

	vec4 h = vec4(0., s.yz, s.y + s.z) + dot(ip, s);

	p -= ip;

	p = p*p*(3. - 2.*p);

	p = p*p*p*(p*(p * 6. - 15.) + 10.);

	h = mix(fract(sin(h)*43758.5453), fract(sin(h + s.x)*43758.5453), p.x);

	h.xy = mix(h.xz, h.yw, p.y);

	return mix(h.x, h.y, p.z);
}

vec3 ro = vec3(0., pc.time * 1.5, pc.time * 1.5);
float map(vec3 p) {
	float or = length(p - ro) - 1.5;
	float noise = noise3D(p * 2.) - .3;
	/* return noise; */
	return smax(-or, noise, -1.1);
	/* return noise3D(p*2.)*.66 + noise3D(p*4.)*.34 - .4; */
}

void main() {
    vec2 uv = (in_uv + -0.5) * 2.0 / vec2(pc.resolution.y / pc.resolution.x, 1);

    vec3 rd = normalize(vec3(uv, (1. - dot(uv, uv) * 0.5) * 0.5));
    rd = normalize(vec3(uv, 1.));
    ro = vec3(0., pc.time * 1.5, pc.time * 1.5);
    /* ro = vec3(0., 0., -5.); */
    vec3 col = vec3(0.);
    vec3 sp = vec3(0.);

    vec2 a = sin(vec2(1.5707963, 0.) + pc.time * 0.375);
    rd.xz = mat2(a, -a.y, a.x) * rd.xz;
    rd.xy = mat2(a, -a.y, a.x) * rd.xy;

    rd *= 0.99 + hash33(rd) * 0.02;

    vec3 rnd = hash33(rd + 311.);

    float t = 0., layers = 0., d = 0.0001, aD = 0.;


    t = .2 * length(hash33(rd));
    {
        vec2 frag_coord = in_uv * pc.resolution;

        float motionJitter = GetBayerFromCoordLevel(frag_coord);
        /* motionJitter = Hash(uint(uv.x * uv.y)).x; */

        uint pseed =
            uint(pc.frame<<16) +(uint(frag_coord.y)) * uint(frag_coord.x);
        vec4 prand = Rand(pseed);

        /* motionJitter = prand.y; */

        t /= motionJitter;
    }

    float thD = .025;// + smoothstep(-0.2, 0.2, sin(pc.time*0.75 - 3.141592 * 0.4)) * 0.025;

    for (float i = 0.; i < 53.; ++i) {
        if (layers > 23. || col.x > 1. || t > 10.)
            break;

        sp = ro + rd * t;

        d = map(sp);

        aD = (thD - abs(d) * 23. / 24.) / thD;

        if (aD > 0.) {
            col += aD*aD*(3.-2.*aD)/(1. + t*t*0.125)*.1 + (fract(rnd + i*27.)-.5)*0.02;
            /* col += aD / (1. + t * t * 0.1) * 0.1 + */
            /*        (fract(rnd + i * 27.) - .5) * 0.01; */

            layers += 1;
        }

        t += max(abs(d) * 0.7, thD * 0.7);
    }

    col = max(col, 0.);

	uv = abs(in_uv * vec2(pc.resolution.y / pc.resolution.x, 1.) - 0.5);
    col = mix(col, vec3(min(col.x * 1.5, 1.), pow(col.x, 2.5), pow(col.x, 12.)),
              min(dot(pow(uv, vec2(5.8)), vec2(1.)) * 8., 1.));

    col = mix(col, col.xzy, dot(sin(rd * 5.), vec3(.166)) + 0.166);

	col = pow(col, vec3(1.4));

    out_color = vec4(col, 1.0);
}
