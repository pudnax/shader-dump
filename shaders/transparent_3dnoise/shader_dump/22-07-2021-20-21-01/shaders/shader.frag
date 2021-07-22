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

	h = mix(fract(sin(h)*43758.5453), fract(sin(h + s.x)*43758.5453), p.x);

	/* h.xy = mix(h.xz, h.yw, p.y); */

	return mix(h.x, h.y, p.z);
}

vec3 ro = vec3(0., pc.time * 1.5, pc.time * 1.5);
float map(vec3 p) {
	float or = length(p - ro) - 3.5;
	float noise = noise3D(p * 2.) - .3;
	/* return noise; */
	return max(-or, noise);
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

    float thD = .025;  // + smoothstep(-0.2, 0.2, sin(pc.time*0.75
                       // - 3.141592 * 0.4)) * 0.025;

    for (float i = 0.; i < 53.; ++i) {
        if (layers > 23. || col.x > 1. || t > 10.)
            break;

        sp = ro + rd * t;

        d = map(sp);

        aD = (thD - abs(d) * 23. / 24.) / thD;

        if (aD > 0.) {
            col += aD / (1. + t * t * 0.1) * 0.1 +
                   (fract(rnd + i * 27.) - .5) * 0.01;

            layers += 1;
        }

        t += max(abs(d) * 0.7, thD * 0.7);
    }

    col = max(col, 0.);

	uv = abs(in_uv * vec2(pc.resolution.y / pc.resolution.x, 1.) - 0.5);
    /* col = mix(col, vec3(min(col.x * 1.5, 1.), pow(col.x, 2.5), pow(col.x, 12.)), */
    /*           min(dot(pow(uv, vec2(4.)), vec2(1.)) * 8., 1.)); */

	/* rd = normalize(vec3(uv, 1.)); */
	/* d = 0.; */
	/* for (int i = 0; i < 100; i++) { */
	/* 	d += map(ro + rd * d); */
	/* } */
	/* col = vec3(fract(d)); */

    col = mix(col, col.xzy, dot(sin(rd * 5.), vec3(.166)) + 0.166);

	col = pow(col, vec3(1.74));

    out_color = vec4(col, 1.0);
}
