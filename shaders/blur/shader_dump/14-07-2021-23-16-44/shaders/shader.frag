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

#define PI acos(-1.)
/* const float time_delta = 8; */

mat2 rot(float a) {
	float c = cos(a), s = sin(a);
	return mat2(c,-s,s,c);
}

float bayer8(ivec2 uv) {
	uv %= 8;
	return texture(float_texture1, in_uv).r;
}

float ramp(float t) {
	t = smoothstep(0.,1.,t);
	t = smoothstep(0.,1.,t);
	return t;
}

float sdf(vec2 p, float time) {
	p *= 2;
	p.y += 0.8;
	p *= rot(PI * (4.17 ));

	const float tau = acos(-1.) * 2.;
	const float duration = 5.;
	const float revolutions = 5;
	float t = ramp(mod(time, duration) / duration) * revolutions;
	float angle = t*tau;

	float a = 2, b = 2, k = 3;
	float r = a + b*cos(k*angle);
	vec2 polars = r * vec2(cos(angle), sin(angle));
	return length(p - polars) - .5;
}

void main() {
    vec2 uv = (in_uv + -0.5) * 2.0 / vec2(pc.resolution.y / pc.resolution.x, 1);
	uv *= 2.;

    float edge = dFdx(uv.x) * 0.5;

    float c = 0.;
    const int steps = 3;
    for (int i = 0; i < steps; ++i) {
        float subsample = bayer8(ivec2(in_uv * pc.resolution));
        float time =
            pc.time + ((float(i) + subsample) / float(steps) - 0.5) * pc.time_delta;
        c += smoothstep(-edge, edge, sdf(uv, time));
    }
    c /= float(steps);

    c = pow(c, .4545);
    /* c = bayer8(ivec2(in_uv)); */

    out_color = vec4(vec3(c), 1.0);
}
