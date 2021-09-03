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

#define WHITE vec3(1.0, 1.0, 1.0)
#define SWEETPEA vec3(1.0, 0.7, 0.75)
#define sat(x) clamp(x, 0., 1.)
#define ASPECT vec2(pc.resolution.y / pc.resolution.x, 1.)

float linearstep(float begin, float end, float t) {
    return sat((t - begin) / (end - begin));
}

float easeInOutCubic(float t) {
	if ((t *= 2.0) < 1.0) {
		return 0.5 * t * t * t;
	} else{
		return 0.5 * ((t -= 2.0) * t * t + 2.0);
	}
}

float easeInOutExpo(float t) {
	if (t == 0.0 || t == 1.0) {
		return t;
	}
	if ((t *= 2.0) < 1.0) {
		return 0.5 * pow(2.0, 10.0 * (t - 1.0));
	} else {
		return 0.5 * (-pow(2.0, -10.0 * (t - 1.0)) + 2.0);
	}
}

float smoothedge(float v, float f) {
	return smoothstep(0.0, f / pc.resolution.x, v);
}

float circle(vec2 p, float radius) {
	return length(p) - radius;
}

float dotPlot(vec2 uv, vec2 p) {
	return 1.0 - smoothedge(circle(uv - p / ASPECT, 0.05), 1.0);
}

void main() {
    vec2 st = (in_uv + -0.5) * 2.0 / vec2(pc.resolution.y / pc.resolution.x, 1.);
	st = in_uv / vec2(pc.resolution.y / pc.resolution.x, 1.);

	float t = fract(pc.time / 2.0), v = 0;

	float t0 = linearstep(0.2, 0.6, t);
	float p0 = easeInOutCubic(t0);
	v = dotPlot(st, vec2(mix(0.2, 0.8, p0), 0.2));

	float t1 = linearstep(0.1, 0.5, t);
	float p1 = easeInOutCubic(t1);
	v = max(v, dotPlot(st, vec2(mix(0.2, 0.8, p1), 0.6)));

	float t2 = linearstep(0.1, 0.5, t);
	float p2 = easeInOutCubic(t2);
	float t3 = linearstep(0.6, 1.0, t);
	float p3 = easeInOutCubic(t3);
	v = max(v, dotPlot(st, vec2(mix(0.2, 0.8, p2 - p3), 0.4)));

	float t4 = linearstep(0.1, 0.5, t);
	float p4 = easeInOutExpo(t4);
	float t5 = linearstep(0.6, 1.0, t);
	float p5 = easeInOutExpo(t5);
	v = max(v, dotPlot(st, vec2(mix(0.2 ,0.8, p4 - p5), 0.8)));

    vec3 col = mix(SWEETPEA, WHITE, v);
	col = pow(col, vec3(2.4));
    out_color = vec4(col, 1.0);
}
