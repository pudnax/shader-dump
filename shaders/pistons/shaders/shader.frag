#version 460

// In the beginning, colours never existed. There's nothing that can be done before you...

#include <prelude.glsl>

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

float scene(vec2 p) {
	return length(p) - 0.5;
}

vec3 shade_dist(float d) {
	float dist = d * 100.;
	float banding = sin(dist);
	float strength = sqrt(1. - exp(-abs(d) * 2.));
	float pattern = mix(strength, banding, (0.6 - abs(strength - 0.5)) * 0.3);

	vec3 color = vec3(pattern);

	color *= d > 0.0 ? vec3(1.0, 0.56, 0.4) : vec3(0.4, 0.9, 1.0);

	return color;
}

void main() {
	vec2 aspect = vec2(pc.resolution.y / pc.resolution.x, 1);
    vec2 uv = (in_uv + -0.5) * 2.0 / aspect;
	vec2 mouse = pc.mouse / aspect;
	uv *= 5.;
	mouse *= 5.0;

	float mouse_dist = scene(mouse);
    vec3 col = shade_dist(scene(uv));

	if (distance(mouse, uv) < abs(mouse_dist) && pc.mouse_pressed)
		col *= 0.5;

	col = pow(col, vec3(2.4));
    out_color = vec4(col, 1.0);
}
