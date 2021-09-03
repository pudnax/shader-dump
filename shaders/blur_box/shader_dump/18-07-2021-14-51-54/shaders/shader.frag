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

#define PI acos(-1.)
#define rot(a) mat2(cos(a+PI*0.5*vec4(0,1,3,0)))

float hash13(vec3 p3) {
    p3 = fract(p3 * 0.1031);
    p3 += dot(p3, p3.yzx + 19.19);
    return fract((p3.x + p3.y) * p3.z);
}

vec3 scene(vec2 uv, float t) {
	uv *= 1.4;
	uv *= rot(t * 10. + (sin(t * 2.0) *0.5 + 0.5) *10.);
	uv = abs(uv);
	float sd = max(uv.x - 0.5, uv.y - 1.5);
	return vec3(smoothstep(0.0, 0.4, sd));
}

void main() {
    vec2 uv = (in_uv + -0.5) * 2.0 / vec2(pc.resolution.y / pc.resolution.x, 1);

	vec3 res = vec3(0.);

	bool motionBlur = true;
	/* motionBlur = false; */
	if (motionBlur) {
#define BLUR 30
		for (int i = 0; i < BLUR; i++) {
			float rnd = hash13(vec3(uv, pc.frame*100+i));
			float time = pc.time + rnd/60.0;
			res += scene(uv, time);
		}
		res /= float(BLUR);
	} else {
		res = scene(uv, pc.time);
	}

    vec3 col = vec3(res);
	col = pow(col, vec3(1.0/2.2));
    out_color = vec4(col, 1.0);
}
