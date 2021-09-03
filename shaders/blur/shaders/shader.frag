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

float bayer8(void) {
	ivec2 uv = ivec2(in_uv * pc.resolution);
	uv %= 8;
	return texelFetch(float_texture1, uv, 0).r;
}

float ramp (float t) {
	t = smoothstep(0., 1., t);
	t = smoothstep(0., 1., t);
	return t;
}

vec4 PRand( uint seed )
{
    return vec4(
        float((seed*0x73494U)&0xfffffU)/float(0x100000),
    	float((seed*0xAF71FU)&0xfffffU)/float(0x100000),
        float((seed*0x67a42U)&0xfffffU)/float(0x100000), // a bit stripey against x and z, but evens out over time
        float((seed*0x95a8cU)&0xfffffU)/float(0x100000) // good vs x or y, not good against z
        );
}

float sdf(vec2 p, float time) {
	const float tau = acos(-1.) * 2.;
	const float duration = 5.;
	const float revolutions = 40.;
	float t = ramp(mod(time, duration) / duration) * revolutions;
	float angle = t * tau;
	return length(p - vec2(sin(angle), cos(angle))) - 0.2;
}

void main() {
    vec2 uv = (in_uv + -0.5) * 2.0 / vec2(pc.resolution.y / pc.resolution.x, 1);
    uv *= 2.;

    float edge = fwidth(uv.x) * 0.5;

	vec2 fragCoord = in_uv*pc.resolution;

    uint pseed =
        uint(int(floor(pc.time)) << 16) + (uint(fragCoord.y) << 8U) + uint(fragCoord.x);

    vec4 prand = PRand(pseed);

    float c = 0.;
    const int steps = 4;
    for (int i = 0; i < steps; ++i) {
        float subsample = prand.y;
        float time = pc.time + ((float(i) + subsample) / float(steps) - 0.5) *
                                   pc.time_delta;
        c += smoothstep(-edge, edge, sdf(uv, time));
    }
    c /= float(steps);

    c = pow(c, 2.4545);

    out_color = vec4(vec3(c), 1.0);
}
