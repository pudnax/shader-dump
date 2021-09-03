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

#define SIZE 25.
#define R_FACTOR 0.8

#define BLACK_OIL vec3(41,55,66)/255.
#define WHITE_OIL vec3(245,248,250)/255.

#define AAstep(thre, val) smoothstep(-.7,.7,(val-thre)/min(0.07,fwidth(val-thre)))

float rand(vec2 co) {
	return fract(sin(dot(co.xy, vec2(12.9898, 78.233))) * 43758.5454);
}

void main() {
    vec2 cuv =
        (in_uv + -0.5) * 2.0 / vec2(pc.resolution.y / pc.resolution.x, 1);

    float t = pc.time * 0.1;
    float ac = cos(t);
    float as = sin(t);
    mat2 rot = mat2(ac, -as, as, ac);
    cuv *= rot;

    float fovThete = 7.55;
    float z = sqrt(0.5 - cuv.x * cuv.x - cuv.y * cuv.y);
    float a = 1. / (z * tan(fovThete * 0.5));
    vec2 uv = cuv * a;

    vec2 ruv = uv * SIZE;
    vec2 id = ceil(ruv);

    ruv.x -= t * 50. * (rand(vec2(id.y)) * 0.5 + 0.5);
    vec2 guv = fract(ruv) - 0.5;

    id = ceil(ruv);

    float m = step(rand(id), R_FACTOR);
    m *= ceil(mod(id.y, 2.));

    for (float i = -1.; i <= 1.; i += 2.) {
        float mn = step(rand(id + vec2(i, 0.)), R_FACTOR);

        float l = length(guv + vec2(0.5, 0.) * i) * +(1. - mn);
        m *= step(l, 0.5);
    }

	float ll = length(cuv);
	m *= step(ll, 0.5);

    vec3 col = vec3(uv, 1.);
	col = mix(WHITE_OIL, BLACK_OIL, m);
    out_color = vec4(col, 1.0);
}
