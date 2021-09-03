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

#define Color vec3
#define Point vec2
#define lerp mix

float hash(float n) {
    return fract(sin(n) * 1e4);
}

float noise(float x) {
    float i = floor(x), f = fract(x);
    float u =
        // 0.5;
        // f;
        f * f * (2.0 - 2.0 * f);
    return 2.0 * mix(hash(i), hash(i + 1.0), u) - 1.0;
}

void getSkyColor(float x, float y, inout Color color) {
    float h = max(0.0, 1.2 - y - pow(abs(x - 0.5), 3.0));
    color.r = pow(h, 3.0);
    color.g = pow(h, 7.0);
    color.b = 0.2 + pow(max(0.0, h - 0.1), 10.0);
}

float terrain(float x) {
    float y = 0.0;

    for (int octave = 0; octave < 10; ++octave) {
        float k = float(1 << octave);
        y += noise(x * k) / k;
    }

    return y * 0.3 + 0.36;
}

float water(float x) {
    return (sin(x * 71.0 - pc.time * 0.7) * 0.5 +
            sin(200. * x - pc.time * 8.0) * 0.002 + 0.25);
}

float tree(float x, float y) {
    if (y < 0.5) {
        return 0.0;
    } else {
        return 0.2 *
               max(0.0, max(max(abs(sin(x * 109.0)), abs(sin(x * 150.0))),
                            abs(sin(x * 117.0))) +
                            noise(37.0 * x) + noise(64.0 * x + 100.0) - 1.6);
    }
}

void main() {
    vec2 uv = (in_uv + -0.5) * 2.0 / vec2(pc.resolution.y / pc.resolution.x, 1);
    uv = (in_uv) / vec2(pc.resolution.y / pc.resolution.x, 1);
    float x = uv.x, y = uv.y;

    vec3 color = vec3(0.);

    getSkyColor(x, y, color);

    float shift = 0.09 * pc.time + 0.2;
    x += shift;

    float h = max(water(x), terrain(x));
    h += tree(x, h);

    if (y < h) {
        color = vec3(0.);
    }

    out_color = vec4(color, 1.);
}
