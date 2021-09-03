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

#define time pc.time
const float PI = 3.141592;

mat2 rot(float a) {
    float c = cos(a), s = sin(a);
    return mat2(c, -s, s, c);
}

float map(vec3 p) {
    p.xy *= rot(time / 6.);
    p.zy *= rot(time / 3.);
    p.z += time / PI;
    p = fract(p) - 0.5;

    float d = length(p) - 0.2;

    return d;
}

vec3 norm(vec3 p) {
    mat3 k = mat3(p, p, p) - mat3(0.001);
    return normalize(map(p) - vec3(map(k[0]), map(k[1]), map(k[2])));
}

void main() {
    /* vec2 uv = (in_uv * pc.resolution - pc.resolution * 0.5) / pc.resolution.y; */
    vec2 uv = (in_uv - 0.5) * vec2(pc.resolution.x / pc.resolution.y, 1.);

    vec3 ro = vec3(0., 0., 0);
    vec3 rd = vec3(uv, 1. / (2. + sin(time)));
    /* rd = vec3(uv, 1.); */

    float d = 0;
    vec3 att = vec3(0.);
    for (int i = 0; i < 400; ++i) {
        vec3 p = rd * d + ro;
        float e = map(p);
        d += max(e * 0.5, abs(d - sin(time) * 5. - 6.) / 3e3 * sqrt(d));
        /* d += e; */
        att += (norm(p) * 0.5 + 0.5).y * exp(-e) / 8e2;
    }
    att -= d / 3e2;

    vec3 col = vec3(0.);
    col = att;

    col = pow(col, vec3(1. / 0.4));

    out_color = vec4(col, 1.);
}
