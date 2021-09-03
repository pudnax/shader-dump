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

#define time pc.time

uvec4 s0, s1;
ivec2 pixel;

void rng_initialize(vec2 p, uint frame) {
    pixel = ivec2(p);

    //white noise seed
    s0 = uvec4(p, uint(frame), uint(p.x) + uint(p.y));

    //blue noise seed
    s1 = uvec4(frame, frame * 15843, frame * 31 + 4566, frame * 2345 + 58585);
}

// https://www.pcg-random.org/
uvec4 pcg4d(inout uvec4 v) {
    v = v * 1664525u + 1013904223u;
    v.x += v.y*v.w; v.y += v.z*v.x; v.z += v.x*v.y; v.w += v.y*v.z;
    v = v ^ (v >> 16u);
    v.x += v.y*v.w; v.y += v.z*v.x; v.z += v.x*v.y; v.w += v.y*v.z;
    return v;
}
float rand() { return float(pcg4d(s0).x) / float(0xffffffffu); }
vec2 rand2() { return vec2(pcg4d(s0).xy) / float(0xffffffffu); }

mat2 rot(float a) {
    float c = cos(a), s = sin(a);
    return mat2(c, -s, s, c);
}

float sphere(vec3 p, float r) {
    return length(p) - r;
}

float box(vec3 p, vec3 s) {
    p = abs(p) - s;
    return length(max(p, 0.0)) + min(max(p.x, max(p.y, p.z)), 0.0);
}

float caps(vec3 p, vec3 p1, vec3 p2, float s) {
    vec3 pa = p - p1;
    vec3 pb = p2 - p1;
    float proj = dot(pa, pb) / dot(pb, pb);
    proj = clamp(proj, 0., 1.);
    return length(p1 + pb * proj - p) - s;
}

int scene = 0;

float map(vec3 p) {
    vec3 p2 = p;
    float t = time * 0.1;
    p2.yz *= rot(t);
    p2.yx *= rot(t * 1.3);
    float d4 = max(box(p2, vec3(3.)), -sphere(p, 1.2));

    float d = box(p, vec3(0.4));

    if (scene == 0) d = d4;

    return d;
}

vec3 norm(vec3 p) {
    mat3 k = mat3(p, p, p) - mat3(0.001);
    return normalize(map(p) - vec3(map(k[0]), map(k[1]), map(k[2])));
}

float rnd(vec2 uv) {
    return fract(dot(sin(uv * 4532.714 + uv.yx *543.524), vec2(352.887)));
}

#define pcount 10
vec3 points[pcount];
int pid = 1;

float atm = 0.;

float map2(vec3 p) {
    float d = map(p);

    float d2 = 10000.;
    for (int i = 0; i < pid - 1; ++i) {
        float d3 = caps(p, points[i], points[i + 1], 0.01);
        atm += 0.013 / (0.05 + abs(d3)) * smoothstep(4., 0.3, d3);

        d2 = min(d2, d3);
    }
    return min(abs(d), d2);
}

vec3 process_hit(inout float dist, vec3 r, vec3 p, inout float side, float ior) {
    vec2 off = vec2(0.01, 0);
    vec3 n = side * normalize(dist - vec3(map(p - off.xyy), map(p - off.yxy), map(p - off.yyx)));
    vec3 r_new = refract(r, n, 1. - side * ( 0.3 + 0.1 * ior));
    if (dot(r_new, r_new) < 0.5) r_new = reflect(r, n);
    side = -side;
    dist = 0.1;
    return r_new;
}

void main() {
    vec2 uv = (in_uv + -0.5) * 2.0 / vec2(pc.resolution.y / pc.resolution.x, 1);

    /* scene = int(mod(floor(time / 10.), 2.)); */

    rng_initialize(in_uv * pc.resolution, pc.frame);

    vec3 s2 = vec3(18., 0., 0.);
    vec3 r2 = normalize(vec3(-1., sin(time) * 0.06, 0));

    float ior = rand() - 0.5;
    float id = ior * 2.;
    vec3 diff = 1.3 - vec3(1. + id, 0.45 + abs(id), 1. - id);

    vec3 p2 = s2;
    points[0] = p2;
    pid = 1;
    float side = 1.;
    for (int i = 0; i < 60; ++i) {
        float d = abs(map(p2));
        if (d < 0.001) {
            points[pid] = p2;
            pid += 1;
            if (pid >= pcount - 1) break;

            r2 = process_hit(d, r2, p2, side, ior);
        }
        if (d > 100.0) break;
        p2 += r2 * d;
    }
    points[pid] = p2 + r2 * 1000.;
    pid += 1;

    vec3 s = vec3(0., 0., -10.);
    vec3 r = normalize(vec3(uv, 1.));

    float rg = rand();
    float mumu = mix(rg, 1., 0.95);
    vec3 p = s;
    float side2 = 1.;
    for (int i = 0; i < 90; ++i) {
        float d = abs(map2(p));
        if (d < 0.001) {
            r = process_hit(d, r, p, side, ior);
        }
        if (d > 100.) break;
        p += r * d * mumu;
    }

    vec3 lazer = diff * atm;
    vec3 n = norm(p);
    vec3 col = vec3(0.);
    vec3 light_dir = normalize(vec3(1.));
    float shade = dot(light_dir, n);
    float amb = 0.2 * (mix(max(shade, 0.), shade * 0.5 + 0.5, .05));
    col = mix(vec3(amb), lazer, vec3(0.9));

    /* col = smoothstep(0.01, 0.9, col); */
    /* col = pow(col, vec3(0.4545)); */
    out_color = vec4(col, 1.0);
}
