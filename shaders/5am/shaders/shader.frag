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

float stg;
vec3 erot(vec3 p, vec3 ax, float ro) {
    return mix(dot(ax, p) * ax, p, cos(ro)) + sin(ro) * cross(ax, p);
}

float smin(float a, float b, float k) {
    float h = max(k - abs(a - b), 0.) / k;
    return min(a, b) - h * h * h * k / 6.;
}

float comp(vec3 p) {
    p = asin(sin(p) * 0.999);
    return dot(p, normalize(vec3(1, 2, 3)));
}

float linedist(vec3 p, vec3 a, vec3 b) {
    float k = dot(p - a, b - a) / dot(b - a, b - a);
    return distance(p, mix(a, b, clamp(k, 0., 1.)));
}

float w,dts;
float scene(vec3 p) {
    float h1 = comp(erot(p, normalize(vec3(1, 2, 3)), 2.));
    float h2 = comp(erot(p, normalize(vec3(3, -1, 2)), stg));
    float h3 = comp(erot(p, normalize(vec3(0, 3, 2)), stg * 2.));
    float cave = (h1 + h2 + h3) / 2.5;
    float lvl1 = -smin(-cave, -p.z, 0.1);
    float lvl2 = -smin(-cave - 0.5, 1. - p.z, 0.1);

    vec3 p2 = p;
    p2.xy = asin(sin(p2.xy * 3.)) / 3.;
    p2.xy = abs(p2.xy);
    if (p2.x < p2.y) p2.xy = p2.yx;
    p2.z += asin(sin(pc.time + stg * 2.) * 0.9) * 0.5 + 2.;
    w = smin(lvl2, lvl1, 0.1);
    dts = linedist(p2, vec3(1, 0, 0), vec3(-1, 0, 0)) - 0.02;
    return min(w, dts);
}

vec3 norm(vec3 p) {
	mat3 k = mat3(p,p,p) - mat3(0.001);
	return normalize(scene(p) - vec3(scene(k[0]), scene(k[1]), scene(k[2])));
}

void main() {
    vec2 uv = (in_uv + -0.5) * 2.0 / vec2(pc.resolution.y / pc.resolution.x, 1);

    stg = floor(pc.time);

    vec3 cam = normalize(vec3(1.5, uv));
    vec3 init = vec3(-15, 0., 0.);
    float roz = pc.time * 0.1;
    if (mod(stg, 2.) < 1.) roz *= -1.;
    float roy = 0.8 + sin(stg) * 0.2;
    cam = erot(cam, vec3(0, 1, 0), roy);
    init = erot(init, vec3(0, 1, 0), roy);
    cam = erot(cam, vec3(0, 0, 1), roz);
    init = erot(init, vec3(0, 0, 1), roz);

    vec3 p = init;
    bool hit = false;
    float dist;
    float glo = 0.;
    for (int i = 0; i < 150 && !hit; ++i) {
        dist = scene(p);
        hit = dist * dist < 1e-6;
        glo += dist / (0.1 + dts * 500.) * 10.;
        p += cam * dist;
    }
    bool isgrid = (dts == dist);
    vec3 n = norm(p);
    vec3 r = reflect(cam, n);
    float fres = 1. - abs(dot(cam, n)) * 0.98;
    r.xy = abs(r.xy);
    n.xy = abs(n.xy);

    float ao = smoothstep(-0.5, 0.5, scene(p + n * 0.5));
    float ro = smoothstep(-0.2, 0.2, scene(p + r * 0.2));
    float spec = length(sin(r * 3.5) * 0.5 + 0.5) / sqrt(3.);
    float diff = length(sin(n * 3.5) * 0.5 + 0.5) / sqrt(3.);
    spec = fres * (spec * 0.1 + pow(spec, 8.) * 10.) * ro;
    vec3 col1 = vec3(1, 0.4, 0.1) * spec;
    vec3 col2 = vec3(diff) * vec3(0.01, 0.03, 0.1) * ao + spec * 0.05;
    float aaa = sin(p.x * 3. + p.y * 2. + p.z * 8.);
    vec3 col = mix(col2, col1, smoothstep(-0.94, -0.96, aaa));
    if (isgrid) col = vec3(0.1, 0.2, 0.5) * 6.;
    col = mix(col, vec3(0.01), smoothstep(0., -7., p.z));
    col += glo * vec3(0.1, 0.2, 0.5) + glo * glo;
    out_color = vec4(col, 1.0);
    out_color *= 1. - dot(uv, uv) * 0.3;
    out_color =
        smoothstep(-0.05, 1.1, sqrt(out_color) + vec4(0.08, 0.02, 0.1, 0));
    float dt = length(uv);
    float f1 = step(0., sin(stg * 0.7));
    float f2 = step(0., sin(stg * 1.7));
    float f3 = step(0., sin(stg * 2.7));
    float circles =
        abs(abs(abs(dt - mix(0.7, 0.5, f3)) - 0.1 * f1) - 0.03 * f2) - 0.002;
    out_color =
        mix(vec4(0.2), out_color, smoothstep(0., fwidth(dt) * 3., circles));
    out_color =
        mix(vec4(1), out_color, smoothstep(0., fwidth(dt) * 1.5, circles));
}
