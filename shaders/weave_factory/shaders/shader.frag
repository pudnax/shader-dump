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
#define TAU (2. * PI)

#define od(p, d) (dot(p, normalize(sign(p))) - d)

#define rot(a) mat2(cos(a), sin(a), -sin(a), cos(a))
#define crep(p, c, l) p = p - c * clamp(round(p / c), -l, l)
#define rep(p, c) p = (mod(p, c) - c * .5)

#define frt(sp, off) fract((pc.time + off) * sp)
#define flt(sp, off) floor((pc.time + off) * sp)

struct obj {
    float d;
    vec3 shadowcol;
    vec3 lightcol;
};

obj minobj(obj a, obj b) {
    if (a.d < b.d)
        return a;
    else
        return b;
}

float box(vec3 p, vec3 c) {
    vec3 q = abs(p) - c;
    return min(0., max(q.x, max(q.y, q.z))) + length(max(q, 0.));
}

obj set(vec3 p, vec3 scol, vec3 lcol) {
    float id = round(p.x / 5.);
    crep(p.x, 5., 2.);
    p.y += sin(p.z * 0.5) * 0.5;
    vec3 pp = p;

    rep(p.z, 5.);
    float speed = 2., offset = id * 0.2;
    p.yz *= rot(PI / 2. * (flt(speed, offset) + pow(frt(speed, offset), 1.5)));

    float d = mix(box(p, vec3(.5)), od(p, 0.6), 0.5);

    p = pp;
    p.z += pc.time;
    p.y += 0.65;
    rep(p.z, 1.1);
    d = min(d, box(p, vec3(1.5, 0.1, 0.5)));

    return obj(d, scol, lcol);
}

obj SDF(vec3 p) {
    p.yz *= rot(-atan(1. / sqrt(2)));
    p.xz *= rot(PI / 4.);

    obj d = set(p, vec3(0.1, 0.4, 0.2), vec3(0.95, 0.9, 0.3));

    p.xy += vec2(-2., 2.5);
    p.xz *= rot(PI / 2.);
    d = minobj(d, set(p, vec3(0.4, 0., 0.1), vec3(0.5, 0.8, 0.9)));

    return d;
}

vec3 getnorm(vec3 p) {
    vec2 eps = vec2(0.001, 0.);
    return normalize(SDF(p).d - vec3(SDF(p - eps.xyy).d, SDF(p - eps.yxy).d,
                                     SDF(p - eps.yyx).d));
}

float AO(float eps, vec3 n, vec3 p) {
    return SDF(p + eps * n).d / eps;
}

void main() {
    vec2 uv = (in_uv + -0.5) * 2.0 / vec2(pc.resolution.y / pc.resolution.x, 1);

    vec3 ro = vec3(uv * 5., -30.);
    vec3 rd = vec3(0., 0., 1.);
    vec3 p = ro;
    vec3 col = vec3(0.);
    vec3 l = normalize(vec3(2., 3., -2.));

    bool hit = false;
    obj O;

    for (float i = 0.; i < 100.; ++i) {
        O = SDF(p);
        if (O.d < 0.001) {
            hit = true;
            break;
        }
        p += O.d * rd * 0.8;
    }

    if (hit) {
        vec3 n = getnorm(p);
        float li = max(dot(n, l), 0.);
        float ao = AO(0.1, n, p) + AO(0.4, n, p) + AO(0.3, n, p);
        col = mix(O.shadowcol, O.lightcol, li) * ao / 2.5;
    }

    out_color = vec4(col, 1.0);
}
