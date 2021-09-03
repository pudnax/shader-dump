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
#define TAU 2*PI
#define EPS vec3(0., 0.001, 0.0001)

void look_at(inout vec3 rd, vec3 ro, vec3 ta, vec3 up) {
    vec3 w = normalize(ta - ro), u = normalize(cross(w, up));
    rd = rd.x * u + rd.y * cross(u, w) + rd.z * w;
}

void rot(inout vec3 p, vec3 a, float t) {
    a = normalize(a);
    vec3 u = cross(a, p), v = cross(a, u);
    p = u * sin(t) + v * cos(t) + a * dot(a, p);
}

#define sabs(p, k) sqrt(p* p + k)

void sfold(inout vec3 p, float n, float k) {
    float c = cos(PI / n), s = sqrt(.75 - c * c);
    vec3 v = vec3(-0.5, -c, s);
    for (int i = 0; i < 5; ++i) {
        p.xy = sabs(p.xy, k);
        float g = dot(p, v);
        p -= (g - sabs(g, k)) * v;
    }
}

void sfold(inout vec2 p, float n, float k) {
    float h = floor(log2(n)), a = TAU * exp2(h) / n;
    for (float i = 0; i < h + 2; ++i) {
        vec2 v = vec2(-cos(a), sin(a));
        float g = dot(p, v);
        p -= (g - sabs(g, k)) * v;
        a *= 0.5;
    }
}

void sfold(inout vec2 p, vec2 v, float k) {
    float g = dot(p, v);
    p -= (g - sabs(g, k)) * v;
}

void sfold90(inout vec2 p, float k) {
    vec2 v = normalize(vec2(1., -1.));
    float g = dot(p, v);
    p -= (g - sabs(g, k)) * v;
}

float de0(vec3 p) {
    float k = 1e-3;
    p.z = sabs(p.z, k);
    sfold(p.xy, 5., k);
    vec3 v = normalize(vec3(2, 1, 3));
    return dot(p, v) - 0.6;
}

float de1(vec3 p) {
    float k = 5e-3;
    p = sabs(p, k);
    sfold90(p.xz, k);
    sfold90(p.yz, k);
    vec3 v = normalize(vec3(1, 1, -1));
    return dot(p, v) - 0.7;
}
float de2(vec3 p) {
    float k = 5e-4;
	float n = 4;
    sfold(p, n, k);
    vec3 v = normalize(vec3(1));
    return dot(p, v) - 1.;
}
float de3(vec3 p) {
    float k = 5e-4;
	float n = 5;
    sfold(p, n, k);
    vec3 v = normalize(vec3(1));
    return dot(p, v) - 1.;
}
float de4(vec3 p) {
    float k = 5e-4;
	float n = 5;
    sfold(p, n, k);
    vec3 v = normalize(vec3(0, 1, 1));
    return dot(p, v) - 1.;
}
float de5(vec3 p) {
    float k = 3e-3;
    p = sabs(p, k);
    sfold90(p.xz, k);
    sfold90(p.yz, k);
    sfold90(p.xy, k);
    vec3 v = normalize(vec3(1, 0, 1));
    return dot(p, v) - 1.;
}
float de6(vec3 p) {
    float k = 1e-3;
    p = sabs(p, k);
    sfold90(p.xz, k);
    sfold90(p.yz, k);
    sfold90(p.xy, k);
    vec3 v = normalize(vec3(2, 3, 1));
    return dot(p, v) - .9;
}
float de7(vec3 p) {
    float k = 3e-3;
    p = sabs(p, k);
    sfold90(p.xz, k);
    sfold90(p.yz, k);
    sfold90(p.xy, k);
    vec3 v = normalize(vec3(1, 2, -1));
    return dot(p, v) - .9;
}

float map(vec3 p) {
	rot(p, vec3(cos(pc.time * 0.3), sin(pc.time * 0.3), 1), pc.time*0.5);
	switch (int(mod(pc.time, 8.))) {
		case 0: return de0(p);
		case 1: return de1(p);
		case 2: return de2(p);
		case 3: return de3(p);
		case 4: return de4(p);
		case 5: return de5(p);
		case 6: return de6(p);
		case 7: return de7(p);
	}
	return 1.0;
}

vec3 normal(vec3 p) {
    mat3 k = mat3(p, p, p) - mat3(EPS.y);
    return normalize(map(p) - vec3(map(k[0]), map(k[1]), map(k[2])));
}

float ray_march(vec3 ro, vec3 rd, float near, float far) {
    float t = near, d = EPS.y;
    for (int i = 0; i < 100; ++i) {
        t += d = map(ro + rd * t);
        if (d < 0.001) return t;
        if (t >= far) return far;
    }
    return far;
}

vec3 do_color(vec3 p) {
	return vec3(0.3, 0.5, 0.8) + cos(p) * 0.5 + 0.5;
}

void main() {
    vec2 uv = (in_uv + -0.5) * 2.0 / vec2(pc.resolution.y / pc.resolution.x, 1);

    vec3 ro = vec3(0, 0, 5);
    vec3 rd = normalize(vec3(uv, 2));
    vec3 ta = vec3(0);
    look_at(rd, ro, ta, vec3(0, 1, 0));
    vec3 col = vec3(0);
    const float maxd = 50.;
    float t = ray_march(ro, rd, 0., maxd);

    if (t < maxd) {
        vec3 p = ro + rd * t;
        col = do_color(p);
        vec3 n = normal(p);
        vec3 light_pos = ro + vec3(2, 5, 2);
        vec3 li = light_pos - p;
        float len = length(li);
        li /= len;
        float dif = clamp(dot(n, li), 0., 1.);

        col *= max(dif, 0.2);

        float rimd = pow(clamp(1. - dot(reflect(-li, n), -rd), 0., 1.), 2.5);
        float frn = rimd + 2.2 * (1. - rimd);

        col *= frn * 0.6;
        col *= max(0.5 * 0.5 * n.y, 0.3);
        col *= exp2(-2. * pow(max(0., 2. - map(p + n * 0.8) / 0.8), 2.0));
        col += vec3(0.8, 0.6, 0.2) *
               pow(clamp(dot(reflect(rd, n), li), 0., 1.), 10.);
    }

	col = pow(col, vec3(0.7544));
	/* col = col / (1 + col); */

    out_color = vec4(col, 1.0);
}
