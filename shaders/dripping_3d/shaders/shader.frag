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
const float PI = 3.14159265359;

float seed = 0.25;

float dripDistance = 0.1;
float density = 0.75;
float bCurve = 1.5;
float bFreq = 3.5;
float bRange = 0.35;
float fallSpeed = 6.0;
float sdfWidth = 0.18;

const vec3 EPS = vec3(0., 0.01, 0.0001);
const float MAX_DIST = 200.;
const uint MAX_STEP = 256;

float rand(vec2 co){ return fract(sin(dot(co.xy ,vec2(12.9898,78.233))) * 43.5453); }
float rand(float x, float y) { return rand(vec2(x, y)); }
float rand(vec3 v) { return fract(sin(dot(v + vec3(-8.5123, 23.2156, 0.0), vec3(12.9898, 6.233, 0.84261))) * 47583.5453123); }
float rand(vec2 v, float x) { return rand(vec3(v, x)); }
float rand(float x, vec2 v) { return rand(vec3(x, v)); }

vec3 erot(vec3 p, vec3 ax, float angle) {
	return mix(dot(ax, p) * ax, p, cos(angle)) + cross(ax, p) * sin(angle);
}
#define sat(x) clamp((x), 0.0, 1.0);

float dripSDF(vec2 uv) {
    float time = pc.time;

    float s = sdfWidth * abs((1.0 - uv.y) - 0.75) + 0.05;
    float o = 1.0;
    float drip = 999.0;

    float x = uv.x - sdfWidth;
    x += dripDistance - mod(x, dripDistance);

    x -= dripDistance;  // ungh... this is dirty... I'll fix it later
    for (int i = 0; i < 1000; i++) {
        if (x > uv.x + sdfWidth)
            break;

        x += dripDistance;
        float isLine = round(rand(x, seed) + density - 0.5);
        if (isLine == 0.0)
            continue;

        float y = rand(seed, x) * 0.8 + 0.1;
        // y *= abs(sin(x*3.0))*0.5 + 0.5;
        float animTime = time + (y * 10.0);
        float bounce = 0.0 - (bCurve * mod(animTime, bFreq)) *
                                 exp(1.0 - bCurve * mod(animTime, bFreq));
        y += bounce * bRange;
        y = min(y, uv.y);

        float f = y + mod(animTime, bFreq) * fallSpeed * bRange;

        // float d = min( distance(vec2(x,y),uv), distance(vec2(x,f),uv) );
        float d = distance(vec2(x, y), uv);

        o *= clamp(d / s, 0.0, 1.0);
        drip = min(drip, distance(vec2(x, f), uv));
    }

    o = min(o, clamp(drip / s, 0.0, 1.0));

    //  s = sin(uv.x*20.0+ time * 0.2)*0.3 + 0.4;
    return o * clamp(distance(0.0, uv.y) / s, 0.0, 1.0);
}

float drip3d(vec3 p) {
	float time = pc.time;

	float s = sdfWidth * abs((1.0 - p.y) - 0.75) + 0.05;
	float o = 1.0;
	float  drip = 999.0;

	vec2 x = p.xz - sdfWidth;
	x += dripDistance - mod(x, dripDistance);

	x -= dripDistance;
	for (int i = 0; i < 2; ++i) {
		if ((x.x > p.x + sdfWidth) || (x.y > p.z + sdfWidth)) break;

		x += dripDistance;
		float isLine = round(rand(x, seed) + density - 0.5);
		if (isLine == 0.0) continue;

		float y = rand(seed, x) * 0.8 + 0.1;
		/* y *= abs(sin(x * 3.0)) * 0.5 + 0.5; */
		float animTime = time + (y * 10.0);
		float bounce = 0.0 - (bCurve * mod(animTime, bFreq)) * exp(1.0 - bCurve * mod(animTime, bFreq));
		y += bounce * bRange;
		y = min(y, p.y);

		float f = y + mod(animTime, bFreq) * fallSpeed * bRange;

		float d = distance(vec3(x,y), p);

		o *= sat(d / s);
		drip = min(drip, distance(vec3(x, f), p));
	}

	o = min(o, clamp(drip / s, 0.0, 1.0));

	s = sin(p.x * 20.0 + time * 0.2) * 0.3 + 0.4;
	return o * clamp(distance(0.0, p.y) / s, 0.0, 1.0);
}

float scene(vec3 p) {
    float res;

    float sphere =
        length(p) - 0.5 +
        sin(p.y * 10 + pc.time) * sin(p.x + pc.time) * 0.1;

    float drip = drip3d(p);

    res = min(res, drip);

	res = drip;


    return res;
}

vec3 norm(vec3 p) {
    mat3 k = mat3(p, p, p) + mat3(EPS.y);
    return vec3(p - vec3(scene(k[0]), scene(k[1]), scene(k[1])));
}

vec3 camera(vec2 uv, vec3 ro, vec3 at, vec3 up, float z) {
    vec3 f = normalize(at - ro),
		r = cross(up, f),
		u = cross(f, r),
        c = ro + f * z,
		i = c + uv.x * r + uv.y * u;
    return i - ro;
}

void main() {
    vec2 uv = (in_uv + -0.5) * 2.0 / vec2(pc.resolution.y / pc.resolution.x, 1);

    /* vec3 ro = vec3(0.0, 0.0, 2); */
    /* vec3 rd = normalize(vec3(uv, -1)); */
	vec3 ro = vec3(-0., 0.2, -1.0);
    vec3 at = vec3(0., -0.0, -0.0);
    vec3 rd = camera(uv, ro, at, vec3(0, 1, 0), 1.);

	ro += pc.pos;

    float dist = EPS.y;
    bool hit = false;
    for (int i = 0; i < MAX_DIST; ++i) {
        if (dist >= MAX_DIST) break;
        float res = scene(ro + rd * dist);
        if (abs(res) < EPS.z) { hit = true; break; }
        dist += res;
    }
    vec3 pos = ro + rd * dist;

    vec3 col = vec3(0.2);
    if (hit) {
        col = norm(pos) * 0.5 + 0.5;
		col = vec3(fract(pos));
    }

	col = pow(col, vec3(2.4));
    out_color = vec4(col, 1.0);
}
