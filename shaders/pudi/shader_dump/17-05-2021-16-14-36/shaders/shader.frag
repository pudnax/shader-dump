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

const vec3 EPS = vec3(0., 0.01, 0.0001);
const uint MAX_STEPS = 256;
const float MAX_DIST = 100.;
const float PI = 3.141592;

float hash(vec3 uv) {
    float f = fract(sin(dot(uv, vec3(.09123898, .0231233, .0532234))) * 1e5);
    return f;
}

// 3d noise function (linear interpolation between hash of integer bounds)
float noise(vec3 uv) {
    vec3 fuv = floor(uv);
    vec4 cell0 = vec4(
        hash(fuv + vec3(0, 0, 0)),
        hash(fuv + vec3(0, 1, 0)),
        hash(fuv + vec3(1, 0, 0)),
        hash(fuv + vec3(1, 1, 0))
    );
    vec2 axis0 = mix(cell0.xz, cell0.yw, fract(uv.y));
    float val0 = mix(axis0.x, axis0.y, fract(uv.x));
    vec4 cell1 = vec4(
        hash(fuv + vec3(0, 0, 1)),
        hash(fuv + vec3(0, 1, 1)),
        hash(fuv + vec3(1, 0, 1)),
        hash(fuv + vec3(1, 1, 1))
    );
    vec2 axis1 = mix(cell1.xz, cell1.yw, fract(uv.y));
    float val1 = mix(axis1.x, axis1.y, fract(uv.x));
    return mix(val0, val1, fract(uv.z));
}

// fractional brownian motion
float fbm(vec3 uv) {
    float f = 0.;
    float r = 1.;
    for (int i = 0; i < 4; ++i) {
        f += noise((uv + 10.) * r) / (r *= 2.);
    }
    return f / (1. - 1. / r);
}

float easeOutBack(float x) {
	x = 1 - abs(x);

    const float c1 = 1.710634661 ;
    const float c3 = c1 + 1;

	float res = 1 + c3 * pow(x - 0.98, 3) + c1 * pow(x - 1, 2);
    /* res = c3 * x * x * x - c1 * x * x; */

    res = max(0., res);
	float crop = 1 - step(1, x);
	res *= crop;

    return res;
}

float opExtrussion( in vec3 p, in float sdf, in float h ) {
    vec2 w = vec2( sdf, abs(p.z) - h );
  	return min(max(w.x,w.y),0.0) + length(max(w,0.0));
}

// http://research.microsoft.com/en-us/um/people/hoppe/ravg.pdf
// { dist, t, y (above the plane of the curve, x (away from curve in the plane of the curve))
float det( vec2 a, vec2 b ) { return a.x*b.y-b.x*a.y; }
vec4 sdBezier(vec3 p, vec3 va, vec3 vb, vec3 vc) {
    vec3 w = normalize(cross(vc - vb, va - vb));
    vec3 u = normalize(vc - vb);
    vec3 v = cross(w, u);

    vec2 m = vec2(dot(va - vb, u), dot(va - vb, v));
    vec2 n = vec2(dot(vc - vb, u), dot(vc - vb, v));
    vec3 q = vec3(dot(p - vb, u), dot(p - vb, v), dot(p - vb, w));

    float mn = det(m, n);
    float mq = det(m, q.xy);
    float nq = det(n, q.xy);

    vec2 g = (nq + mq + mn) * n + (nq + mq - mn) * m;
    float f = (nq - mq + mn) * (nq - mq + mn) + 4.0 * mq * nq;
    vec2 z = 0.5 * f * vec2(-g.y, g.x) / dot(g, g);
    // float t = clamp(0.5+0.5*(det(z,m+n)+mq+nq)/mn, 0.0 ,1.0 );
    float t = clamp(0.5 + 0.5 * (det(z - q.xy, m + n)) / mn, 0.0, 1.0);
    vec2 cp = m * (1.0 - t) * (1.0 - t) + n * t * t - q.xy;

    float d2 = dot(cp, cp);
    return vec4(sqrt(d2 + q.z * q.z), t, q.z, -sign(f) * sqrt(d2));
}

vec3 DomainRotateSymmetry(const in vec3 vPos, const in float fSteps) {
    float angle = atan(vPos.x, vPos.z);

    float fScale = fSteps / (PI * 2.0);
    float steppedAngle = (floor(angle * fScale + 0.5)) / fScale;

    float s = sin(-steppedAngle);
    float c = cos(-steppedAngle);

    return vec3(c * vPos.x + s * vPos.z, vPos.y, -s * vPos.x + c * vPos.z);
}

#define sabs(x, k) sqrt((x) * (x) + k)

float smin(float a, float b, float k) {
    float h = max(k - abs(a - b), 0.0);
    return min(a, b) - h * h * 0.25 / k;
}

float sdTorus(vec3 p, vec2 t) {
    vec2 q = vec2(length(p.xz) - t.x, p.y);
    return length(q) - t.y;
}

float smax(float a, float b, float k) {
    k *= 1.4;
    float h = max(k - abs(a - b), 0.0);
    return max(a, b) + h * h * h / (6.0 * k * k);
}

vec3 erot(vec3 p, vec3 d, float ro) {
	return mix(dot(p, d) * d, p, cos(ro)) + cross(p, d) * sin(ro);
}

vec2 rot(vec2 p, float a) {
	float c = cos(a), s = sin(a);
	return mat2(c, -s, s, c) * p;
}

float GetDistanceGear(const in vec3 vPos) {
    vec3 vToothDomain = DomainRotateSymmetry(vPos, 16.0);
    vToothDomain.xz = abs(vToothDomain.xz);
    float fGearDist = dot(vToothDomain.xz, normalize(vec2(1.0, 0.55))) - 0.55;
    float fSlabDist = abs(vPos.y) - 1.;

    float fResult = max(fGearDist, fSlabDist);
    return fResult;
}

vec3 sdgEllipse(vec2 p, in vec2 ab) {
    vec2 sp = sign(p);
    p = abs(p);

    bool s = dot(p / ab, p / ab) > 1.0;
    float w = atan(p.y * ab.x, p.x * ab.y);
    if (!s)
        w = (ab.x * (p.x - ab.x) < ab.y * (p.y - ab.y)) ? 1.570796327 : 0.0;

    for (int i = 0; i < 4; i++) {
        vec2 cs = vec2(cos(w), sin(w));
        vec2 u = ab * vec2(cs.x, cs.y);
        vec2 v = ab * vec2(-cs.y, cs.x);
        w = w + dot(p - u, v) / (dot(p - u, u) + dot(v, v));
    }
    vec2 q = ab * vec2(cos(w), sin(w));

    float d = length(p - q);
    return vec3(d, sp * (p - q) / d) * (s ? 1.0 : -1.0);
}

float sdRoundedCylinder(vec3 p, float ra, float rb, float h) {
    vec2 d = vec2(length(p.xz) - 2.0 * ra + rb, abs(p.y) - h);
    return min(max(d.x, d.y), 0.0) + length(max(d, 0.0)) - rb;
}

float sdEllipsoid(vec3 p, vec3 r) {
    float k0 = length(p / r);
    float k1 = length(p / (r * r));
    return k0 * (k0 - 1.0) / k1;
}

float sdCone(vec3 p, vec2 c, float h) {
    float q = length(p.xz);
    return max(dot(c.xy, vec2(q, p.y)), -h - p.y);
}

float sdCappedCone(vec3 p, vec3 a, vec3 b, float ra, float rb) {
    float rba = rb - ra;
    float baba = dot(b - a, b - a);
    float papa = dot(p - a, p - a);
    float paba = dot(p - a, b - a) / baba;
    float x = sqrt(papa - paba * paba * baba);
    float cax = max(0.0, x - ((paba < 0.5) ? ra : rb));
    float cay = abs(paba - 0.5) - 0.5;
    float k = rba * rba + baba;
    float f = clamp((rba * (x - ra) + paba * baba) / k, 0.0, 1.0);
    float cbx = x - ra - f * rba;
    float cby = paba - f;
    float s = (cbx < 0.0 && cay < 0.0) ? -1.0 : 1.0;
    return s * sqrt(min(cax * cax + cay * cay * baba,
                        cbx * cbx + cby * cby * baba));
}

float sdCappedCylinder(vec3 p, vec3 a, vec3 b, float r) {
    vec3 ba = b - a;
    vec3 pa = p - a;
    float baba = dot(ba, ba);
    float paba = dot(pa, ba);
    float x = length(pa * baba - ba * paba) - r * baba;
    float y = abs(paba - baba * 0.5) - baba * 0.5;
    float x2 = x * x;
    float y2 = y * y * baba;
    float d = (max(x, y) < 0.0)
                  ? -min(x2, y2)
                  : (((x > 0.0) ? x2 : 0.0) + ((y > 0.0) ? y2 : 0.0));
    return sign(d) * sqrt(abs(d)) / baba;
}

float sdTriangleWave(in vec3 p, in float f, in float a, in float d, in float t) {
    float pw = 1.0 / f, qw = 0.25 * pw;
    vec2 sc = vec2(2.0 * a, pw);
    float l = length(sc);
    p.x = abs(mod(p.x + qw, pw) - 0.5 * pw) - qw;
    p.xy *= mat2(sc, -sc.y, sc.x) / l;
    return length(vec3(p.x, max(abs(p.yz) - vec2(0.25 * l, d), 0.0))) - t;
}

float logPolarPolka(vec2 pos) {
    // Apply the forward log-polar map
    pos = vec2(log(length(pos)), atan(pos.y, pos.x));

    // Scale everything so tiles will fit nicely in the ]-pi,pi] interval
    pos *= 6.0 / PI;

    // Convert pos to single-tile coordinates
    pos = fract(pos) - 0.5;

    // Return color depending on whether we are inside or outside a disk
    return 1.0 - smoothstep(0.3, 0.31, length(pos));
}

float topping2(in vec3 p) {
	vec3 a = vec3(0., 0.0, -0.0);
	vec3 bb = vec3(-0.1, 1.2, 0.2);
	vec3 c = vec3(-0.2, 1.5, -0.0);

	vec4 b = sdBezier(p, a, bb, c);

	vec2 q = rot(b.zw, PI);

	vec2 id2 = round(q / 0.1);
	id2 = clamp(id2, vec2(0), vec2(2, 1));
	/* q -= 0.1 * id2; */

	float id = 11. * id2.x + id2.y * 13.;

	/* q += smoothstep(0.5, 0.8, b.y) * 002 * vec2(0.4, 1.5) * cos(23.0 * b.y + id*vec2(13, 17)); */
	vec2 ooc_id;
	ooc_id.x = clamp(length(q) * 8.0 - 0.2, 0.0, 1.0);
	vec4 res = vec4(99, q, b.y);
	for (int i = 0; i < 3; ++i) {
		vec2 tmp = q + 0.01 * cos(id + 180 * b.y + vec2(2*i, 6 - 2 * i));
		float lt = length(tmp) - 0.02;
		if ( lt < res.x) {
			ooc_id.y = id + float(i);
			res.x = lt;
			res.yz = tmp;
		}
	}

	/* vec4 bbb = sdBezier(p, a, bb, c); */
	/* res.x = bbb.x - 0.1; */

	return res.x;
}

float topping3(vec3 p) {
	vec3 pp = p;
	pp.y -= easeOutBack(pp.y);
    float sphere = length(pp + vec3(0.0, 1.1399993, 0.0)) - 0.5;

    float ring = sdTorus(p, vec2(0.45, 0.1));

    vec3 q = p * 100;
	q = erot(q, normalize(vec3(-0.11, 0.2, 0.0)), PI);
    vec3 q2 = p * 100;
	q2 = erot(q2, normalize(pc.pos), PI);
    /* float bumps = cos(q2.y) - sin(q2.y); // * cos(q.z) * sin(q.x + PI * 3 / 4); */
    float bumps = cos(q.y) - sin(q.y); // * cos(q.z) * sin(q.x + PI * 3 / 4);
	bumps = bumps * 0.5 + 0.5;
	bumps *= 0.01;

	vec3 pc = p + vec3(0.0, -0.54, 0.0);
	float cone = sdCone(pc, vec2(0.5, 0.57), 0.1);

    float res;
	res = smin(sphere, cone, 0.01);
    res = res + bumps;
	res = max(res, -p.y);
	res = smin(res, ring, 0.1);

    return res;
}

float sdRect(vec2 p, vec2 r) {
    p = abs(p) - r;
	return min(max(p.x, p.y), 0.) + length(max(p, 0.));
}

float opSU(float a, float b, float k) {
    float h = clamp(.5 + .5 * (b - a) / k, 0., 1.);
    return mix(b, a, h) - k * h * (1. - h);
}

float sdSphere(vec3 p, float r) {
	return length(p) - r;
}

float topping(vec3 p) {
    float ring = sdTorus(p + vec3(0., 0.34, 0.0), vec2(0.38, 0.1));

	vec3 q = erot(p, normalize(vec3(-0.2, 1., 0.1)),  3.);
	float arg = q.x * q.z * 100 + pc.time;
    float bumps = cos(q.y * 100); + sin(arg)*cos(arg);
	bumps = bumps * 0.5 + 0.5;
	bumps *= 0.01;
	ring += bumps;

	p = erot(p, vec3(0., 1., 0.), p.y * 3.);
	float rect = sdRect(p.xz, vec2(0.5));

	p = erot(p, vec3(0., 1., 0.), PI / 4.);
	float neigh = sdRect(p.xz, vec2(0.4));

	float res = opSU(rect, neigh, 0.1);
	res += p.y + 0.2;
	res = -opSU(-res * 0.5, -sdSphere(p, 0.5), 0.1);
	res = smin(ring, res, 0.01);
	res += fbm(p*10) * 0.03;

    return res;
}

float scene(vec3 p) {
    float t = pc.time;
    float wiggle = cos((p.y + t) * 11) * 0.02;
    float wigglee = cos((p.y + t + 0.05) * 11) * 0.02;

    float res = length(p) + 0.5;
    {
        vec3 pa = vec3(0.0, -0.4, 0.0);
        float body = sdCappedCone(p, pa, vec3(0.), 0.3, 0.2) - 0.1;
        res = body +wiggle;
    }

    {
        vec3 pp = p + vec3(0, .5, 0);
        pp.yz = rot(pp.yz, PI / 2);
        float r_plate = 0.5;
        float plate = sdEllipsoid(pp, vec3(r_plate, r_plate, 0.07));
        res = min(res, plate);
    }

    {
        vec3 pe = vec3(abs(p.x), p.yz);
        pe += vec3(-0.12, 0.05, 0.3);
        float r_eye = 0.08;
        float eyes = sdEllipsoid(pe, vec3(r_eye, r_eye, 0.05));
        res = smin(res, eyes, 0.02);
    }

    {
        vec3 pb = vec3(abs(p.x), p.yz);
        pb += vec3(-0.17, 0.19, 0.36);
        vec2 r_blush = vec2(0.04, 0.02);
        float blush = sdEllipsoid(pb, vec3(r_blush, 0.02));
        blush = max(max(blush, pb.z), -pb.z + 0.001);
        res = min(res, blush);
    }

    {
        vec3 po = p - vec3(0., 0.2, 0.);
        vec4 ahoge =
            sdBezier(po, vec3(0., -0.12 + wiggle, 0.), vec3(-0.03, 0.08, -0.01),
                     vec3(-0.08, 0.01 + wigglee * 0.4, 0.));
        ahoge.x -= 0.027 * ahoge.y + (1 - ahoge.y) * 0.008;
        res = min(res, ahoge.x + 0.003);
    }

    {
		p.y += wiggle;
        vec3 pt = vec3(abs(p.x), p.yz) + vec3(-0.17, -0.15899997, 0.0);
        float scale = .1;
        float top = topping(pt / scale) * scale;
        res = smin(res, top, 0.03);
    }

	{
		vec3 pm = p + vec3(0., 0.38, 0.22);
		float exterrior = sdEllipsoid(pm, vec3(0.19000003, 0.24000011, 0.24000005));
		exterrior = max(exterrior, -(p.y + 0.424));
		vec4 umuer =
            sdBezier(pm + vec3(-0.02, -0.18, 0.089999974), vec3(0.01, -0.12, 0.55), vec3(0.06, 0.19, 0.),
                     vec3(0.01, 0.01, -0.04));
		umuer.x -= 0.03 + (1 - umuer.y)*0.04;
		exterrior = smax(exterrior, -umuer.x, 0.08);
		res = max(res, -exterrior);
	}

    return res;
}

vec3 norm(vec3 p) {
    mat3 k = mat3(p, p, p) - mat3(0.001);
    return normalize(scene(p) - vec3(scene(k[0]), scene(k[1]), scene(k[2])));
}

vec3 camera(vec2 uv, vec3 ro, vec3 at, vec3 up, float z) {
    vec3 f = normalize(at - ro),
		r = cross(up, f),
		u = cross(f, r),
        c = ro + f * z,
		i = c + uv.x * r + uv.y * u;
    return i - ro;
}

float height = 0.01;
vec3 color(in vec3 p) {
	p.y += 0.43;
    vec3 top = vec3(0.9921875, 0.984375, 0.25578125);
    vec3 ring = vec3(0.6, 0.04, 0.0);
    vec3 bottom = vec3(0.3, 0.3, 0.3);
    bottom = mix(vec3(0.0), bottom, min(1.0, -1.0 / (p.y * 20.0 - 0.0)));
    vec3 side = mix(bottom, ring, smoothstep(-height - 0.001, -height, p.y));
    return mix(side, top, smoothstep(-0.01, 0.0, p.y));
}

void main() {
    vec2 uv = (in_uv + -0.5) * 2.0 / vec2(pc.resolution.y / pc.resolution.x, 1);

	float t = pc.time * 0.5;

	vec3 ro = vec3(-0., 0.5, -1.9);
    vec3 at = vec3(0., 0.0, -0.0);
    vec3 rd = camera(uv, ro, at, vec3(0, 1, 0), 1.);

	ro += pc.pos;
	/* ro += vec3(0.01, -0.5099998, 1.1099993); */

    ro.xz = rot(ro.xz, sin(t)*0.3);
	rd.xz = rot(rd.xz, sin(t)*0.3);

	float dist = EPS.y;
	bool hit = false;
	for (int i = 0; i < MAX_STEPS; ++i) {
		if (dist >= MAX_DIST) break;
		vec3 pos = ro + rd * dist;
		float res = scene(pos);
		if (abs(res) < 0.001) { hit = true; break; }
		dist += res * 0.9;
	}
	vec3 pos = ro + rd * dist;

	vec3 col = vec3(0.1);
    if (hit) {
        vec3 nor = norm(pos);
        float dif = clamp(dot(nor, vec3(0.57703)), 0.0, 1.0);
        float amb = 0.5 + 0.5 * dot(nor, vec3(0.0, 1.0, 0.0));
        col = color(pos) * amb + color(pos) * dif;

		col = nor * 0.5 + 0.5;
    }

    out_color = vec4(col, 1.0);
}
