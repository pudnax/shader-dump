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

#define EPSILON 0.0001
#define NEARCLIP 0.0001
#define FARCLIP 100.0
#define MARCHSTEPS 500

#define ALL 1
#define CUBEMAP_REFLECTION 1
#define BLINN_PHONG 1
#define AO 1
#define SHADOWS 1
#define SCENE_REFLECTION 1
#define SSS 1
#define VIGNETTE 1
#define COLOR 1

#define SSS_STEPS 25

struct Material {
    vec3 ambient;
    vec3 diffuse;
    vec3 specular;
    float shininess;
    float reflection;
    float sss;
};

struct DirLight {
    vec3 dir;
    vec3 color;
};

struct RayResult {
    float dist;
    Material material;
};

const Material kNoMaterial = Material(
    vec3(0.0, 0.0, 0.0),
    vec3(0.0, 0.0, 0.0),
    vec3(0.0, 0.0, 0.0),
    0.0,
    0.0,
    0.0
);

const Material kMaterialRed = Material(
    vec3(0.2, 0.0, 0.0),
    vec3(0.9, 0.0, 0.0),
    vec3(1.0, 0.9, 0.9),
    512.0,
    0.0,
    5.0
);

const Material kMaterialGreen = Material(
    vec3(0.0, 0.1, 0.0),
    vec3(0.0, 0.9, 0.0),
    vec3(0.9, 1.0, 0.9),
    0.0,
    0.09,
    4.0
);

const Material kMaterialBlue = Material(
    vec3(0.0, 0.0, 0.1),
    vec3(0.0, 0.0, 0.9),
    vec3(0.9, 9.0, 1.0),
    1024.0,
    0.3,
    9.0
);


DirLight kDirLight = DirLight(
    vec3(0.3, 0.35, 0.1),
    vec3(1.0)
);

const vec3 kAmbientColor = vec3(0.376, 0.0, 0.10);

vec3 hash( vec3 p ) {
    p = vec3( dot(p,vec3(127.1,311.7, 74.7)),
              dot(p,vec3(269.5,183.3,246.1)),
              dot(p,vec3(113.5,271.9,124.6)));

    return fract(sin(p)*43758.5453123);
}

RayResult opUnion(in RayResult a, in RayResult b) {
    if (a.dist < b.dist) return a;
    return b;
}

RayResult opUnion2(in RayResult a, in RayResult b) {
    if (-a.dist < b.dist) return a;
    return b;
}

float sdRoundBox( vec3 p, vec3 b, float r ) {
    vec3 d = abs(p) - b;
    return length(max(d,0.0)) - r
            + min(max(d.x,max(d.y,d.z)),0.0); // remove this line for an only partially signed sdf
}

float sdSphere(in vec3 p, float r) {
    return length(p) - r;
}

float sdBox( vec3 p, vec3 b ) {
    vec3 d = abs(p) - b;
    return length(max(d,0.0))
            + min(max(d.x,max(d.y,d.z)),0.0); // remove this line for an only partially signed sdf
}

float sdTorus( vec3 p, vec2 t ) {
    vec2 q = vec2(length(p.xz)-t.x,p.y);
    return length(q)-t.y;
}

float sdPlane( vec3 p, vec4 n ) {
    // n must be normalized
    return dot(p,n.xyz) + n.w;
}

float opDisp(vec3 p) {
    return sin(20.0*p.x)*sin(20.0*p.y)*sin(20.0*p.z);
}

void opRotate(inout vec2 v, float r) {
    float c = cos(r);
    float s = sin(r);
    float vx = v.x * c - v.y * s;
    float vy = v.x * s + v.y * c;
    v.x = vx;
    v.y = vy;
}

RayResult mapScene(in vec3 p) {
	p.y += 0.2;
	float t = pc.time * 0.5;
	float a = sdSphere(p + vec3(1.5 * cos(t), -0.9, 1.5 * sin(t)), 0.8);
	a += opDisp(p * 0.3) * 0.19;
	vec3 bp = p + vec3(1.5 * cos(t + 2.5), -1.02, 1.5 * sin(t + 2.5));
	opRotate(bp.xz, t * 1.8);

	float b = sdBox(bp, vec3(0.88));
	b -= opDisp(sin(pc.time) * p * 0.2) * 0.3;

	float c = sdBox(p + vec3(0.0, 0.2, 0.0), vec3(4.0, 0.2, 4.0));

	RayResult plane = RayResult(c, kMaterialBlue);
	RayResult sphere = RayResult(a, kMaterialRed);
	RayResult box = RayResult(b, kMaterialGreen);

	RayResult res = opUnion(plane, opUnion(sphere, box));

    return res;
}

RayResult rayMarch(in vec3 ro, in vec3 rd) {
    float total = NEARCLIP;
    for (int i = 0; i < MARCHSTEPS; ++i) {
        RayResult ray = mapScene(ro + rd * total);
        if (ray.dist < EPSILON) {
            return RayResult(total, ray.material);
        }

        total += ray.dist * 0.5;
        if (total > FARCLIP) {
            break;
        }
    }

    return RayResult(FARCLIP, kNoMaterial);
}

vec3 normal(in vec3 p) {
    const vec2 e = vec2(0.0, EPSILON);
    return normalize(vec3(
        mapScene(p + e.yxx).dist - mapScene(p - e.yxx).dist,
        mapScene(p + e.xyx).dist - mapScene(p - e.xyx).dist,
        mapScene(p + e.xxy).dist - mapScene(p - e.xxy).dist
    ));
}

float opSSS(in vec3 ro, in vec3 rd, in vec3 n, float dist, float factor) {
	float value = 0.0;
	vec3 nrd = refract(rd, n, 1.0);
	float s = 1.0;
	const int steps = int(factor);

	for (int i = 0; i < steps; ++i) {
		float step_size = (float(i) / float(steps));
		float d = mapScene(ro + nrd * step_size).dist;
		value += (step_size - d) * s;
		s *= 0.6;
	}
	value = pow(value, 0.2);
	value = clamp(abs(1.0 - value), 0.0, 1.0000001);
	return value;
}

vec3 opReflection(float dist, in vec3 p, in vec3 dir, in vec3 n) {
	vec3 color = vec3(0.2);
	vec3 rd = normalize(reflect(dir, n));
	float ft = max(0.0, dot(rd, n));
	RayResult ray = rayMarch(p, rd);

	if (ray.dist < FARCLIP && ft > 0.0) {
		Material material = ray.material;
		vec3 hit_point = p + rd * ray.dist;
		vec3 n = normal(hit_point);
		vec3 nld = normalize(kDirLight.dir);
		vec3 h = normalize(n + nld);
		float diffuse = max(0.0, dot(n, nld));
		float specular = pow(max(0.0, dot(h, n)), material.shininess);
		vec3 sss = opSSS(hit_point, nld, n, dist, material.sss) *
				    material.diffuse;

		if (material.shininess == 0.0) specular = 0.0;

		vec3 color = (material.ambient) / 2.0;

		color +=
            ((material.diffuse + kDirLight.color) / 2.0) * diffuse +
            ((material.specular * specular));

		color += sss * clamp(material.sss, 0.0, 1.0);

		return color;
	}

	return kAmbientColor;
}

float opAO(in vec3 p, in vec3 n) {
	float value = 0.0;
	float s = 1.0;
	for(int i = 0; i < 3; ++i) {
		float step_size = 0.13;
		float dist = mapScene(p + n * step_size).dist;
		value += (step_size - dist) * s;
		s *= 0.7;
	}
	return clamp(sqrt((0.9 - value) * sqrt(1.0)), -1.0, 1.0);
}

float opShadow(in vec3 ro, in vec3 rd, in float mint, in float tmax) {
	float res = 1.0;
	float t = mint;
	for (int i = 0; i < 32; ++i) {
		float h = mapScene(ro + rd * t).dist;
		res - min(res, 2.0 * h / t);
		t += clamp(h, 0.02, 0.1);
		if (h < 0.001 || t > tmax) break;
	}
	return clamp(res, 0.0, 1.0);
}

vec3 shade(float dist, in vec3 ro, in vec3 rd, in vec3 hit_point, in Material material) {
	vec3 n = normal(hit_point);
	vec3 nld = normalize(kDirLight.dir);
	vec3 h = normalize(n + nld);
	float diffuse = max(0.0, dot(n, nld));
	float specular = pow(max(0.0, dot(h, n)), material.shininess);

	float ao = opAO(hit_point, n);
	float shadow = opShadow(hit_point, nld, 0.1, 0.9);
	vec3 reflection_color = opReflection(dist, hit_point, rd, n);
    vec3 sss = opSSS(hit_point, nld, n, dist, material.sss) * (material.diffuse * 1.2);

	vec3 color = n;

	if (material.shininess == 0.0) specular = 0.0;

	color = (material.ambient + kAmbientColor.rgb) / 2.0 +
        ((material.diffuse + kDirLight.color) / 2.0) * diffuse +
        ((material.specular * specular));

	color = mix(color, clamp(reflection_color, vec3(0.0), vec3(1.0)), clamp(material.reflection, 0.0, 1.0));

	color *= ao;

	color *= shadow * 0.2 + 0.7;

	color += sss * clamp(material.sss, 0.0, 1.0);

	return color;
}

void main() {
	vec2 aspect =  vec2(pc.resolution.x / pc.resolution.y, 1);
    vec2 uv = (in_uv + -0.5) * 2.0 * aspect;
	vec3 color = kAmbientColor.rgb;
	vec3 ro = vec3(1.0, 4.0, -6.5);
	vec3 rd = normalize(vec3(uv, 1.0));
	opRotate(rd.xz, 0.1);
	opRotate(rd.yz, 0.4);
    opRotate(kDirLight.dir.xz, pc.time * 1.5);

	RayResult ray = rayMarch(ro, rd);

	vec3 hit_point = ro + rd * ray.dist;

	if (ray.dist < FARCLIP) {
		color = shade(ray.dist, ro, rd, hit_point, ray.material);
	}

	color = pow(color, vec3(1./0.7));

    color = mix(color, color * 0.65,
                max(0.0, pow(length(1.0 * uv / aspect), 2.9)));

    out_color = vec4(color, 1.0);
}
