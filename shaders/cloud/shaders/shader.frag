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

mat3 m = mat3( 0.00,  0.80,  0.60,
              -0.80,  0.36, -0.48,
              -0.60, -0.48,  0.64);

float hash(float n) {
    return fract(sin(n) * 43758.5453);
}

float noise(in vec3 x) {
    vec3 p = floor(x);
    vec3 f = fract(x);

    f = f * f * (3.0 - 2.0 * f);

    float n = p.x + p.y * 57.0 + 113.0 * p.z;

    float res = mix(mix(mix(hash(n + 0.0), hash(n + 1.0), f.x),
                        mix(hash(n + 57.0), hash(n + 58.0), f.x), f.y),
                    mix(mix(hash(n + 113.0), hash(n + 114.0), f.x),
                        mix(hash(n + 170.0), hash(n + 171.0), f.x), f.y),
                    f.z);
    return res;
}

float fbm(vec3 p) {
    float f;
    f  = 0.5000 * noise(p); p = m * p * 2.02;
    f += 0.2500 * noise(p); p = m * p * 2.03;
    f += 0.1250 * noise(p);
    return f;
}

float scene(in vec3 pos) {
    return 0.1 - length(pos) * 0.05 + fbm(pos * 0.3);
}

vec3 getNormal(in vec3 p) {
    const float e = 0.01;
    return normalize(vec3(scene(vec3(p.x + e, p.y, p.z)) - scene(vec3(p.x - e, p.y, p.z)),
                          scene(vec3(p.x, p.y + e, p.z)) - scene(vec3(p.x, p.y - e, p.z)),
                          scene(vec3(p.x, p.y, p.z + e)) - scene(vec3(p.x, p.y, p.z - e))));
}

mat3 camera(vec3 ro, vec3 ta) {
    vec3 cw = normalize(ta - ro);
    vec3 cp = vec3(0.0, 1.0, 0.0);
    vec3 cu = cross(cw, cp);
    vec3 cv = cross(cu, cw);
    return mat3(cu, cv, cw);
}

void main() {
    vec2 uv = (in_uv + -0.5) * 2.0 / vec2(pc.resolution.y / pc.resolution.x, 1);

	vec2 mo = vec2(pc.time * 0.1, cos(pc.time * 0.25) * 3.0);
    float camDist = 35.0;
    vec3 ta = vec3(0.0, 1.0, 0.0);

    vec3 ro = camDist * normalize(vec3(cos(2.75 - 3.0 * mo.x), 0.7 - 1.0 * (mo.y - 1.0), sin(2.75 - 3.0 * mo.x)));
    float targetDepth = 1.0;
    mat3 c = camera(ro, ta);
    vec3 dir = c * normalize(vec3(uv, targetDepth));

    const int sampleCount = 64;
    const int sampleLightCount = 6;
    const float eps = 0.01;
    float zMax = 40.0;
    float zstep = zMax / float(sampleCount);
    float zMaxl = 20.0;
    float zstepl = zMaxl / float(sampleLightCount);
    vec3 p = ro;
    float T = 1.0;
    float absorption = 100.0;
    vec3 sun_direction = normalize(vec3(1.0, 0.0, 0.0));
    vec3 color = vec3(0.0);

	for (int i = 0; i < sampleCount; ++i) {
		float density = scene(p);
		if (density > 0.0) {
			float tmp = density / float(sampleCount);
			T *= 1.0 - (tmp * absorption);
			if (T <= 0.01) break;
			float Tl = 1.0;
			vec3 lp = p;
			for (int j = 0; j < sampleLightCount; ++j) {
				float density_light = scene(lp);
				if (density_light > 0.0) {
					float tmpl = density_light / float(sampleCount);
					Tl *= 1.0 - (tmpl * absorption);
				}
				if (Tl <= 0.01) break;

				lp += sun_direction * zstepl;
			}

			float opacity = 50.0;
			float k = opacity * tmp * T;
			vec3 cloud_color = vec3(1.0);
			vec3 col1 = cloud_color * k;

			float opacityl = 30.;
			float kl = opacityl * tmp * T * Tl;
			vec3 light_color = vec3(1.0, 0.7, 0.9);
			vec3 col2 = light_color * kl;

			color += col1 + col2;
		}
		p += dir * zstep;
	}

	vec3 bg = mix(vec3(0.3, 0.1, 0.8), vec3(0.7, 0.7, 1.0), 1.0 - (uv.y + 1.0) * 0.5);
    color += bg;

    out_color = vec4(color, 1.0);
}
