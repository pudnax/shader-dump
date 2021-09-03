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

vec3 hash3(float n) {
    return fract(sin(vec3(n, n + 1.0, n + 2.0)) * 43758.5453123);
}

vec4 hash4(float n) {
    return fract(sin(vec4(n, n + 1.0, n + 2.0, n + 3.0)) * 43758.5453123);
}

const float speed = 5.0;
vec2 getPos(float time, vec4 id) {
    return vec2(0.9 * sin((speed * (0.75 + 0.5 * id.z)) * time + 20.0 * id.x),
                0.75 * cos(speed * (0.75 + 0.5 * id.w) * time + 20.0 * id.y));
}

vec2 getVelocity(float time, vec4 id) {
    return vec2(
        speed * 0.9 * cos((speed * (0.75 + 0.5 * id.z)) * time + 20.0 * id.x) *
            (0.75 + 0.5 * id.z),
        -speed * 0.75 * sin(speed * (0.75 + 0.5 * id.w) * time + 20.0 * id.y) *
            (0.75 + 0.5 * id.w));
}

vec3 disk(vec3 col, in vec2 uv, in vec3 sph, in vec3 sphcol) {
    vec2 xc = uv - sph.xy;
    float h = (sph.z - sqrt(dot(xc, xc))) / sph.z;
    return mix(col, sphcol, clamp(2.0 * h, 0.0, 1.0));
}

vec3 diskWithMotionBlur(vec3 col,
                        in vec2 uv,
                        in vec3 sph,
                        in vec2 cd,
                        in vec3 sphcol) {
    vec2 xc = uv - sph.xy;
    float a = dot(cd, cd);
    float b = dot(cd, xc);
    float c = dot(xc, xc) - sph.z * sph.z;
    float h = b * b - a * c;
    if (h > 0.0)  // determine if there is a real root
    {
        h = sqrt(h);
        float ta = max(0.0, (-b - h) / a);  // Solve with the rooting formula
        float tb = min(1.0, (-b + h) / a);
        if (ta < tb)
            col = mix(col, sphcol, clamp(2.0 * (tb - ta), 0.0, 1.0));
    }
    return col;
}

const int ballNum = 15;
void main() {
    vec2 uv = (in_uv + -0.5) * 2.0 / vec2(pc.resolution.y / pc.resolution.x, 1);
    vec2 p = uv;

    vec3 col = vec3(0.2) + 0.05 * p.y;

    for (int i = 0; i < ballNum; i++) {
        vec4 off = hash4(float(i) * 13.13);
        vec3 sph = vec3(getPos(pc.time, off), 0.02 + 0.1 * off.x);
        vec2 cd = getVelocity(pc.time, off) / 24.0;

        vec3 sphcol = 0.7 + 0.3 * sin(3.0 * off.z + vec3(4.0, 0.0, 2.0));

        col = diskWithMotionBlur(col, p, sph, cd, sphcol);
    }

    col += (1.0 / 255.0) * hash3(p.x + 13.0 * p.y);
    out_color = vec4(col, 1.0);
}
