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

#define EPS vec3(0., 0.001, 0.0001)
#define PI acos(-1.)

vec3 erot(vec3 p, vec3 ax, float ro) {
    return mix(dot(ax,p)*ax, p, cos(ro)) + sin(ro)*cross(ax,p);
}

float smin2( float a, float b, float k ) {
    float h = max( k-abs(a-b), 0.0 )/k;
    return min( a, b ) - h*h*k*(1.0/4.0);
}

float smin( float a, float b, float k ) {
    float res = exp2( -k*a ) + exp2( -k*b );
    return -log2( res )/k;
}

float smax(float a, float b, float k) {
    return smin(a, b, -k);
}

float random(vec3 p) {
    float x = dot(p, vec3(4371.321, -9137.327, 235235.23423));
    return 2.0 * fract(sin(x) * 17381.94472) - 1.0;
}

float noise(in vec3 p) {
    vec3 id = floor(p);
    vec3 f = fract(p);

    vec3 u = f * f * (3.0 - 2.0 * f);

    return mix(mix(random(id + vec3(0.0)), random(id + vec3(1.0, 0.0, 0.0)), u.x),
             mix(random(id + vec3(0.0, 1.0, 1.0)), random(id + vec3(1.0, 1.0, 1.0)), u.x),
             u.y);
}

float fbm( in vec3 x) {
    float G = exp2(-1);
    float f = 1.0;
    float a = 1.0;
    float t = 0.0;
    for( int i=0; i< 4; i++ )
    {
        t += a*noise(f*x);
        f *= 2.0;
        a *= G;
    }
    return t;
}

float noise_fbm(vec3 p) {
    float h =
        fbm(0.09 * pc.time + p + fbm(0.065 * pc.time + 2.0 * p - 5.0 * fbm(4.0 * p)));
    return h;
}

float sdSphere(in vec3 p, in float d) {
	return length(p) - d;
}

float sdTorus(in vec3 p, in vec2 t) {
    vec2 q = vec2( length( p.xz ) - t.x, p.y );
    return length( q ) - t.y;
}

vec2 rotate(in vec2 p, in float t) {
    return p * cos(-t) + vec2(p.y, -p.x) * sin(-t);
}

vec3 rotate(in vec3 p, in vec3 t) {
  p.yz = rotate(p.yz, t.x);
  p.zx = rotate(p.zx, t.y);
  p.xy = rotate(p.xy, t.z);
  return p;
}

float map(vec3 p) {
	float size = 2.0;
	vec3 ps = p;
	ps = rotate(ps, vec3(0.2, pc.time*0.2, -pc.time*0.3));
	float sphere = sdSphere(ps, size);
	sphere = max(
			sphere,
            abs(fract(ps.x + pc.time * 0.3)) - .8
			);

	float size_tor = size * 1.5;
	vec3 p1 = p;
	p1.xy = rotate(p1.xy, -pc.time * 0.2);
	float torus1 = sdTorus(p1, vec2(size_tor, 0.2));

	vec3 p2 = p;
	p2.xy = rotate(p2.xy, pc.time * 0.2);
	p2.xy = rotate(p2.xy, PI/2);
	float torus2 = sdTorus(p2, vec2(size_tor, 0.2));

	float torus = smin(torus1, torus2, 1.7);

	float res = min(sphere, torus);

	return res;
}

vec3 norm(vec3 p) {
	mat3 k = mat3(p,p,p) - mat3(EPS.y);
	return normalize(map(p) - vec3(map(k[0]),map(k[0]),map(k[0])));
}

vec3 hsv(float h, float s, float v) {
  return mix( vec3( 1.0 ), clamp( ( abs( fract(
    h + vec3( 3.0, 2.0, 1.0 ) / 3.0 ) * 6.0 - 3.0 ) - 1.0 ), 0.0, 1.0 ), s ) * v;
}

void main() {
    vec2 uv = (in_uv + -0.5) * 2.0 / vec2(pc.resolution.y / pc.resolution.x, 1);
	vec2 mouse = pc.mouse / vec2(pc.resolution.y / pc.resolution.x, 1);
	if (distance(mouse, uv) < 0.5) {
        uv *= 0.5 + dot(uv * .2, uv * .3);
	}

	vec3 rd = normalize(vec3(uv, -1.8));
    vec3 ro = vec3(0., 0., 8.);

    vec3 rot = vec3(1.4, 0.2, 0.1);
    rot = vec3(0.5 + pc.time * 0.2, 0.3 + pc.time * 0.2, 0.2);
	ro = rotate(ro, rot);
    rd = rotate(rd, rot);

    float s = fract(abs(sin(pc.time * 0.1))) * 8 + 2;
	s = 10;
    ro *= s;

	vec3 grid = floor(ro);
	vec3 grid_step = sign( rd );
	vec3 delta = ( -fract(ro) + 0.5 * ( grid_step + 1.0 ) ) / rd;
	vec3 delta_step =  1.0 / abs( rd );
	vec3 mask = vec3(.0);

	float dist = EPS.y;
	bool hit = false;
	vec3 pos = vec3(EPS.y);
	for (int i = 0; i < 200; ++i) {
	    pos = (grid + 0.5) / s;
	    if (map(pos) < 0.0) {
	        hit = true; break;
	    }
	    vec3 c = step(delta, delta.yzx);
	    mask = c * (1.0 - c.zxy);
	    grid += grid_step * mask;
	    delta += delta_step * mask;
	}
    vec3 col = vec3(0.3 + 0.15 * uv.x);
	if (hit) {
		vec3 normal =  norm(pos);
		normal = erot(normal, normalize(vec3(-1, 1, 0)), .96);
		float light = dot(max(normal, 0.), vec3(.4));

		vec3 p = pos;
		col = hsv(0.3 * length(pos), 0.5, 0.8);
		col += vec3(light) * 0.2;

		float br = dot(vec3(0.5,0.9,0.7), mask);

		float depth = dot(delta - delta_step, mask);
		float fog = min(1.0, 1000.0 / depth / depth);

		vec3 uvw = fract(ro + rd * depth);
		vec2 uu = vec2(dot(uvw.yzx, mask), dot(uvw.zxy, mask));
		uu = abs(uu - vec2(0.5));
		float gr = 1.0 - 0.1 * smoothstep(0.4, 0.5, max(uu.x, uu.y));

		col *= br * fog * gr;
	}

    out_color = vec4(col, 1.0);
}
