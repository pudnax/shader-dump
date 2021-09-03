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

vec3 standardObserver1931[] =
    vec3[] (
    vec3( 0.001368, 0.000039, 0.006450 ), // 380 nm
    vec3( 0.002236, 0.000064, 0.010550 ), // 385 nm
    vec3( 0.004243, 0.000120, 0.020050 ), // 390 nm
    vec3( 0.007650, 0.000217, 0.036210 ), // 395 nm
    vec3( 0.014310, 0.000396, 0.067850 ), // 400 nm
    vec3( 0.023190, 0.000640, 0.110200 ), // 405 nm
    vec3( 0.043510, 0.001210, 0.207400 ), // 410 nm
    vec3( 0.077630, 0.002180, 0.371300 ), // 415 nm
    vec3( 0.134380, 0.004000, 0.645600 ), // 420 nm
    vec3( 0.214770, 0.007300, 1.039050 ), // 425 nm
    vec3( 0.283900, 0.011600, 1.385600 ), // 430 nm
    vec3( 0.328500, 0.016840, 1.622960 ), // 435 nm
    vec3( 0.348280, 0.023000, 1.747060 ), // 440 nm
    vec3( 0.348060, 0.029800, 1.782600 ), // 445 nm
    vec3( 0.336200, 0.038000, 1.772110 ), // 450 nm
    vec3( 0.318700, 0.048000, 1.744100 ), // 455 nm
    vec3( 0.290800, 0.060000, 1.669200 ), // 460 nm
    vec3( 0.251100, 0.073900, 1.528100 ), // 465 nm
    vec3( 0.195360, 0.090980, 1.287640 ), // 470 nm
    vec3( 0.142100, 0.112600, 1.041900 ), // 475 nm
    vec3( 0.095640, 0.139020, 0.812950 ), // 480 nm
    vec3( 0.057950, 0.169300, 0.616200 ), // 485 nm
    vec3( 0.032010, 0.208020, 0.465180 ), // 490 nm
    vec3( 0.014700, 0.258600, 0.353300 ), // 495 nm
    vec3( 0.004900, 0.323000, 0.272000 ), // 500 nm
    vec3( 0.002400, 0.407300, 0.212300 ), // 505 nm
    vec3( 0.009300, 0.503000, 0.158200 ), // 510 nm
    vec3( 0.029100, 0.608200, 0.111700 ), // 515 nm
    vec3( 0.063270, 0.710000, 0.078250 ), // 520 nm
    vec3( 0.109600, 0.793200, 0.057250 ), // 525 nm
    vec3( 0.165500, 0.862000, 0.042160 ), // 530 nm
    vec3( 0.225750, 0.914850, 0.029840 ), // 535 nm
    vec3( 0.290400, 0.954000, 0.020300 ), // 540 nm
    vec3( 0.359700, 0.980300, 0.013400 ), // 545 nm
    vec3( 0.433450, 0.994950, 0.008750 ), // 550 nm
    vec3( 0.512050, 1.000000, 0.005750 ), // 555 nm
    vec3( 0.594500, 0.995000, 0.003900 ), // 560 nm
    vec3( 0.678400, 0.978600, 0.002750 ), // 565 nm
    vec3( 0.762100, 0.952000, 0.002100 ), // 570 nm
    vec3( 0.842500, 0.915400, 0.001800 ), // 575 nm
    vec3( 0.916300, 0.870000, 0.001650 ), // 580 nm
    vec3( 0.978600, 0.816300, 0.001400 ), // 585 nm
    vec3( 1.026300, 0.757000, 0.001100 ), // 590 nm
    vec3( 1.056700, 0.694900, 0.001000 ), // 595 nm
    vec3( 1.062200, 0.631000, 0.000800 ), // 600 nm
    vec3( 1.045600, 0.566800, 0.000600 ), // 605 nm
    vec3( 1.002600, 0.503000, 0.000340 ), // 610 nm
    vec3( 0.938400, 0.441200, 0.000240 ), // 615 nm
    vec3( 0.854450, 0.381000, 0.000190 ), // 620 nm
    vec3( 0.751400, 0.321000, 0.000100 ), // 625 nm
    vec3( 0.642400, 0.265000, 0.000050 ), // 630 nm
    vec3( 0.541900, 0.217000, 0.000030 ), // 635 nm
    vec3( 0.447900, 0.175000, 0.000020 ), // 640 nm
    vec3( 0.360800, 0.138200, 0.000010 ), // 645 nm
    vec3( 0.283500, 0.107000, 0.000000 ), // 650 nm
    vec3( 0.218700, 0.081600, 0.000000 ), // 655 nm
    vec3( 0.164900, 0.061000, 0.000000 ), // 660 nm
    vec3( 0.121200, 0.044580, 0.000000 ), // 665 nm
    vec3( 0.087400, 0.032000, 0.000000 ), // 670 nm
    vec3( 0.063600, 0.023200, 0.000000 ), // 675 nm
    vec3( 0.046770, 0.017000, 0.000000 ), // 680 nm
    vec3( 0.032900, 0.011920, 0.000000 ), // 685 nm
    vec3( 0.022700, 0.008210, 0.000000 ), // 690 nm
    vec3( 0.015840, 0.005723, 0.000000 ), // 695 nm
    vec3( 0.011359, 0.004102, 0.000000 ), // 700 nm
    vec3( 0.008111, 0.002929, 0.000000 ), // 705 nm
    vec3( 0.005790, 0.002091, 0.000000 ), // 710 nm
    vec3( 0.004109, 0.001484, 0.000000 ), // 715 nm
    vec3( 0.002899, 0.001047, 0.000000 ), // 720 nm
    vec3( 0.002049, 0.000740, 0.000000 ), // 725 nm
    vec3( 0.001440, 0.000520, 0.000000 ), // 730 nm
    vec3( 0.001000, 0.000361, 0.000000 ), // 735 nm
    vec3( 0.000690, 0.000249, 0.000000 ), // 740 nm
    vec3( 0.000476, 0.000172, 0.000000 ), // 745 nm
    vec3( 0.000332, 0.000120, 0.000000 ), // 750 nm
    vec3( 0.000235, 0.000085, 0.000000 ), // 755 nm
    vec3( 0.000166, 0.000060, 0.000000 ), // 760 nm
    vec3( 0.000117, 0.000042, 0.000000 ), // 765 nm
    vec3( 0.000083, 0.000030, 0.000000 ), // 770 nm
    vec3( 0.000059, 0.000021, 0.000000 ), // 775 nm
    vec3( 0.000042, 0.000015, 0.000000 )  // 780 nm
);
float standardObserver1931_w_min = 380.0f;
float standardObserver1931_w_max = 780.0f;
int standardObserver1931_length = 81;
vec3 WavelengthToXYZLinear(float fWavelength) {
    float fPos = (fWavelength - standardObserver1931_w_min) /
                 (standardObserver1931_w_max - standardObserver1931_w_min);
    float fIndex = fPos * float(standardObserver1931_length);
    float fFloorIndex = floor(fIndex);
    float fBlend = clamp(fIndex - fFloorIndex, 0.0, 1.0);
    int iIndex0 = int(fFloorIndex);
    int iIndex1 = iIndex0 + 1;
    iIndex1 = min(iIndex1, standardObserver1931_length - 1);

    return mix(standardObserver1931[iIndex0], standardObserver1931[iIndex1],
               fBlend);
}

float DispersionLaw(float wv) {
    wv = wv * 1e-3;  // micro meters
    const float A = 1.7280;
    const float B = 0.2;
    return A + B / (wv * wv);
}

float BlackBody(float t, float w_nm) {
    float h = 6.6e-34;
    float k = 1.4e-23;
    float c = 3e8;

    float w = w_nm / 1e9;

    float w5 = w * w * w * w * w;
    float o = 2. * h * (c * c) / (w5 * (exp(h * c / (w * k * t)) - 1.0));

    return o;
}

mat2 rot(float a) {
    float c = cos(a), s = sin(a);
    return mat2(c, -s, s, c);
}

float box(vec3 p, vec3 s) {
    p = abs(p) - s;
    /* return max(p.x, max(p.y, p.z)); */
    return length(max(p, 0.0)) + min(max(p.x, max(p.y, p.z)), 0.0);
}

float caps(vec3 p, vec3 p1, vec3 p2, float s) {
    vec3 pa = p - p1;
    vec3 pb = p2 - p1;
    float prog = dot(pa, pb) / dot(pb, pb);
    prog = clamp(prog, 0., 1.);
    return length(p1 + pb * prog - p) - s;
}

int scene = 0;

float map(vec3 p) {
    vec3 p2 = p;
    float t = pc.time * 0.1;
    p2.yz *= rot(t);
    p2.yx *= rot(t * 1.3);
    float d4 = max(box(p2, vec3(3)), 1.2 - length(p));

	vec3 pb = p;
	float d2 = 10000.;
	for (float i = 0.; i < 3.; ++i) {
		float t = pc.time * 0.03 + i;
		p.yz *= rot(t + i);
		p.yx *= rot(t * 1.3);
		d2 = min(d2, length(p) - 0.47);
		p = abs(p);
		p -= 0.9;
	}

    float d = box(p, vec3(0.4, 0.4, 0.4));

	d = min(d, d2);

	if (scene == 0) d = d4;

    return d;
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


#define PATTERN_RAND

const float PI = acos(-1.);
const float TWO_PI = 2. * PI;

vec4 PatternRand( uint seed )
{
    return vec4(
        float((seed*0x73494U)&0xfffffU)/float(0x100000),
    	float((seed*0xAF71FU)&0xfffffU)/float(0x100000),
        float((seed*0x67a42U)&0xfffffU)/float(0x100000), // a bit stripey against x and z, but evens out over time
        float((seed*0x95a8cU)&0xfffffU)/float(0x100000) // good vs x or y, not good against z
    );
}

uvec4 Hash(uint seed) {
    // integer hash from Hugo Elias
    seed = (seed << 13U) ^ seed;
    seed = seed * (seed * seed * 15731U + 789221U) + 1376312589U;
    return seed * uvec4(seed, seed * 16807U, seed * 48271U, seed * 31713U);
}

vec4 Rand(uint seed) {
#if defined(PATTERN_RAND)
    return PatternRand(seed);
#else
    return vec4(Hash(seed) & 0x7fffffffU) / float(0x7fffffffU);
#endif
}

#define MAX_LEVEL 5
#define bayer2x2(a) (4-(a).x-((a).y<<1))%4
float GetBayerFromCoordLevel(vec2 pixelpos)
{
    ivec2 ppos = ivec2(pixelpos);
    int sum = 0;
    for(int i=0; i<MAX_LEVEL; i++)
    {
        sum += bayer2x2(ppos>>(MAX_LEVEL-1-i)&1)<<(2*i);
    }

    return float(sum) / float(2<<(MAX_LEVEL*2-1));
}

float rnd(vec2 uv) {
    return fract(dot(sin(uv * 452.714 + uv.yx * 547.524), vec2(352.887)));
}

vec3 norm(vec3 p) {
  mat3 k = mat3(p,p,p) - mat3(0.001);
  return normalize(map(p) - vec3(map(k[0]),map(k[1]),map(k[2])));
}

vec3 scol = vec3(0.);

void main() {
    vec2 uv = (in_uv + -0.5) * 2.0 / vec2(pc.resolution.y / pc.resolution.x, 1);

    /* scene = int(mod(floor(pc.time / 10.), 2.)); */

    vec3 s2 = vec3(10, 0., 0.);
    vec3 r2 = normalize(vec3(-1, sin(pc.time) * 0.1, 0));

	vec2 frag_coord = in_uv*pc.resolution * 2.1;

    float motionJitter = GetBayerFromCoordLevel(in_uv*pc.resolution);
    /* motionJitter = Hash(uint(uv.x * uv.y)).x; */

    uint pseed = /*uint(pc.frame<<16)*/ + (uint(frag_coord.y)) * uint(frag_coord.x);
    vec4 prand = Rand( pseed );
	/* motionJitter = prand.y; */

    float rnd1 = rnd(uv + fract(pc.time * 0.1));
    float ior = (rnd1 - 0.5);
    ior = motionJitter - 0.5;
	/* ior = pc.pos.x; */
    float id = ior * 2.0;
    vec3 diff = 1.3 - vec3(1. + id, 0.454 + abs(id), 1. - id);

    float wv = mix(420.0, 750.0, rnd1 - 0.5);
    wv = mix(420.0, 750.0, motionJitter + 0.5);
    float id1 = DispersionLaw(wv);
    float blackb = 1e-13 * BlackBody(6000.0, wv);
    scol = blackb * WavelengthToXYZLinear(wv);

    ior = id1;
    diff = 1.3 - vec3(1. + id, 0.454 + abs(id), 1. - id);

    vec3 p2 = s2;
    points[0] = p2;
    pid = 1;
    float side = 1;
    for (int i = 0; i < 60; ++i) {
        float d = abs(map(p2));
        if (d < 0.001) {
            points[pid] = p2;
            pid += 1;
            if (pid >= pcount - 1) break;

            vec2 off = vec2(0.01, 0);
            vec3 n2 =
                side * normalize(d - vec3(map(p2 - off.xyy), map(p2 - off.yxy),
                                          map(p2 - off.yyx)));

            vec3 r3 = refract(r2, n2, 1. - side * (0.3 + 0.1 * ior));
            if (dot(r3, r3) < 0.5) r3 = reflect(r2, n2);
            r2 = r3;
            side = -side;
            d = 0.1;
        }
        if (d > 100.0) break;
        p2 += r2 * d;
    }
    points[pid] = p2 + r2 * 1000.;
    pid += 1;

    vec3 s = vec3(0., 0., -10.);
    vec3 r = normalize(vec3(uv, 1.));

    float mumu = mix(rnd(-uv + fract(pc.time * 0.1)), 1., 0.9);
	mumu = mix(motionJitter, 1., 0.9);
    vec3 p = s;
    float side2 = 1.;
    for (int i = 0; i < 90; ++i) {
        float d = abs(map2(p));
        if (d < 0.001) {
            vec2 off = vec2(0.01, 0);
            vec3 n =
                side2 * normalize(d - vec3(map(p - off.xyy), map(p - off.yxy),
                                           map(p - off.yyx)));
            vec3 r3 = refract(r, n, 1. - side2 * (0.3 + 0.1 * ior));
            if (dot(r3, r3) < 0.5) r3 = reflect(r, n);
            r = r3;

            side2 = -side2;
            d = 0.1;
        }
        if (d > 100.) break;
        p += r * d * mumu;
    }

    vec3 col = vec3(0.);
    vec3 lazer = diff * atm;

    vec3 light_dir = normalize(vec3(1));
    float shade = dot(light_dir, norm(p));
    float amb = 0.2 * (mix(max(shade, 0.), shade * 0.5 + 0.5, .05));

    /* col = mix(lazer, vec3(amb), motionJitter); */
    col = mix(lazer, vec3(amb), 0.);

    /* col = vec3(amb); */

    col = smoothstep(0.01, 0.9, col);
    /* col = pow(col, vec3(0.4545)); */

    out_color = vec4(col, 1.0);
}
