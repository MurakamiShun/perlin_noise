#pragma once
#include <cmath>

unsigned long xorshift32(unsigned long x) {
	x = x ^ (x << 13);
	x = x ^ (x >> 17);
	x = x ^ (x << 15);
	return x;
}


// 1D

double wavelet1d(double x, double grad) {
	double x_2 = x * x;
	double x_4 = x_2 * x_2;
	return (1 - (6 * std::abs(x_4*x) - 15 * x_4 + 10 * std::abs(x_2*x))) // C(x)
		* (grad * x); // L(x)|ax
}

double pos2grad1d(unsigned long x, unsigned long seed = 0) {
	double grad;
	unsigned long r = xorshift32(x + seed);
	// grad_x
	r = xorshift32(r + 0) % 10000UL;
	grad = (double)r / 5000.0 - 1.0; // -1.0~1.0
	return grad;
}

double perlin1d(double x, unsigned long seed = 0) {
	double x_p = x - (unsigned long)x;
	unsigned long x_l = (unsigned long)x;
	double wav[2];
	wav[0] = wavelet1d(x_p, pos2grad1d(x_l, seed));
	wav[1] = wavelet1d(x_p - 1.0, pos2grad1d(x_l + 1, seed));
	double result = wav[0] + x_p * (wav[1] - wav[0]);
	return result;
}

// 2D

template<typename T>
struct vec2 {
	T x;
	T y;
};

double wavelet2d(vec2<double> pos, vec2<double> grad) {
	double x_2 = pos.x * pos.x;
	double x_4 = x_2 * x_2;
	double y_2 = pos.y * pos.y;
	double y_4 = y_2 * y_2;
	return (1 - (6 * std::abs(x_4*pos.x) - 15 * x_4 + 10 * std::abs(x_2*pos.x))) // C(x)
		*(1 - (6 * std::abs(y_4*pos.y) - 15 * y_4 + 10 * std::abs(y_2*pos.y))) // C(y)
		*(grad.x*pos.x + grad.y * pos.y); // L(x,y)|ax,ay
}

vec2<double> pos2grad2d(vec2<unsigned long> pos, unsigned long seed = 0) {
	vec2<double> grad;
	unsigned long r = xorshift32(pos.x + seed);
	r = xorshift32(r + pos.y);
	// grad_x
	unsigned long r_dash = xorshift32(r + 0) % 10000UL;
	grad.x = (double)r_dash / 5000.0 - 1.0; // -1.0~1.0
	// grad_y
	r_dash = xorshift32(r + 1) % 10000UL;
	grad.y = (double)r_dash / 5000.0 - 1.0; // -1.0~1.0
	return grad;
}

double perlin2d(vec2<double> pos, unsigned long seed = 0) {
	double x_p = pos.x - (unsigned long)pos.x, y_p = pos.y - (unsigned long)pos.y;
	unsigned long x_l = (unsigned long)pos.x, y_l = (unsigned long)pos.y;
	double wav[2][2];
	wav[0][0] = wavelet2d({ x_p, y_p }, pos2grad2d({ x_l, y_l }, seed));
	wav[0][1] = wavelet2d({ x_p - 1.0, y_p }, pos2grad2d({ x_l + 1,y_l }, seed));
	wav[1][0] = wavelet2d({ x_p, y_p - 1.0 }, pos2grad2d({ x_l,y_l + 1 }, seed));
	wav[1][1] = wavelet2d({ x_p - 1.0, y_p - 1.0 }, pos2grad2d({ x_l + 1,y_l + 1 }, seed));
	double wav_y[2];
	wav_y[0] = wav[0][0] + x_p * (wav[0][1] - wav[0][0]);
	wav_y[1] = wav[1][0] + x_p * (wav[1][1] - wav[1][0]);
	double result = wav_y[0] + y_p * (wav_y[1] - wav_y[0]);
	return result;
}

// 3D

template<typename T>
struct vec3 {
	T x;
	T y;
	T z;
};

double wavelet3d(vec3<double> pos, vec3<double> grad) {
	double x_2 = pos.x * pos.x;
	double x_4 = x_2 * x_2;
	double y_2 = pos.y * pos.y;
	double y_4 = y_2 * y_2;
	double z_2 = pos.z * pos.z;
	double z_4 = z_2 * z_2;
	return (1 - (6 * std::abs(x_4*pos.x) - 15 * x_4 + 10 * std::abs(x_2*pos.x))) // C(x)
		*(1 - (6 * std::abs(y_4*pos.y) - 15 * y_4 + 10 * std::abs(y_2*pos.y))) // C(y)
		*(1 - (6 * std::abs(z_4*pos.z) - 15 * z_4 + 10 * std::abs(z_2*pos.z))) // C(z)
		*(grad.x*pos.x + grad.y * pos.y + grad.z * pos.z); // L(x,y,z)|ax,ay.az
}

vec3<double> pos2grad3d(vec3<unsigned long> pos, unsigned long seed = 0) {
	vec3<double> grad;
	unsigned long r = xorshift32(pos.x + seed);
	r = xorshift32(r + pos.y);
	r = xorshift32(r + pos.z);
	// grad_x
	unsigned long r_dash = xorshift32(r + 0) % 10000UL;
	grad.x = (double)r_dash / 5000.0 - 1.0; // -1.0~1.0
	// grad_y
	r_dash = xorshift32(r + 1) % 10000UL;
	grad.y = (double)r_dash / 5000.0 - 1.0; // -1.0~1.0
	// grad_y
	r_dash = xorshift32(r + 2) % 10000UL;
	grad.z = (double)r_dash / 5000.0 - 1.0; // -1.0~1.0
	return grad;
}

double perlin3d(vec3<double> pos, unsigned long seed = 0) {
	// After decimal point
	double x_p = pos.x - (unsigned long)pos.x;
	double y_p = pos.y - (unsigned long)pos.y;
	double z_p = pos.z - (unsigned long)pos.z;
	// integer
	unsigned long x_l = (unsigned long)pos.x;
	unsigned long y_l = (unsigned long)pos.y;
	unsigned long z_l = (unsigned long)pos.z;
	double wav[2][2][2];
	wav[0][0][0] = wavelet3d({ x_p, y_p, z_p }, pos2grad3d({ x_l, y_l, z_l }, seed));
	wav[0][0][1] = wavelet3d({ x_p - 1.0, y_p, z_p }, pos2grad3d({ x_l + 1,y_l,z_l }, seed));
	wav[0][1][0] = wavelet3d({ x_p, y_p - 1.0, z_p }, pos2grad3d({ x_l,y_l + 1,z_l }, seed));
	wav[0][1][1] = wavelet3d({ x_p - 1.0, y_p - 1.0, z_p }, pos2grad3d({ x_l + 1,y_l + 1,z_l }, seed));
	wav[1][0][0] = wavelet3d({ x_p, y_p, z_p - 1.0 }, pos2grad3d({ x_l, y_l, z_l + 1 }, seed));
	wav[1][0][1] = wavelet3d({ x_p - 1.0, y_p, z_p - 1.0 }, pos2grad3d({ x_l + 1,y_l,z_l + 1 }, seed));
	wav[1][1][0] = wavelet3d({ x_p, y_p - 1.0, z_p - 1.0 }, pos2grad3d({ x_l,y_l + 1,z_l + 1 }, seed));
	wav[1][1][1] = wavelet3d({ x_p - 1.0, y_p - 1.0, z_p - 1.0 }, pos2grad3d({ x_l + 1,y_l + 1,z_l + 1 }, seed));
	double wav_y[2][2];
	wav_y[0][0] = wav[0][0][0] + x_p * (wav[0][0][1] - wav[0][0][0]);
	wav_y[0][1] = wav[0][1][0] + x_p * (wav[0][1][1] - wav[0][1][0]);
	wav_y[1][0] = wav[1][0][0] + x_p * (wav[1][0][1] - wav[1][0][0]);
	wav_y[1][1] = wav[1][1][0] + x_p * (wav[1][1][1] - wav[1][1][0]);
	double wav_z[2];
	wav_z[0] = wav_y[0][0] + y_p * (wav_y[0][1] - wav_y[0][0]);
	wav_z[1] = wav_y[1][0] + y_p * (wav_y[1][1] - wav_y[1][0]);
	double result = wav_z[0] + z_p * (wav_z[1] - wav_z[0]);
	return result;
}