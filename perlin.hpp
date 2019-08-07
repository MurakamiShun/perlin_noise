#pragma once
#include <cmath>
#include <array>

constexpr uint64_t xorshift(uint64_t x) noexcept{
	x = (x + 88172645463325252ULL) ^ (x << 13);
	x = x ^ (x >> 7);
	return x ^ (x << 17);
}
double fade(double t) noexcept {
	t = std::abs(t);
	return 1-(t * t * t * (t * (t * 6 - 15) + 10));
}
constexpr double lerp(const double t, const double a1, const double a2) noexcept {
	return a1 + t * (a2 - a1);
}

// 1D
double wave1d(double x, double grad) noexcept{
	return fade(x) * (grad * x);
}
double pos2grad1d(int64_t x, uint64_t seed = 0)  noexcept {
	return (double)xorshift(x + seed) / (UINT16_MAX / 2.0) - 1.0; // grad // -1.0~1.0
}
double perlin1d(double x, uint64_t seed = 0)  noexcept {
	double x_decimal = x - (int64_t)x;
	int64_t x_int = (int64_t)x;
	return lerp(x_decimal, 
		wave1d(x_decimal, pos2grad1d(x_int, seed)), 
		wave1d(x_decimal - 1.0, pos2grad1d(x_int + 1, seed))
	);
}

// 2D
double wave2d(const std::array<double,2>&& pos, const std::array<double, 2>&& grad) noexcept {
	return fade(pos[0]) * fade(pos[1]) * (grad[0] * pos[0] + grad[1] * pos[1]);
}
std::array<double, 2> pos2grad2d(const std::array<int64_t, 2>&& pos, const uint64_t seed = 0) noexcept {
	constexpr auto mask = UINT16_MAX;
	constexpr auto scale = (double)UINT16_MAX / 2.0;
	uint64_t r = xorshift(pos[0] + seed) + pos[1];
	return {
		(double)(xorshift(r) & mask) / scale - 1.0, // grad_x // -1.0~1.0
		(double)(xorshift(r + 1) & mask) / scale - 1.0 // grad_y // -1.0~1.0
	};
}
double perlin2d(std::array<double, 2> pos, uint64_t seed = 0) noexcept {
	int64_t x_int = (int64_t)pos[0],
		    y_int = (int64_t)pos[1];
	double x_decimal = pos[0] - (int64_t)pos[0],
		   y_decimal = pos[1] - (int64_t)pos[1];
	double wav[2];
	
	wav[0] = lerp(x_decimal,
		wave2d({ x_decimal      , y_decimal }, pos2grad2d({ x_int    , y_int }, seed)),
		wave2d({ x_decimal - 1.0, y_decimal }, pos2grad2d({ x_int + 1, y_int }, seed))
	);
	wav[1] = lerp(x_decimal,
		wave2d({ x_decimal      , y_decimal - 1.0 }, pos2grad2d({ x_int    , y_int + 1 }, seed)),
		wave2d({ x_decimal - 1.0, y_decimal - 1.0 }, pos2grad2d({ x_int + 1, y_int + 1 }, seed))
	);
	return lerp(y_decimal,
		wav[0],
		wav[1]
	);
}

// 3D
double wave3d(const std::array<double,3>&& pos, const std::array<double, 3>&& grad) {
	return fade(pos[0])*fade(pos[1])*fade(pos[2])*
		(grad[0]*pos[0] + grad[1] * pos[1] + grad[2] * pos[2]);
}
std::array<double, 3> pos2grad3d(const std::array<int64_t, 3> pos, uint64_t seed = 0) {
	constexpr auto mask = UINT16_MAX;
	constexpr auto scale = (double)UINT16_MAX / 2.0;
	uint64_t r = xorshift(xorshift(pos[0] + seed) + pos[1]) + pos[2];
	return {
		(double)(xorshift(r) & mask) / scale - 1.0, // grad_x // -1.0~1.0
		(double)(xorshift(r+1) & mask) / scale - 1.0, // grad_y // -1.0~1.0
		(double)(xorshift(r+2) & mask) / scale - 1.0 // grad_z // -1.0~1.0
	};
}
double perlin3d(std::array<double, 3> pos, uint64_t seed = 0) {
	// integer
	int64_t x_int = (int64_t)pos[0];
	int64_t y_int = (int64_t)pos[1];
	int64_t z_int = (int64_t)pos[2];
	// After decimal point
	double x_decimal = pos[0] - x_int;
	double y_decimal = pos[1] - y_int;
	double z_decimal = pos[2] - z_int;

	double wav_y[2][2];
	wav_y[0][0] = lerp(x_decimal,
		wave3d({ x_decimal, y_decimal, z_decimal }, pos2grad3d({ x_int, y_int, z_int }, seed)),
		wave3d({ x_decimal - 1.0, y_decimal, z_decimal }, pos2grad3d({ x_int + 1,y_int,z_int }, seed))
	);
	wav_y[0][1] = lerp(x_decimal,
		wave3d({ x_decimal, y_decimal - 1.0, z_decimal }, pos2grad3d({ x_int,y_int + 1,z_int }, seed)),
		wave3d({ x_decimal - 1.0, y_decimal - 1.0, z_decimal }, pos2grad3d({ x_int + 1,y_int + 1,z_int }, seed))
	);
	wav_y[1][0] = lerp(x_decimal, 
		wave3d({ x_decimal, y_decimal, z_decimal - 1.0 }, pos2grad3d({ x_int, y_int, z_int + 1 }, seed)),
		wave3d({ x_decimal - 1.0, y_decimal, z_decimal - 1.0 }, pos2grad3d({ x_int + 1,y_int,z_int + 1 }, seed))
	);
	wav_y[1][1] = lerp(x_decimal, 
		wave3d({ x_decimal, y_decimal - 1.0, z_decimal - 1.0 }, pos2grad3d({ x_int,y_int + 1,z_int + 1 }, seed)),
		wave3d({ x_decimal - 1.0, y_decimal - 1.0, z_decimal - 1.0 }, pos2grad3d({ x_int + 1,y_int + 1,z_int + 1 }, seed))
	);
	double wav_z[2];
	wav_z[0] = lerp(y_decimal, wav_y[0][0], wav_y[0][1]);
	wav_z[1] = lerp(y_decimal, wav_y[1][0], wav_y[1][1]);
	return lerp(z_decimal, wav_z[0], wav_z[1]);
}


constexpr uint64_t pow_int(const uint64_t t, const uint64_t n) { //t^n
	uint64_t rtn = 1;
	for (uint64_t i = 0; i < n; i++) {
		rtn *= t;
	}
	return rtn;
}
template<uint64_t N>
double wave(const std::array<double, N>& pos, const std::array<double, N>& grad) {
	double mul = fade(pos[0]);
	double dot = pos[0] * grad[0];	
	for (uint64_t n = 1; n < N; n++) {
		mul *= fade(pos[n]);
		dot += pos[n] * grad[n];
	}
	return mul * dot;
}
template<uint64_t N>
std::array<double, N> pos2grad(const std::array<int64_t, N>& pos, uint64_t seed = 0) {
	constexpr auto mask = UINT16_MAX;
	constexpr auto scale = (double)UINT16_MAX / 2.0;
	uint64_t r = seed + pos[0];
	for (uint64_t n = 1; n < N; n++) {
		r = xorshift(r) + pos[n];
	}
	std::array<double, N> grad;
	grad[0] = (double)(r & mask) / scale - 1.0;
	for (uint64_t n = 1; n < N; n++) {
		grad[n] = (double)(xorshift(r + n) & mask) / scale - 1.0;
	}
	return std::move(grad);
}
template<uint64_t N>
double perlin(std::array<double, N> pos_decimal, uint64_t seed = 0) {
	std::array<int64_t, N> pos_int;
	for (uint64_t n = 0; n < N; n++) {
		// integer
		pos_int[n] = (int64_t)pos_decimal[n];
		// After decimal point
		pos_decimal[n] = pos_decimal[n] - pos_int[n];
	}

	// wave
	constexpr uint64_t wave_count = pow_int(2, N - 1);
	std::array<double, wave_count> wavs;
	for (uint64_t i = 0; i < wave_count; i++) {
		//initialize positions
		std::array<double, N> p_d = pos_decimal;
		std::array<int64_t, N> p_i = pos_int;
		uint64_t num = 2*i;
		for (uint64_t n = 1; n < N; n++) {
			uint64_t b = ((num & ((uint64_t)1 << n)) != 0); // does i stand n bit
			p_d[n] -= (double)b;
			p_i[n] += b;
		}
		//do wave
		double wav = wave<N>(p_d, pos2grad<N>(p_i, seed));
		//initialize positions
		p_d = pos_decimal;
		p_i = pos_int;
		p_d[0] -= 1;
		p_i[0] += 1;
		for (uint64_t n = 1; n < N; n++) {
			uint64_t b = ((num & ((uint64_t)1 << n)) != 0); // does i stand n bit
			p_d[n] -= (double)b;
			p_i[n] += b;
		}
		//do wave
		wavs[i] = lerp(pos_decimal[0],
			wav,
			wave<N>(p_d, pos2grad<N>(p_i, seed))
		);
	}
	// lerp
	for (uint64_t n = 1; n < N; n++) {
		const uint64_t loop_count = pow_int(2, N - n);
		for (uint64_t i = 0; i < loop_count; i+=2) {
			wavs[i] = lerp(pos_decimal[n], wavs[i], wavs[i + 1]);
		}
	}
	return wavs[0];
}