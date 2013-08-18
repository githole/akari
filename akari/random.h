#ifndef _RANDOM_H_
#define _RANDOM_H_

#include <climits>


// Xor-Shiftによる乱数ジェネレータ
class XorShift {
	unsigned int seed_[4];

	const float coefficient_;
public:
	unsigned int next(void) { 
		const unsigned int t = seed_[0] ^ (seed_[0] << 11);
		seed_[0] = seed_[1]; 
		seed_[1] = seed_[2];
		seed_[2] = seed_[3];
		return seed_[3] = (seed_[3] ^ (seed_[3] >> 19)) ^ (t ^ (t >> 8)); 
	}

	float next01(void) {
		return (float)next() * coefficient_;
	}

	XorShift(const unsigned int initial_seed) :
	coefficient_(1.0f / (1.0f + UINT_MAX)){
		unsigned int s = initial_seed;
		for (int i = 1; i <= 4; i++){
			seed_[i-1] = s = 1812433253U * (s^(s>>30)) + i;
		}
	}
};

typedef XorShift Random;

#endif
