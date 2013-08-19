#ifndef	_VEC_H_
#define	_VEC_H_

#include <cmath>

#include "constant.h"

struct Vec {
	float x_, y_, z_;
	Vec(const float x = 0, const float y = 0, const float z = 0) : x_(x), y_(y), z_(z) {}
	inline Vec operator+(const Vec &b) const {
		return Vec(x_ + b.x_, y_ + b.y_, z_ + b.z_);
	}
	inline Vec operator-(const Vec &b) const {
		return Vec(x_ - b.x_, y_ - b.y_, z_ - b.z_);
	}
	inline Vec operator*(const float b) const {
		return Vec(x_ * b, y_ * b, z_ * b);
	}
	inline Vec operator/(const float b) const {
		return Vec(x_ / b, y_ / b, z_ / b);
	}
	inline const float length_squared() const { 
		return x_*x_ + y_*y_ + z_*z_; 
	}
	inline const float length() const { 
		return sqrt(length_squared()); 
	}
	
	Vec operator-() const { return Vec(-x_, -y_, -z_); }
	Vec& operator+=(const Vec &v) {
		x_ += v.x_; y_ += v.y_; z_ += v.z_;
		return *this;
	}
	Vec& operator-=(const Vec &v) {
		x_ -= v.x_; y_ -= v.y_; z_ -= v.z_;
		return *this;
	}

	float operator[](const int i) const {
		return (&x_)[i];
	}
	float &operator[](const int i) {
		return (&x_)[i];
	}
};
inline Vec operator*(float f, const Vec &v) { 
	return v * f; 
}
inline Vec normalize(const Vec &v) {
	return v * (1.0f/ v.length()); 
}
inline const Vec multiply(const Vec &v1, const Vec &v2) {
	return Vec(v1.x_ * v2.x_, v1.y_ * v2.y_, v1.z_ * v2.z_);
}
inline const float dot(const Vec &v1, const Vec &v2) {
	return v1.x_ * v2.x_ + v1.y_ * v2.y_ + v1.z_ * v2.z_;
}
inline const Vec cross(const Vec &v1, const Vec &v2) {
	return Vec(
		(v1.y_ * v2.z_) - (v1.z_ * v2.y_),
		(v1.z_ * v2.x_) - (v1.x_ * v2.z_),
		(v1.x_ * v2.y_) - (v1.y_ * v2.x_));
}


void put(const Vec &v) {
	std::cout << "<" << v.x_ << "," << v.y_ << "," << v.z_ <<  ">";
}

bool valid(const Vec &v) {
	return (-kINF < v.x_ && v.x_ < kINF && !_isnan(v.x_)) &&
		(-kINF < v.y_ && v.y_ < kINF && !_isnan(v.y_)) &&
		(-kINF < v.z_ && v.z_ < kINF && !_isnan(v.z_));
}

#endif
