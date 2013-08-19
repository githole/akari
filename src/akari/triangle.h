#ifndef _TRIANGLE_H_
#define _TRIANGLE_H_

#include "vec.h"
#include "constant.h"
#include "bbox.h"
#include "intersection.h"

// どこか別のメモリ上に存在する頂点データに対する参照を保持する三角形データ構造
class RefTriangle {
private:
	float innerT;
public:
	Vec *p_[3];
	RefTriangle(Vec *p1, Vec *p2, Vec *p3) {
		p_[0] = p1;
		p_[1] = p2;
		p_[2] = p3;
	}

	BBox ObjectBound() const {
		return Union(BBox(*p_[0], *p_[1]), *p_[2]);
	}
	
	bool intersect(const Ray &ray, Hitpoint* hitpoint) {    // Get triangle vertices in _p1_, _p2_, and _p3_
		const Vec &p1 = *p_[0];
		const Vec &p2 = *p_[1];
		const Vec &p3 = *p_[2];
		const Vec e1 = p2 - p1;
		const Vec e2 = p3 - p1;
		const Vec s1 = cross(ray.dir_, e2);
		const float divisor = dot(s1, e1);
	
		if (divisor == 0.0f)
			return false;
		float invDivisor = 1.0f / divisor;

		// Compute first barycentric coordinate
		const Vec d = ray.org_ - p1;
		float b1 = dot(d, s1) * invDivisor;
		if (b1 < 0. || b1 > 1.)
			return false;

		// Compute second barycentric coordinate
		Vec s2 = cross(d, e1);
		float b2 = dot(ray.dir_, s2) * invDivisor;
		if (b2 < 0. || b1 + b2 > 1.)
			return false;

		// Compute _t_ to intersection point
		float t = dot(e2, s2) * invDivisor;
		if (t < kEPS || t > kINF)
			return false;

		hitpoint->normal_ = normalize(cross(e2, e1));
		hitpoint->distance_ = t;
		hitpoint->position_ = ray.org_ + hitpoint->distance_ * ray.dir_;
		return true;
	}
};


#endif