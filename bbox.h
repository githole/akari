#ifndef _BBOX_H_
#define _BBOX_H_

#include <algorithm>

#include "vec.h"
#include "constant.h"
#include "ray.h"

class BBox {
public:
	Vec pmin_, pmax_;

	BBox() {
		pmin_ = Vec(kINF, kINF, kINF);
		pmax_ = -1.0 * pmin_;
	}

	BBox(const Vec &p) : pmin_(p), pmax_(p) {}

	BBox(const Vec &p1, const Vec &p2) {
		pmin_ = Vec(std::min(p1.x_, p2.x_), std::min(p1.y_, p2.y_), std::min(p1.z_, p2.z_));
		pmax_ = Vec(std::max(p1.x_, p2.x_), std::max(p1.y_, p2.y_), std::max(p1.z_, p2.z_));
	}

	bool inside(const Vec &pt) const {
		return (pmin_.x_ <= pt.x_ && pt.x_ <= pmax_.x_ &&
				pmin_.y_ <= pt.y_ && pt.y_ <= pmax_.y_ &&
				pmin_.z_ <= pt.z_ && pt.z_ <= pmax_.z_);
	}

	Vec &operator[](int i) {
		if (i == 0)
			return pmin_;
		return pmax_;
	}	
	const Vec &operator[](int i) const {
		if (i == 0)
			return pmin_;
		return pmax_;
	}

	void expand(float delta) {
		const Vec v(delta, delta, delta);
		pmin_ = pmin_ - v;
		pmax_ = pmax_ + v;
	}

	float surface_area() {
		const Vec d = pmax_ - pmin_;
		return 2.0f * (d.x_ * d.y_ + d.x_ * d.z_ + d.y_ * d.z_);
	}
	float volume() {
		const Vec d = pmax_ - pmin_;
		return d.x_ * d.y_ * d.z_;
	}

	enum LongestAxis {
		AxisX,
		AxisY,
		AxisZ,
	};

	LongestAxis maximum_extent() const {
		const Vec diag = pmax_ - pmin_;
		if (diag.x_ > diag.y_ && diag.x_ > diag.z_)
			return AxisX;
		else if (diag.y_ > diag.z_)
			return AxisY;
		else
			return AxisZ;
	}

	bool check_intersect(const Ray &ray, float *hitt0, float *hitt1) {
		float t0 = 0.0, t1 = kINF;
		for (int i = 0; i < 3; ++i) {
			// Update interval for _i_th bounding box slab
			float invRayDir = 1.f / ray.dir_[i];
			float tNear = (pmin_[i] - ray.org_[i]) * invRayDir;
			float tFar  = (pmax_[i] - ray.org_[i]) * invRayDir;

			// Update parametric interval from slab intersection $t$s
			if (tNear > tFar) std::swap(tNear, tFar);
			t0 = tNear > t0 ? tNear : t0;
			t1 = tFar  < t1 ? tFar  : t1;
			if (t0 > t1) return false;
		}
		if (hitt0) *hitt0 = t0;
		if (hitt1) *hitt1 = t1;
		return true;
	}
};


inline BBox Union(const BBox &b, const Vec &p) {
	BBox ret = b;
	ret.pmin_.x_ = std::min(b.pmin_.x_, p.x_);
	ret.pmin_.y_ = std::min(b.pmin_.y_, p.y_);
	ret.pmin_.z_ = std::min(b.pmin_.z_, p.z_);
		
	ret.pmax_.x_ = std::max(b.pmax_.x_, p.x_);
	ret.pmax_.y_ = std::max(b.pmax_.y_, p.y_);
	ret.pmax_.z_ = std::max(b.pmax_.z_, p.z_);

	return ret;
}

inline BBox Union(const BBox &b1, const BBox &b2) {
	BBox ret = b1;
	ret.pmin_.x_ = std::min(b1.pmin_.x_, b2.pmin_.x_);
	ret.pmin_.y_ = std::min(b1.pmin_.y_, b2.pmin_.y_);
	ret.pmin_.z_ = std::min(b1.pmin_.z_, b2.pmin_.z_);
	
	ret.pmax_.x_ = std::max(b1.pmax_.x_, b2.pmax_.x_);
	ret.pmax_.y_ = std::max(b1.pmax_.y_, b2.pmax_.y_);
	ret.pmax_.z_ = std::max(b1.pmax_.z_, b2.pmax_.z_);

	return ret;
}

#endif