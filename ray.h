#ifndef _RAY_H_
#define _RAY_H_

#include "vec.h"

struct Ray {
	Vec org_, dir_;
	Ray(const Vec &org, const Vec &dir) : org_(org), dir_(dir) {}
	Ray() {}
};


#endif
