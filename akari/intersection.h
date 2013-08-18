#ifndef _INTERSECTION_H_
#define _INTERSECTION_H_

#include "vec.h"
#include "constant.h"
#include "material.h"

struct Hitpoint {
	float distance_;
	Vec normal_;
	Vec position_;
	Color color_;
	
	int triangleID;

	Hitpoint() : distance_(kINF), normal_(), color_(), position_(), triangleID(-1) {}
};

struct Intersection {
	Hitpoint hitpoint_;

	int objectID;
	Material *objectMaterial;

	Intersection() : objectID(-1), objectMaterial(NULL) {}
};

#endif
