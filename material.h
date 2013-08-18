#ifndef _MATERIAL_H_
#define _MATERIAL_H_

#include "color.h"

enum MaterialType {
	MT_Null,
	MT_Diffuse,
	MT_Glass,
	MT_Mirror,
	MT_Glossy
};

struct Material {
	enum MaterialType type_;

	Material() : type_(MT_Null) {}
};

struct MaterialDiffuse : public Material {
	MaterialDiffuse() {
		type_ = MT_Diffuse;
	}
};

struct MaterialMirror : public Material {
	MaterialMirror() {
		type_ = MT_Mirror;
	}
};

struct MaterialGlossy : public Material {
	float lobe_;

	MaterialGlossy() {
		type_ = MT_Glossy;
	}
};

struct MaterialGlass : public Material {
	float ior_;
	Color surface_color_;
	Color sigma_t_;
	
	MaterialGlass() {
		type_ = MT_Glass;
	}
};

#endif