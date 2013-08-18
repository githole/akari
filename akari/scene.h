#ifndef _SCENE_H_
#define _SCENE_H_

#include <iostream>
#include <string>
#include <fstream>
#include <map>

#include "vec.h"
#include "ray.h"
#include "material.h"
#include "triangle_mesh_vvv.h"

class Scene {
private:
	struct Object {
		VVVTriangleMesh *mesh_;
		Material *material_;
	};

	std::vector<Object> objects_;
public:

	Scene(const std::string &filename) {
		std::ifstream ifs(filename);

		if (!ifs)
			return;

		while (!ifs.eof()) {
			std::string filename;
			ifs >> filename;

			float scale = 0.0f;
			ifs >> scale;

			Vec pos;
			ifs >> pos.x_ >> pos.y_ >> pos.z_;

			std::string material_type;
			ifs >> material_type;

			Material *material;

			if (material_type == "diffuse") {
				material = new MaterialDiffuse;
			} else if (material_type == "mirror") {
				material = new MaterialMirror;
			} else if (material_type == "glossy") {
				MaterialGlossy *p = new MaterialGlossy;
				material = p;

				ifs >> p->lobe_;
			} else if (material_type == "glass") {
				MaterialGlass *p = new MaterialGlass;
				material = p;

				ifs >> p->ior_;
				ifs >> p->surface_color_.x_  >> p->surface_color_.y_ >> p->surface_color_.z_;
				ifs >> p->sigma_t_.x_  >> p->sigma_t_.y_ >> p->sigma_t_.z_;
			}

			VVVTriangleMesh *mesh = new VVVTriangleMesh(filename, scale, pos);

			Object o;
			o.mesh_ = mesh;
			o.material_ = material;

			objects_.push_back(o);
		}
	}

	bool intersect(const Ray &ray, const Intersection &old_intersection, Intersection *intersection) {
		intersection->objectID = -1;

		for (int i = 0; i < objects_.size(); ++i) {
			Hitpoint hitpoint;
			if (i == old_intersection.objectID) {
				hitpoint = old_intersection.hitpoint_;
			} else {
				hitpoint = Hitpoint();
			}

			if (objects_[i].mesh_->intersect(ray, hitpoint, &intersection->hitpoint_, objects_[i].material_->type_ == MT_Glass)) {
				intersection->objectID = i;
				intersection->objectMaterial = objects_[i].material_;
			}
		}

		if (intersection->objectID == -1)
			return false;

		return true;
	}

	virtual ~Scene() {
		for (int i = 0; i < objects_.size(); ++i) {
			delete objects_[i].mesh_;
			delete objects_[i].material_;
		}
	}
};

#endif