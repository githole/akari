#ifndef _RENDER_H_
#define _RENDER_H_

#include "scene.h"
#include "setting.h"
#include "vec.h"
#include "hdr.h"
#include "random.h"
#include "ibl.h"
#include "triangle_mesh_vvv.h"

class Render {
private:
	Scene *scene_;
	IBL *ibl_;
	Setting *setting_;

	int depth_max_;
public:

	Render(Setting *setting) : setting_(setting) {
		scene_ = new Scene(setting->string_value("scene", "scene.txt"));
		ibl_   = new IBL(setting->string_value("ibl", "ibl.hdr"), true); 
		depth_max_ = setting->int_value("depth_max_", 32);
	}
	
	virtual ~Render() {
		delete scene_;
		delete ibl_;
	}
	
	Color radiance_indirect_light_env(const Ray &ray, Random *rnd, const Intersection &old_intersection, const Material *around_material, const int depth, bool inside = false) {
		if (depth == depth_max_)
			return Color(0, 0, 0);

		// シーンと交差判定
		Intersection intersection;
		if (!scene_->intersect(ray, old_intersection, &intersection)) {
			// 直接光除外
			if (depth >= 2)
				return ibl_->sample(ray.dir_);
			else
				return Color();
		}

		const Hitpoint &hitpoint = intersection.hitpoint_;
		const Vec orienting_normal = dot(hitpoint.normal_ , ray.dir_) < 0.0 ? hitpoint.normal_: (-1.0 * hitpoint.normal_); // 交差位置の法線（物体からのレイの入出を考慮）

		Color incoming_radiance;
		Color weight;


		switch (intersection.objectMaterial->type_) {
		case MT_Diffuse:
			{
				Vec w, u, v;
				w = orienting_normal;
				if (fabs(w.x_) > kEPS) // ベクトルwと直交するベクトルを作る。w.xが0に近い場合とそうでない場合とで使うベクトルを変える。
					u = normalize(cross(Vec(0.0, 1.0, 0.0), w));
				else
					u = normalize(cross(Vec(1.0, 0.0, 0.0), w));
				v = cross(w, u);
				// コサイン項を使った重点的サンプリング
				const float r1 = 2 * kPI * rnd->next01();
				const float r2 = rnd->next01(), r2s = sqrt(r2);
				Vec dir = normalize((
					u * cos(r1) * r2s +
					v * sin(r1) * r2s +
					w * sqrt(1.0f - r2)));

				incoming_radiance = radiance_indirect_light_env(Ray(hitpoint.position_, dir), rnd, intersection, around_material, depth + 1, inside);
				weight = hitpoint.color_;
			}
			break;
		case MT_Mirror:
			{
				const Ray reflection_ray = Ray(hitpoint.position_, ray.dir_ - hitpoint.normal_ * 2.0 * dot(hitpoint.normal_, ray.dir_));
				incoming_radiance = radiance_indirect_light_env(reflection_ray, rnd, intersection, around_material, depth + 1, inside);
				weight = hitpoint.color_;
			}
			break;
		case MT_Glossy:
			{
				MaterialGlossy *material = (MaterialGlossy*)intersection.objectMaterial;
				const Ray reflection_ray = Ray(hitpoint.position_, ray.dir_ - hitpoint.normal_ * 2.0 * dot(hitpoint.normal_, ray.dir_));

				Vec w, u, v;
				Vec dir;
				w = reflection_ray.dir_;
				if (fabs(w.x_) > kEPS)
					u = normalize(cross(Vec(0.0, 1.0, 0.0), w));
				else
					u = normalize(cross(Vec(1.0, 0.0, 0.0), w));
				v = cross(w, u);
				do {
					const float r1 = 2 * kPI * rnd->next01();
					const float r2 = -material->lobe_ * 2.0f * rnd->next01() + 1.0f;
					const float r2s = sqrt(1.0 - r2 * r2);
					dir = normalize((
						u * cos(r1) * r2s +
						v * sin(r1) * r2s +
						w * r2));

					if (dot(orienting_normal, dir) >= 0.0)
						break;
				} while(1);

				incoming_radiance = radiance_indirect_light_env(Ray(hitpoint.position_, dir), rnd, intersection, around_material, depth + 1, inside);
				weight = hitpoint.color_ * dot(orienting_normal, dir);
			}
			break;
		case MT_Glass: 
			{
				MaterialGlass *material = (MaterialGlass*)intersection.objectMaterial;

				const Ray reflection_ray = Ray(hitpoint.position_, ray.dir_ - hitpoint.normal_ * 2.0 * dot(hitpoint.normal_, ray.dir_));
				const bool into = dot(hitpoint.normal_, orienting_normal) > 0.0; // レイがオブジェクトから出るのか、入るのか

				// Snellの法則
				const float nc = 1.0f; // 真空の屈折率
				const float nt = material->ior_; // オブジェクトの屈折率
				const float nnt = into ? nc / nt : nt / nc;
				const float ddn = dot(ray.dir_, orienting_normal);
				const float cos2t = 1.0f - nnt * nnt * (1.0f - ddn * ddn);

				if (cos2t < 0.0) { // 全反射
					incoming_radiance = radiance_indirect_light_env(reflection_ray, rnd, intersection, around_material, depth + 1, inside);
					weight = material->surface_color_;
					break;
				}
				// 屈折の方向
				const Ray refraction_ray = Ray(hitpoint.position_,
					normalize(ray.dir_ * nnt - hitpoint.normal_ * (into ? 1.0f : -1.0f) * (ddn * nnt + sqrt(cos2t))));

				// SchlickによるFresnelの反射係数の近似を使う
				const float a = nt - nc, b = nt + nc;
				const float R0 = (a * a) / (b * b);

				const float c = 1.0f - (into ? -ddn : dot(refraction_ray.dir_, -1.0f * orienting_normal));
				const float Re = R0 + (1.0f - R0) * pow(c, 5.0f); // 反射方向の光が反射してray.dirの方向に運ぶ割合。同時に屈折方向の光が反射する方向に運ぶ割合。
				const float nnt2 = pow(into ? nc / nt : nt / nc, 2.0f); // レイの運ぶ放射輝度は屈折率の異なる物体間を移動するとき、屈折率の比の二乗の分だけ変化する。
				const float Tr = (1.0f - Re) * nnt2; // 屈折方向の光が屈折してray.dirの方向に運ぶ割合

				// 一定以上レイを追跡したら屈折と反射のどちらか一方を追跡する。（さもないと指数的にレイが増える）
				// ロシアンルーレットで決定する。
				const float probability  = 0.25f + 0.5f * Re;
				if (rnd->next01() < probability) { // 反射
					incoming_radiance = radiance_indirect_light_env(reflection_ray, rnd, intersection, around_material, depth + 1, inside) * Re;
					weight = material->surface_color_ / probability;
				} else { // 屈折;
					Material *next_around_material = NULL;
					if (around_material == NULL) {
						next_around_material = intersection.objectMaterial;
					}
					inside = !inside;
					incoming_radiance = radiance_indirect_light_env(refraction_ray, rnd, intersection, next_around_material, depth + 1, inside) * Tr;
					weight = material->surface_color_ / (1.0f - probability);
				}
			}
			break;
		}

		if (around_material != NULL && around_material->type_ == MT_Glass) {
			MaterialGlass *m = (MaterialGlass*)around_material;
			const float length = (ray.org_ - hitpoint.position_).length();
			weight = multiply(weight, Vec(exp(-length * m->sigma_t_.x_), exp(-length * m->sigma_t_.y_), exp(-length * m->sigma_t_.z_)));
		}

		return multiply(weight, incoming_radiance);
	}

	Color radiance_direct_light_env(const Ray &ray, Random *rnd, const int usamples, const int vsamples, const Intersection &old_intersection, Intersection *intersection, const int depth) {
		if (depth >= 2)
			return Color();

		// シーンと交差判定
		if (!scene_->intersect(ray, old_intersection, intersection)) {
			return ibl_->sample(ray.dir_);
		}
	
		const Hitpoint &hitpoint = intersection->hitpoint_;
		const Vec orienting_normal = dot(hitpoint.normal_ , ray.dir_) < 0.0 ? hitpoint.normal_: (-1.0 * hitpoint.normal_); // 交差位置の法線（物体からのレイの入出を考慮）
	
		Color incoming_radiance;
		Color weight;
		
		if (depth >= 1)
			return Color();

		switch (intersection->objectMaterial->type_) {
		case MT_Diffuse:
			{
				std::vector<IBL::Sample> sample;
				sample.reserve(usamples * vsamples);

				ibl_->importance_sampling_unsafe_fast(usamples, vsamples, orienting_normal, sample, rnd);
				for (int i = 0; i < sample.size(); ++i) {
					Intersection next_intersection;
					if (!scene_->intersect(Ray(hitpoint.position_, sample[i].dir_), *intersection, &next_intersection)) {
						const float w = (1.0f / sample[i].weight_) * dot(sample[i].dir_, orienting_normal) / kPI;
						if (w > kINF || w < -kINF || w <= 0.0f) {
							continue;
						}
			
						Color col = ibl_->sample(sample[i].dir_);
						incoming_radiance = incoming_radiance + ibl_->sample(sample[i].dir_) * w / (float)sample.size();
					}
				}
				weight = hitpoint.color_;
			}
			break;
		case MT_Mirror:
			{
				Intersection next_intersection;
				const Ray reflection_ray = Ray(hitpoint.position_, ray.dir_ - hitpoint.normal_ * 2.0 * dot(hitpoint.normal_, ray.dir_));
				incoming_radiance = radiance_direct_light_env(reflection_ray, rnd, usamples, vsamples, *intersection, &next_intersection, depth + 1);
				weight = hitpoint.color_;
			}
			break;
		case MT_Glossy:
			{
				MaterialGlossy *material = (MaterialGlossy*)intersection->objectMaterial;
				const Ray reflection_ray = Ray(hitpoint.position_, ray.dir_ - hitpoint.normal_ * 2.0 * dot(hitpoint.normal_, ray.dir_));

				Vec w, u, v;
				Vec dir;
				w = reflection_ray.dir_;
				if (fabs(w.x_) > kEPS)
					u = normalize(cross(Vec(0.0, 1.0, 0.0), w));
				else
					u = normalize(cross(Vec(1.0, 0.0, 0.0), w));
				v = cross(w, u);
				for (int i = 0; i < usamples * vsamples; ++i) {
					do {
						const float r1 = 2 * kPI * rnd->next01();
						const float r2 = -material->lobe_ * 2.0f * rnd->next01() + 1.0f;
						const float r2s = sqrt(1.0 - r2 * r2);
						dir = normalize((
							u * cos(r1) * r2s +
							v * sin(r1) * r2s +
							w * r2));

						if (dot(orienting_normal, dir) >= 0.0)
							break;
					} while(1);

					Intersection next_intersection;
					incoming_radiance = incoming_radiance + radiance_direct_light_env(Ray(hitpoint.position_, dir), rnd, usamples, vsamples, *intersection, &next_intersection, depth + 1) / 
						(usamples * vsamples) * dot(orienting_normal, dir);
				}

				weight = hitpoint.color_;
			}
			break;
		case MT_Glass:
			{
				MaterialGlass *material = (MaterialGlass*)intersection->objectMaterial;
				weight = material->surface_color_;

				const Ray reflection_ray = Ray(hitpoint.position_, ray.dir_ - hitpoint.normal_ * 2.0 * dot(hitpoint.normal_, ray.dir_));
				const bool into = dot(hitpoint.normal_, orienting_normal) > 0.0; // レイがオブジェクトから出るのか、入るのか

				// Snellの法則
				const float nc = 1.0f; // 真空の屈折率
				const float nt = material->ior_; // オブジェクトの屈折率
				const float nnt = into ? nc / nt : nt / nc;
				const float ddn = dot(ray.dir_, orienting_normal);
				const float cos2t = 1.0f - nnt * nnt * (1.0f - ddn * ddn);

				if (cos2t < 0.0f) { // 全反射
					Intersection next_intersection;
					incoming_radiance = radiance_direct_light_env(reflection_ray, rnd, usamples, vsamples, *intersection, &next_intersection, depth + 1);
					break;
				}
				// 屈折の方向
				const Ray refraction_ray = Ray(hitpoint.position_,
					normalize(ray.dir_ * nnt - hitpoint.normal_ * (into ? 1.0f : -1.0f) * (ddn * nnt + sqrt(cos2t))));

				// SchlickによるFresnelの反射係数の近似を使う
				const float a = nt - nc, b = nt + nc;
				const float R0 = (a * a) / (b * b);

				const float c = 1.0f - (into ? -ddn : dot(refraction_ray.dir_, -1.0f * orienting_normal));
				const float Re = R0 + (1.0f - R0) * pow(c, 5.0f); // 反射方向の光が反射してray.dirの方向に運ぶ割合。同時に屈折方向の光が反射する方向に運ぶ割合。
				const float nnt2 = pow(into ? nc / nt : nt / nc, 2.0f); // レイの運ぶ放射輝度は屈折率の異なる物体間を移動するとき、屈折率の比の二乗の分だけ変化する。
				const float Tr = (1.0f - Re) * nnt2; // 屈折方向の光が屈折してray.dirの方向に運ぶ割合
			
				Intersection next_intersection0, next_intersection1;
				incoming_radiance = 
					radiance_direct_light_env(reflection_ray, rnd, usamples, vsamples, *intersection, &next_intersection0, depth + 1) * Re +
					radiance_direct_light_env(refraction_ray, rnd, usamples, vsamples, *intersection, &next_intersection1, depth + 1) * Tr;
			}
			break;
		}
	
		return multiply(weight, incoming_radiance);
	}


};

#endif