#ifndef _QBVH_
#define _QBVH_

#include <vector>
#include <algorithm>
#include <queue>

/*
#ifdef _MSC_VER
#ifndef _M_IX86
#error "No SSE!"
#endif
#else
#ifndef __SSE__
#error "No SSE!"
#endif
#endif
*/

#include <xmmintrin.h>

#include "vec.h"
#include "bbox.h"
#include "triangle.h"

static const int QOrderTable[] = {
	0x44444,0x44444,0x44444,0x44444,0x44444,0x44444,0x44444,0x44444,
	0x44440,0x44440,0x44440,0x44440,0x44440,0x44440,0x44440,0x44440,
	0x44441,0x44441,0x44441,0x44441,0x44441,0x44441,0x44441,0x44441,
	0x44401,0x44401,0x44410,0x44410,0x44401,0x44401,0x44410,0x44410,
	0x44442,0x44442,0x44442,0x44442,0x44442,0x44442,0x44442,0x44442,
	0x44402,0x44402,0x44402,0x44402,0x44420,0x44420,0x44420,0x44420,
	0x44412,0x44412,0x44412,0x44412,0x44421,0x44421,0x44421,0x44421,
	0x44012,0x44012,0x44102,0x44102,0x44201,0x44201,0x44210,0x44210,
	0x44443,0x44443,0x44443,0x44443,0x44443,0x44443,0x44443,0x44443,
	0x44403,0x44403,0x44403,0x44403,0x44430,0x44430,0x44430,0x44430,
	0x44413,0x44413,0x44413,0x44413,0x44431,0x44431,0x44431,0x44431,
	0x44013,0x44013,0x44103,0x44103,0x44301,0x44301,0x44310,0x44310,
	0x44423,0x44432,0x44423,0x44432,0x44423,0x44432,0x44423,0x44432,
	0x44023,0x44032,0x44023,0x44032,0x44230,0x44320,0x44230,0x44320,
	0x44123,0x44132,0x44123,0x44132,0x44231,0x44321,0x44231,0x44321,
	0x40123,0x40132,0x41023,0x41032,0x42301,0x43201,0x42310,0x43210,
};

template <class Triangle>
class QBVHSSE {
private:
	int maxPrimsInNode_;

	
	struct BVHPrimitiveInfo {
		int primitiveNumber_;
		Vec centroid_;
		BBox bounds_;

		BVHPrimitiveInfo(const int pn, const BBox& b) :
			primitiveNumber_(pn), bounds_(b){
				centroid_ = 0.5f * b.pmin_ + 0.5f * b.pmax_;
		}
	};


	// 四つの三角形をパックする
	struct SIMDTrianglePack {
		__m128 x_[3];
		__m128 y_[3];
		__m128 z_[3];
		int idx_[4];
	};

	struct BVHBuildNode  {
		BBox bounds_;
		BVHBuildNode* children_[2];
		int splitAxis_, firstPrimOffset_, nPrimitives_;

		int simdTrisIdx_;

		BVHBuildNode() {
		}

		void InitLeaf(int first, int n, const BBox& b, const int simdTrisIdx) {
			firstPrimOffset_ = first;
			nPrimitives_ = n;
			bounds_ = b;
			simdTrisIdx_ = simdTrisIdx;
			splitAxis_ = 0;
		}

		void InitInterior(int axis, BVHBuildNode* c0, BVHBuildNode* c1) {
			children_[0] = c0;
			children_[1] = c1;
			bounds_ = Union(c0->bounds_, c1->bounds_);
			splitAxis_ = axis;
			firstPrimOffset_ = -1;
			nPrimitives_ = 0;
		}
	};

	struct ComparePoints {
		int dim_;
		ComparePoints(const int d) : dim_(d){}
		bool operator()(const BVHPrimitiveInfo &a, const BVHPrimitiveInfo &b) const {
			return a.centroid_[dim_] < b.centroid_[dim_];
		}
	};

	struct CompareToBucket {
		int splitBucket_, nBuckets_, dim_;
		const BBox &centroidBounds_;

		CompareToBucket(int split, int num, int d, const BBox &b)
			: centroidBounds_(b) {
				splitBucket_ = split;
				nBuckets_ = num;
				dim_ = d;
		}

		bool operator()(const BVHPrimitiveInfo &p) const {
			int b = (int)(nBuckets_ * ((p.centroid_[dim_] - centroidBounds_.pmin_[dim_]) / 
								(centroidBounds_.pmax_[dim_] - centroidBounds_.pmin_[dim_])));
			if (b == nBuckets_) 
				b = nBuckets_ - 1;
			return b <= splitBucket_;
		};
	};

	struct Children {
		union {
			struct Node {
				unsigned flag_ : 1;
				unsigned index_: 31;
			} node_;
			struct Leaf {
				unsigned flag_ : 1;
				unsigned nPrimitives_: 3;
				unsigned index_: 28;
			} leaf_;

			unsigned int raw_;
		};
	};

	struct SIMDBVHNode{
		__m128 bboxes[2][3];
		Children children[4];
		int axis_top;
		int axis_left;
		int axis_right;
		int reserved;
	};

public:
	std::vector<Triangle*> orderedTris_;
	std::vector<Triangle*> tris_;
	std::vector<SIMDTrianglePack*> simdTris_;
	std::vector<SIMDBVHNode*> simdNodes_;;
	
	int stats[16];

	__m128 zero_;
	__m128 one_;
	__m128 inf_;
	__m128 keps_;

	QBVHSSE() {
		for (int i = 0; i < 16; ++i)
			stats[i] = 0;
		std::cout << sizeof(SIMDTrianglePack) << std::endl;
		__declspec(align(16)) float one_f[4] = {1.0f, 1.0f, 1.0f, 1.0f};
		__declspec(align(16)) float inf_f[4] = {kINF, kINF, kINF, kINF};
		__declspec(align(16)) float zero_f[4] = {0.0f, 0.0f, 0.0f, 0.0f};
		__declspec(align(16)) float keps_f[4] = {kEPS, kEPS, kEPS, kEPS};
		
		zero_ = _mm_load_ps(zero_f);
		one_ = _mm_load_ps(one_f);
		inf_ = _mm_load_ps(inf_f);
		keps_ = _mm_load_ps(keps_f);
	};

	BVHBuildNode* recursiveBuild(std::vector<BVHPrimitiveInfo> &buildData, int start, int end, int *totalNodes, std::vector<Triangle*> &orderedPrims) {
		(*totalNodes) ++;
		BVHBuildNode* node = new BVHBuildNode;

		BBox bbox;
		for (int i = start; i < end; ++i)
			bbox = Union(bbox, buildData[i].bounds_);

		int nPrimitives = end - start;
		if (nPrimitives <= 4) {
			// 葉
			const int firstPrimOffset = (int)orderedPrims.size();

			SIMDTrianglePack *simdt = (SIMDTrianglePack*)_aligned_malloc(sizeof(SIMDTrianglePack), 16);
			
			__declspec(align(16)) float x[4 * 3] = {0};
			__declspec(align(16)) float y[4 * 3] = {0};
			__declspec(align(16)) float z[4 * 3] = {0};

			int cnt = 0;
			for (int i = start; i < end; ++i , ++cnt){
				const int idx = buildData[i].primitiveNumber_;
				orderedPrims.push_back(tris_[idx]);
				
				const int t = cnt % 4;

				simdt->idx_[t] = firstPrimOffset + cnt;
				x[t] = tris_[idx]->p_[0]->x_;
				x[4 + t] = tris_[idx]->p_[1]->x_;
				x[8 + t] = tris_[idx]->p_[2]->x_;
				y[t] = tris_[idx]->p_[0]->y_;
				y[4 + t] = tris_[idx]->p_[1]->y_;
				y[8 + t] = tris_[idx]->p_[2]->y_;
				z[t] = tris_[idx]->p_[0]->z_;
				z[4 + t] = tris_[idx]->p_[1]->z_;
				z[8 + t] = tris_[idx]->p_[2]->z_;
			}
			for (; cnt < 4; ++cnt) {
				simdt->idx_[cnt%4] = -1;
			}

			for (int i = 0; i < 3; ++i) {
				simdt->x_[i] = _mm_load_ps(x + 4 * i);
				simdt->y_[i] = _mm_load_ps(y + 4 * i);
				simdt->z_[i] = _mm_load_ps(z + 4 * i);
			}

			simdTris_.push_back(simdt);
			node->InitLeaf(firstPrimOffset, nPrimitives, bbox, (int)simdTris_.size() - 1);
		} else {
			// 中間ノード
			BBox centroidBounds;
			for (int i = start; i < end; ++i)
				centroidBounds = Union(centroidBounds, buildData[i].centroid_);
			int dim = centroidBounds.maximum_extent();
			int mid = (start + end) / 2;

			if (nPrimitives <= 16) {
				std::nth_element(&buildData[start], &buildData[mid], &buildData[end-1] + 1, ComparePoints(dim));
			} else {
				// SAH
				const int nBuckets = 12;

				struct BucketInfo {
					BucketInfo() { count = 0; }
					int count;
					BBox bounds;
				};
				BucketInfo buckets[nBuckets];
				for (int i = start; i < end;  ++i ) {
					int b = (int)(nBuckets * ((buildData[i].centroid_[dim] - centroidBounds.pmin_[dim]) /
										(centroidBounds.pmax_[dim] - centroidBounds.pmin_[dim])));
					if (b == nBuckets) b = nBuckets - 1;
					buckets[b].count ++;
					buckets[b].bounds = Union(buckets[b].bounds, buildData[i].bounds_);
				}
				float cost[nBuckets - 1];
				for (int i = 0; i < nBuckets - 1; i ++) {
					BBox b0, b1;
					int count0 = 0, count1  = 0;
					for (int j = 0; j <= i; j ++) {
						b0 = Union(b0, buckets[j].bounds);
						count0 += buckets[j].count;
					}
					for (int j = i + 1; j < nBuckets; j ++) {
						b1 = Union(b1, buckets[j].bounds);
						count1 += buckets[j].count;
					}
					cost[i] += 0.125f + (count0 * b0.surface_area() + count1 * b1.surface_area()) / bbox.surface_area();
				}

				float minCost = cost[0];
				int minCostSplit = 0;
				for (int i = 1; i < nBuckets - 1; ++ i) {
					if (cost[i] < minCost) {
						minCost = cost[i];
						minCostSplit = i;
					}
				}

				if (nPrimitives > 0 || minCost < nPrimitives) {
					BVHPrimitiveInfo* pmid = std::partition(&buildData[start], &buildData[end-1] + 1, CompareToBucket(minCostSplit, nBuckets, dim, centroidBounds));
					mid = (int)(pmid - &buildData[0]);
				}
			}

			node->InitInterior(dim,
				recursiveBuild(buildData, start, mid, totalNodes, orderedPrims),
				recursiveBuild(buildData, mid, end, totalNodes, orderedPrims));
		}

		return node;
	}
	void collapse2QBVH(BVHBuildNode* node) {
		BVHBuildNode *lc = node->children_[0];
		BVHBuildNode *rc = node->children_[1];

		BVHBuildNode *c[4] = {0};
		
		SIMDBVHNode *n;
		n = (SIMDBVHNode*)_aligned_malloc(sizeof(SIMDBVHNode), 16);
		simdNodes_.push_back(n);
		n->axis_top = node->splitAxis_;
		n->axis_left = n->axis_right = 0;

		if (lc != NULL) {
			n->axis_left = lc->splitAxis_;
			if (lc->nPrimitives_ == 0) {
				c[0] = lc->children_[0];
				c[1] = lc->children_[1];
			} else {
				c[0] = lc;
			}
		}
		if (rc != NULL) {
			n->axis_right = rc->splitAxis_;
			if (rc->nPrimitives_ == 0) {
				c[2] = rc->children_[0];
				c[3] = rc->children_[1];
			} else {
				c[2] = rc;
			}
		}
		__declspec(align(16)) float bboxes[2][3][4];
		for (int j = 0; j < 3; ++j) {
			for (int k = 0; k < 4; ++k) {
				if (c[k] != NULL) {
					bboxes[0][j][k] = c[k]->bounds_.pmin_[j];
					bboxes[1][j][k] = c[k]->bounds_.pmax_[j];
				}
			}
		}
		for (int m = 0; m < 2; m++) {
			for (int j = 0; j < 3; j++) {
				n->bboxes[m][j] = _mm_load_ps(bboxes[m][j]);
			}
		}

		for (int i = 0; i < 4; ++i) {
			if (c[i] == NULL) {
				n->children[i].leaf_.flag_ = 1;
				n->children[i].leaf_.nPrimitives_ = 0;
				n->children[i].leaf_.index_ = 0;
			} else {
				if (c[i]->nPrimitives_ == 0) {
					n->children[i].node_.flag_ = 0;
					n->children[i].node_.index_= simdNodes_.size();
					collapse2QBVH(c[i]);
				} else {
					n->children[i].leaf_.flag_ = 1;
					n->children[i].leaf_.nPrimitives_ = c[i]->nPrimitives_;
					n->children[i].leaf_.index_ = c[i]->simdTrisIdx_;
				}
			}
		}

		return;
	}

	BVHBuildNode *rootNode_;
	void CreateBVHFromTriangle2s(const std::vector<Triangle*>& tris) {
		tris_ = tris;
		orderedTris_.clear();
		maxPrimsInNode_ = 32;

		std::vector<BVHPrimitiveInfo> buildData;
		for (int i = 0; i < tris_.size(); i ++) {
			BBox b = tris_[i]->ObjectBound();

			buildData.push_back(BVHPrimitiveInfo(i, b));
		}
		int totalNodes = 0;
		orderedTris_.reserve(tris_.size());

		rootNode_ = recursiveBuild(buildData, 0, (int)tris.size(), &totalNodes, orderedTris_);

		// collapse
		collapse2QBVH(rootNode_);

		SIMDBVHNode *root = simdNodes_[0];
		std::cout << "QBVH: " << simdNodes_.size() << std::endl;
	}

	inline int test_AABB(
		const __m128 bboxes[2][3],

		const __m128 org[3],
		const __m128 idir[3],
		const int sign[3],
		__m128 tmin, __m128 tmax 
		)
	{
		// x coordinate
		tmin = _mm_max_ps(
			tmin,
			_mm_mul_ps(_mm_sub_ps(bboxes[sign[0]][0],org[0]), idir[0])
			);
		tmax = _mm_min_ps(
			tmax,
			_mm_mul_ps(_mm_sub_ps(bboxes[1 - sign[0]][0], org[0]), idir[0])
			);

		// y coordinate
		tmin = _mm_max_ps(
			tmin,
			_mm_mul_ps(_mm_sub_ps(bboxes[sign[1]][1],org[1]), idir[1])
			);
		tmax = _mm_min_ps(
			tmax,
			_mm_mul_ps(_mm_sub_ps(bboxes[1 - sign[1]][1], org[1]), idir[1])
			);

		// z coordinate
		tmin = _mm_max_ps(
			tmin,
			_mm_mul_ps(_mm_sub_ps(bboxes[sign[2]][2],org[2]), idir[2])
			);
		tmax = _mm_min_ps(
			tmax,
			_mm_mul_ps(_mm_sub_ps(bboxes[1 - sign[2]][2], org[2]), idir[2])
			);
		return _mm_movemask_ps(_mm_cmpge_ps(tmax, tmin));
	}
	
	bool intersect(const Ray &ray, const Hitpoint &old_hitpoint, Hitpoint* hitpoint, const bool using_interpolated_normal = false) {
		__m128 sseOrg[3];
		__m128 sseiDir[3];
		int sign[3];
		Vec idir(1.0f / ray.dir_.x_, 1.0f / ray.dir_.y_, 1.0f / ray.dir_.z_);
		__declspec(align(16)) float r_idir_x[4] = {idir.x_, idir.x_, idir.x_, idir.x_};
		__declspec(align(16)) float r_idir_y[4] = {idir.y_, idir.y_, idir.y_, idir.y_};
		__declspec(align(16)) float r_idir_z[4] = {idir.z_, idir.z_, idir.z_, idir.z_};
					
		__declspec(align(16)) float r_org_x[4] = {ray.org_.x_, ray.org_.x_, ray.org_.x_, ray.org_.x_};
		__declspec(align(16)) float r_org_y[4] = {ray.org_.y_, ray.org_.y_, ray.org_.y_, ray.org_.y_};
		__declspec(align(16)) float r_org_z[4] = {ray.org_.z_, ray.org_.z_, ray.org_.z_, ray.org_.z_};

		__declspec(align(16)) float r_dir_x[4] = {ray.dir_.x_, ray.dir_.x_, ray.dir_.x_, ray.dir_.x_};
		__declspec(align(16)) float r_dir_y[4] = {ray.dir_.y_, ray.dir_.y_, ray.dir_.y_, ray.dir_.y_};
		__declspec(align(16)) float r_dir_z[4] = {ray.dir_.z_, ray.dir_.z_, ray.dir_.z_, ray.dir_.z_};
					
		__m128 dir_x = _mm_load_ps(r_dir_x);
		__m128 dir_y = _mm_load_ps(r_dir_y);
		__m128 dir_z = _mm_load_ps(r_dir_z);

		sseOrg[0] = _mm_load_ps(r_org_x);
		sseOrg[1] = _mm_load_ps(r_org_y);
		sseOrg[2] = _mm_load_ps(r_org_z);

		sseiDir[0] = _mm_load_ps(r_idir_x);
		sseiDir[1] = _mm_load_ps(r_idir_y);
		sseiDir[2] = _mm_load_ps(r_idir_z);
		
		sign[0] = idir[0] < 0;
		sign[1] = idir[1] < 0;
		sign[2] = idir[2] < 0;

		const SIMDBVHNode* nodes = simdNodes_[0];

		Children nodeStack[40];
		int todoNode = 0;

		nodeStack[0].raw_ = 0;

		bool hit = false;
		int triangle_index = -1;

		int cnt = 0;

		float rb1, rb2;

		while (todoNode >= 0) {
			Children item = nodeStack[todoNode];
			todoNode--;//pop stack

			if(item.node_.flag_ == 0){
				const SIMDBVHNode& node = *(simdNodes_[item.node_.index_]);
				__declspec(align(16)) float now_distance_f[4] = {hitpoint->distance_, hitpoint->distance_, hitpoint->distance_, hitpoint->distance_};
				__m128 now_distance = _mm_load_ps(now_distance_f);
				const int HitMask = test_AABB(node.bboxes, sseOrg, sseiDir, sign, zero_, now_distance);

				if (HitMask) {
					const int nodeIdx = (sign[node.axis_top] << 2) | (sign[node.axis_left] << 1) | (sign[node.axis_right]);
					int bboxOrder = QOrderTable[HitMask * 8 + nodeIdx];
					

					for (int i = 0; i < 4; ++i) {
						if (bboxOrder & 0x4)
							break;
						++todoNode;
						nodeStack[todoNode] = node.children[bboxOrder & 0x3];
						bboxOrder >>= 4;
					}
				}

			} else {
				__declspec(align(16)) float t_f[4];
				__declspec(align(16)) float hit_b1[4], hit_b2[4];
				int nohitmask;
				SIMDTrianglePack *s = simdTris_[item.leaf_.index_];

				__m128 e1_x = _mm_sub_ps(s->x_[1], s->x_[0]);
				__m128 e1_y = _mm_sub_ps(s->y_[1], s->y_[0]);
				__m128 e1_z = _mm_sub_ps(s->z_[1], s->z_[0]);
					
				__m128 e2_x = _mm_sub_ps(s->x_[2], s->x_[0]);
				__m128 e2_y = _mm_sub_ps(s->y_[2], s->y_[0]);
				__m128 e2_z = _mm_sub_ps(s->z_[2], s->z_[0]);

				__m128 s1_x = _mm_sub_ps(_mm_mul_ps(dir_y, e2_z), _mm_mul_ps(dir_z, e2_y));
				__m128 s1_y = _mm_sub_ps(_mm_mul_ps(dir_z, e2_x), _mm_mul_ps(dir_x, e2_z));
				__m128 s1_z = _mm_sub_ps(_mm_mul_ps(dir_x, e2_y), _mm_mul_ps(dir_y, e2_x));

				__m128 divisor = _mm_add_ps(_mm_add_ps(_mm_mul_ps(s1_x, e1_x), _mm_mul_ps(s1_y, e1_y)), _mm_mul_ps(s1_z, e1_z));
				__m128 no_hit  = _mm_cmpeq_ps(divisor, zero_);	

				__m128 invDivisor = _mm_rcp_ps(divisor);

				__m128 d_x = _mm_sub_ps(sseOrg[0], s->x_[0]);
				__m128 d_y = _mm_sub_ps(sseOrg[1], s->y_[0]);
				__m128 d_z = _mm_sub_ps(sseOrg[2], s->z_[0]);

				__m128 b1 = _mm_mul_ps(_mm_add_ps(_mm_add_ps(_mm_mul_ps(d_x, s1_x), _mm_mul_ps(d_y, s1_y)), _mm_mul_ps(d_z, s1_z)),
										invDivisor);
				no_hit = _mm_or_ps(no_hit, _mm_or_ps(_mm_cmplt_ps(b1, zero_), _mm_cmpgt_ps(b1, one_)));
					
				__m128 s2_x = _mm_sub_ps(_mm_mul_ps(d_y, e1_z), _mm_mul_ps(d_z, e1_y));
				__m128 s2_y = _mm_sub_ps(_mm_mul_ps(d_z, e1_x), _mm_mul_ps(d_x, e1_z));
				__m128 s2_z = _mm_sub_ps(_mm_mul_ps(d_x, e1_y), _mm_mul_ps(d_y, e1_x));

				__m128 b2 = _mm_mul_ps(_mm_add_ps(_mm_add_ps(_mm_mul_ps(dir_x, s2_x), _mm_mul_ps(dir_y, s2_y)), _mm_mul_ps(dir_z, s2_z)),
										invDivisor);
				no_hit = _mm_or_ps(no_hit, _mm_or_ps(_mm_cmplt_ps(b2, zero_), _mm_cmpgt_ps(_mm_add_ps(b1, b2), one_)));

				__m128 t = _mm_mul_ps(_mm_add_ps(_mm_add_ps(_mm_mul_ps(e2_x, s2_x), _mm_mul_ps(e2_y, s2_y)), _mm_mul_ps(e2_z, s2_z)),
										invDivisor);
				
				no_hit = _mm_or_ps(no_hit, _mm_cmplt_ps(t, keps_));
				
				nohitmask = _mm_movemask_ps(no_hit);
				_mm_store_ps(t_f, t);
				_mm_store_ps(hit_b1, b1);
				_mm_store_ps(hit_b2, b2);
				
				for (int i = 0; i < 4; ++i) {
					if ((nohitmask & (1 << i)) == 0 && hitpoint->distance_ > t_f[i] && old_hitpoint.triangleID != s->idx_[i]) {
						hit = true;
						triangle_index = s->idx_[i];
						hitpoint->distance_ = t_f[i];

						rb1 = hit_b1[i];
						rb2 = hit_b2[i];
					}
				}
			}
		}
		if (hit) {
			Triangle *t = orderedTris_[triangle_index];
			hitpoint->position_ = ray.org_+ hitpoint->distance_ * ray.dir_;
			hitpoint->triangleID = triangle_index;

			if (using_interpolated_normal) {
				hitpoint->normal_ = normalize(rb2 / (1 - rb1) * ((1 - rb1) * t->v_[2]->normal_ + rb1 * t->v_[1]->normal_) + 
					(1 - rb1 - rb2) / (1 - rb1) * ((1 - rb1) * t->v_[0]->normal_ + rb1 * t->v_[1]->normal_)) ;
			} else {
				hitpoint->normal_   = normalize(cross((*t->p_[2]) - (*t->p_[0]), (*t->p_[1]) - (*t->p_[0])));
			}
			hitpoint->color_ = rb2 / (1 - rb1) * ((1 - rb1) * t->v_[2]->color_ + rb1 * t->v_[1]->color_) + 
				(1 - rb1 - rb2) / (1 - rb1) * ((1 - rb1) * t->v_[0]->color_ + rb1 * t->v_[1]->color_) ;
		}

		return hit;
	}

	virtual ~QBVHSSE() {
	}

};

#endif