#include <errno.h>
#include <locale.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "common.h"

#include "bsdf.h"
#include "main.h"
#include "parser.h"

#include "bsdf.cpp"
#include "parser.cpp"

struct BVHNode {
	Vec3 aabb[2];
	i32 index, object_count;
};

struct Ray {
	Vec3 origin, direction;
};

struct Hit {
	Vec3 point;
	Vec3 normal;
	i32 obj_idx;
};


static float LinearToGamma(float color, float exposure)
{
	float m = 1 - expf(-color * exposure);

	m = Clamp(m, 0, 1);

	if (m <= .0031308f)
		return m * 12.92f;
	else
		return 1.055f * powf(m, 1 / 2.4f) - .055f;
}

static void WriteBMP(const char *file_name, int w, int h, byte *data)
{
	int filesize = 54 + 3*w*h;

	byte bmpfileheader[14] = {'B','M', 0,0,0,0, 0,0, 0,0, 54,0,0,0};
	byte bmpinfoheader[40] = {40,0,0,0, 0,0,0,0, 0,0,0,0, 1,0, 24,0};
	byte bmppad[3] = {0,0,0};

	bmpfileheader[ 2] = (byte)(filesize    );
	bmpfileheader[ 3] = (byte)(filesize>> 8);
	bmpfileheader[ 4] = (byte)(filesize>>16);
	bmpfileheader[ 5] = (byte)(filesize>>24);

	bmpinfoheader[ 4] = (byte)(       w    );
	bmpinfoheader[ 5] = (byte)(       w>> 8);
	bmpinfoheader[ 6] = (byte)(       w>>16);
	bmpinfoheader[ 7] = (byte)(       w>>24);
	bmpinfoheader[ 8] = (byte)(       h    );
	bmpinfoheader[ 9] = (byte)(       h>> 8);
	bmpinfoheader[10] = (byte)(       h>>16);
	bmpinfoheader[11] = (byte)(       h>>24);

	FILE *f = fopen(file_name, "wb");

	fwrite(bmpfileheader,1,14,f);
	fwrite(bmpinfoheader,1,40,f);
	for(int i=0; i<h; i++)
	{
		fwrite(data+(w*(h-i-1)*3),3,w,f);
		fwrite(bmppad,1,(4-(w*3)%4)%4,f);
	}

	fclose(f);
}

static u32 _rand(u32 *rng_state)
{
	// Xorshift32
	u32 x = *rng_state;
	x ^= x << 13;
	x ^= x >> 17;
	x ^= x << 5;
	*rng_state = x;

	return x;
}

float Rand(u32 *rng_state)
{
	return (float)(_rand(rng_state) >> 8) / (1<<24);
}

static u32 Rand(u32 *rng_state, u32 min, u32 max)
{
	return _rand(rng_state) % (max - min + 1) + min;
}

static Vec2 RandomDisk(u32 *rng_state) {
	for (;;) {
		float x = Rand(rng_state) * 2 - 1;
		float y = Rand(rng_state) * 2 - 1;

		if (x * y < 1)
			return {x, y};
	}
}

static Vec3 RandomTriangle(u32 *rng_state)
{
	float r1 = Rand(rng_state);
	float r2 = Rand(rng_state);

	float u, v;

	if (r1 < r2) {
		u = r1 / 2;
		v = r2 - u;
	} else {
		v = r2 / 2;
		u = r1 - v;
	}

	return {u, v, 1 - u - v};
}

static void UpdateBounds(BVHTree *tree, i32 index)
{
	BVHNode *node = &tree->nodes[index];
	Assert(node->object_count > 0);

	node->aabb[0] = V3(INFINITY);
	node->aabb[1] = V3(-INFINITY);

	i32 last = node->index + node->object_count;

	for (i32 i = node->index; i < last; i++) {
		Vec3 aabb_min = V3(INFINITY);
		Vec3 aabb_max = V3(-INFINITY);

		Object *obj = &tree->objects[i];

		switch (obj->type) {
		case OBJ_SPHERE:
			aabb_min = obj->sphere.center - V3(obj->sphere.radius);
			aabb_max = obj->sphere.center + V3(obj->sphere.radius);
			break;

		case OBJ_TRIANGLE:
			aabb_min = obj->triangle.v0;
			aabb_min = Min(aabb_min, obj->triangle.v1);
			aabb_min = Min(aabb_min, obj->triangle.v2);

			aabb_max = obj->triangle.v0;
			aabb_max = Max(aabb_max, obj->triangle.v1);
			aabb_max = Max(aabb_max, obj->triangle.v2);
		}

		node->aabb[0] = Min(node->aabb[0], aabb_min);
		node->aabb[1] = Max(node->aabb[1], aabb_max);
	}
}

static void BuildBVH(BVHTree *tree, i32 object_count)
{
	i32 node_count = 1;
	BVHNode *node = tree->nodes;
	node->index = 0;
	node->object_count = object_count;
	UpdateBounds(tree, 0);

	BVHNode *stack[64];
	int stack_ptr = 0;

	for (;;) {
		if (stack_ptr >= (int)CountOf(stack)) {
			node = stack[--stack_ptr];
			continue;
		}

		if (node->object_count <= 2) {
			if (stack_ptr == 0)
				break;
			node = stack[--stack_ptr];
			continue;
		}

		Vec3 extent = node->aabb[1] - node->aabb[0];
		int axis = 0;
		if (extent.y > extent.x)
			axis = 1;
		if (extent.z > extent[axis])
			axis = 2;
		float split_pos = node->aabb[0][axis] + extent[axis] * .5f;

		i32 i = node->index;
		i32 j = i + node->object_count - 1;
		while (i <= j) {
			Object *obj = &tree->objects[i];
			Vec3 centroid = {};
			switch (obj->type) {
			case OBJ_SPHERE:
				centroid = obj->sphere.center;
				break;
			case OBJ_TRIANGLE:
				centroid = (obj->triangle.v0 + obj->triangle.v1 + obj->triangle.v2) * (1.f / 3.f);
			}

			if (centroid[axis] < split_pos)
				i++;
			else {
				Object temp = tree->objects[i];
				tree->objects[i] = tree->objects[j];
				tree->objects[j] = temp;
				j--;
			}
		}

		i32 left_count = i - node->index;
		if (left_count == 0 || left_count == node->object_count) {
			if (stack_ptr == 0)
				break;
			node = stack[--stack_ptr];
			continue;
		}

		i32 left_child = node_count++;
		i32 right_child = node_count++;
		tree->nodes[left_child].index = node->index;
		tree->nodes[left_child].object_count = left_count;
		tree->nodes[right_child].index = i;
		tree->nodes[right_child].object_count = node->object_count - left_count;
		node->index = left_child;
		node->object_count = 0;

		UpdateBounds(tree, left_child);
		UpdateBounds(tree, right_child);

		stack[stack_ptr++] = &tree->nodes[right_child];
		node = &tree->nodes[left_child];
	}
}

static float HitTriangle(Triangle *tri, Ray *ray, float *_u, float *_v)
{
	Vec3 e0 = tri->v0 - tri->v2;
	Vec3 e1 = tri->v1 - tri->v2;
	Vec3 pvec = Cross(ray->direction, e1);
	float det = Dot(e0, pvec);

	if (det > -.0001f && det < .0001f)
		return INFINITY;

	Vec3 tvec = ray->origin - tri->v2;
	float u = Dot(tvec, pvec) / det;
	if (u < 0 || u > 1)
		return INFINITY;

	Vec3 qvec = Cross(tvec, e0);
	float v = Dot(ray->direction, qvec) / det;
	if (v < 0 || u + v > 1)
		return INFINITY;

	float t = Dot(e1, qvec) / det;

	*_u = u;
	*_v = v;

	if (t > .0001f)
		return t;

	return INFINITY;
}

static float HitSphere(Sphere *sphere, Ray *ray)
{
	Vec3 oc = sphere->center - ray->origin;
	float h = Dot(ray->direction, oc);
	float c = Length2(oc) - sphere->radius*sphere->radius;
	float delta = h*h - c;

	if (delta < .001f)
		return INFINITY;

	float sqd = sqrtf(delta);
	float distance = h - sqd;

	if (distance < .001f) {
		distance = h + sqd;
		if (distance < .001f)
			return INFINITY;
	}

	return distance;
}

static float IntersectAABB(Ray *ray, Vec3 *aabb, float max_distance)
{
	float tmin = 0.f;
	float tmax = INFINITY;

	bool sign = signbit(ray->direction.x);
	float tx1 = (aabb[sign].x - ray->origin.x) / ray->direction.x;
	float tx2 = (aabb[!sign].x - ray->origin.x) / ray->direction.x;

	tmin = Max(tx1, tmin);
	tmax = Min(tx2, tmax);

	sign = signbit(ray->direction.y);
	float ty1 = (aabb[sign].y - ray->origin.y) / ray->direction.y;
	float ty2 = (aabb[!sign].y - ray->origin.y) / ray->direction.y;

	tmin = Max(ty1, tmin);
	tmax = Min(ty2, tmax);

	sign = signbit(ray->direction.z);
	float tz1 = (aabb[sign].z - ray->origin.z) / ray->direction.z;
	float tz2 = (aabb[!sign].z - ray->origin.z) / ray->direction.z;

	tmin = Max(tz1, tmin);
	tmax = Min(tz2, tmax);

	return tmax >= tmin && tmin < max_distance && tmax > 0 ? tmin : INFINITY;
}

static bool NearestHit(BVHTree *tree, Ray *ray, Hit *hit)
{
	float min_distance = INFINITY;
	Vec3 normal, point;
	i32 obj_idx = -1;

	BVHNode *node = tree->nodes;
	BVHNode *stack[64];
	int stack_ptr = 0;

	for (;;) {
		if (node->object_count > 0) {
			i32 last = node->index + node->object_count;
			for (i32 i = node->index; i < last; i++) {
				Object *obj = &tree->objects[i];

				float t;
				switch (obj->type) {
				case OBJ_SPHERE:
					t = HitSphere(&obj->sphere, ray);
					if (t < min_distance) {
						point = ray->origin + ray->direction * t;
						normal = (point - obj->sphere.center) / obj->sphere.radius;

						obj_idx = i;
						min_distance = t;
					}
					break;

				case OBJ_TRIANGLE:
					float u, v;
					t = HitTriangle(&obj->triangle, ray, &u, &v);
					if (t < min_distance) {
						point = ray->origin + ray->direction * t;

						float w = 1 - u - v;
						normal = Normalize(obj->triangle.n0 * u + obj->triangle.n1 * v + obj->triangle.n2 * w);

						obj_idx = i;
						min_distance = t;
					}
				}
			}

			if (stack_ptr == 0)
				break;
			node = stack[--stack_ptr];
		} else {

			BVHNode *child1 = &tree->nodes[node->index];
			BVHNode *child2 = &tree->nodes[node->index + 1];

			float dist1 = IntersectAABB(ray, child1->aabb, min_distance);
			float dist2 = IntersectAABB(ray, child2->aabb, min_distance);

			BVHNode *near_child;
			BVHNode *far_child;
			float near_dist;
			float far_dist;
			if (dist1 < dist2) {
				near_child = child1;
				near_dist = dist1;
				far_child = child2;
				far_dist = dist2;
			} else {
				near_child = child2;
				near_dist = dist2;
				far_child = child1;
				far_dist = dist1;
			}

			if (near_dist < INFINITY) {
				node = near_child;
				if (far_dist < INFINITY)
					stack[stack_ptr++] = far_child;
			} else if (far_dist < INFINITY) {
				node = far_child;
			} else {
				if (stack_ptr == 0)
					break;
				node = stack[--stack_ptr];
			}
		}
	}

	hit->point = point;
	hit->normal = normal;
	hit->obj_idx = obj_idx;

	return obj_idx >= 0;
}

static bool Occluded(BVHTree *tree, Ray *ray, float distance)
{
	BVHNode *node = tree->nodes;
	BVHNode *stack[64];
	int stack_ptr = 0;

	for (;;) {
		if (node->object_count > 0) {
			i32 last = node->index + node->object_count;
			for (i32 i = node->index; i < last; i++) {
				Object *obj = &tree->objects[i];

				float t;
				switch (obj->type) {
				case OBJ_SPHERE:
					t = HitSphere(&obj->sphere, ray);
					if (t < distance)
						return true;
					break;

				case OBJ_TRIANGLE:
					float u, v;
					t = HitTriangle(&obj->triangle, ray, &u, &v);
					if (t < distance)
						return true;
				}
			}

			if (stack_ptr == 0)
				break;
			node = stack[--stack_ptr];
		} else {

			BVHNode *child1 = &tree->nodes[node->index];
			BVHNode *child2 = &tree->nodes[node->index + 1];

			bool isect1 = IntersectAABB(ray, child1->aabb, distance) < INFINITY;
			bool isect2 = IntersectAABB(ray, child2->aabb, distance) < INFINITY;

			if (isect1) {
				node = child1;
				if (isect2)
					stack[stack_ptr++] = child2;
			} else if (isect2) {
				node = child2;
			} else {
				if (stack_ptr == 0)
					break;
				node = stack[--stack_ptr];
			}
		}
	}

	return false;
}

static float PowerHeuristic(float f_pdf, float g_pdf)
{
	return f_pdf*f_pdf / (f_pdf*f_pdf + g_pdf*g_pdf);
}

static float TrianglePDF(Triangle *triangle, Vec3 point, Vec3 triangle_point, Vec3 triangle_normal)
{
	Vec3 e0 = triangle->v1 - triangle->v0;
	Vec3 e1 = triangle->v2 - triangle->v0;
	float area = Length(Cross(e0, e1)) / 2;

	Vec3 direction = Normalize(point - triangle_point);
	float length2 = Length2(point - triangle_point);
	return length2 / Dot(triangle_normal, direction) / area;
}

static Vec3 RayTrace(Scene *scene, Ray *_ray, u32 *rng_state)
{
	Ray ray = *_ray;
	bool sample_lights = scene->light_count > 0;

	int max_bounces = 10;

	Vec3 color = {};
	Vec3 throughput = {1, 1, 1};
	float bsdf_pdf = 1;

	for (int i = 0; i < max_bounces; i++) {
		Hit hit;
		bool did_hit = NearestHit(&scene->bvh, &ray, &hit);

		if (!did_hit) {
			color += throughput * scene->sky_box_color;
			break;
		}

		bool facing_forward = true;
		if (Dot(ray.direction, hit.normal) > 0) {
			hit.normal = -hit.normal;
			facing_forward = false;
		}

		Object *obj = &scene->objects[hit.obj_idx];
		Material mat = scene->materials[obj->material_idx];
		Light *light = obj->light_idx >= 0 ? &scene->lights[obj->light_idx] : nullptr;

		Mat3 global_basis = OrthoNormalBasis(hit.normal);
		Mat3 local_basis = Transpose(global_basis);

		Vec3 v = local_basis * -ray.direction;

		if (facing_forward && light != nullptr) {
			if (i == 0 || !sample_lights) {
				color += throughput * light->color;
			} else {
				float pmf = 1.f / (float)scene->light_count;
				float light_pdf = pmf * TrianglePDF(&obj->triangle, ray.origin, hit.point, hit.normal); 
				float mis_weight = PowerHeuristic(bsdf_pdf, light_pdf);

				Assert(light_pdf > 0);

				color += throughput * light->color * mis_weight;
			}
		}

		if (sample_lights) {
			u32 rand_idx = Rand(rng_state, 0, scene->light_count-1);
			Light *sampled_light = &scene->lights[rand_idx];
			Triangle *tri = &scene->objects[sampled_light->object_idx].triangle;

			Vec3 uvw = RandomTriangle(rng_state);
			Vec3 light_point = uvw.x * tri->v0 + uvw.y * tri->v1 + uvw.z * tri->v2;
			Vec3 light_normal = Normalize(uvw.x * tri->n0 + uvw.y * tri->n1 + uvw.z * tri->n2);

			Vec3 light_direction = light_point - hit.point;
			float light_distance = Length(light_direction);
			light_direction /= light_distance;
			Ray shadow_ray = {hit.point, light_direction};
			Vec3 l = local_basis * light_direction;

			if (Dot(light_direction, light_normal) < 0 &&
			   !Occluded(&scene->bvh, &shadow_ray, light_distance - .0001f)) {
				float pmf = 1.f / (float)scene->light_count;
				float light_pdf = pmf * TrianglePDF(tri, hit.point, light_point, light_normal); 
				float b_pdf = BSDFPDF(v, l, &mat);
				float mis_weight = PowerHeuristic(light_pdf, b_pdf);

				Assert(light_pdf > 0);

				color += throughput * sampled_light->color * BSDF(v, l, &mat) * mis_weight / light_pdf;
			}
		}

		Sample sample = SampleBSDF(v, &mat, rng_state);

		throughput *= sample.bsdf / sample.pdf;

		if (i > 3) {
			float roulette_prob = Max(throughput.r, Max(throughput.g, throughput.b));

			if (Rand(rng_state) < (1 - roulette_prob))
				break;

			throughput /= roulette_prob;
		}

		ray.origin = hit.point;
		ray.direction = global_basis * sample.l;
		bsdf_pdf = sample.pdf;
	}

	return color;
}

int main(int argc, char **argv)
{
	setlocale(LC_NUMERIC, "C");
	
	static Scene scene;

	if (argc < 2) {
		Log("Usage: %s [FILE]\n", argv[0]);
		return 1;
	}

	LoadScene(&scene, argv[1]);

	static BVHNode nodes[16384];
	BVHTree tree = {
		nodes,
		16384,
		scene.objects
	};

	BuildBVH(&tree, scene.object_count);
	scene.bvh = tree;

	for (i32 i = 0; i < scene.object_count; i++) {
		Object *obj = &scene.objects[i];
		if (obj->light_idx >= 0) {
			scene.lights[obj->light_idx].object_idx = i;
		}
	}

	float focus_dist = Length(scene.look_at - scene.camera);
	float fov_radians = scene.fov * PI / 180;
	float aspect_ratio = (float)scene.width/scene.height;
	float viewport_height = 2 * tanf(fov_radians/2) * focus_dist;
	float viewport_width = viewport_height * aspect_ratio;

	Vec3 forward = Normalize(scene.look_at - scene.camera);
	Vec3 right = Normalize(Cross(forward, scene.up));
	Vec3 up = Cross(right, forward);

	Vec3 viewport_u = right * viewport_width;
	Vec3 viewport_v = -up * viewport_height;

	Vec3 upper_left = scene.camera +
		forward * focus_dist
		- viewport_u/2
		- viewport_v/2;

	Vec3 du = viewport_u / scene.width;
	Vec3 dv = viewport_v / scene.height;

	float defocus_angle_radians = scene.defocus_angle * PI / 180;
	float defocus_radius = focus_dist * tanf(defocus_angle_radians / 2);
	Vec3 defocus_disk_u = right * defocus_radius;
	Vec3 defocus_disk_v = up * defocus_radius;

	byte *image_data = (byte *)malloc(scene.width * scene.height * 3);

	double percent_row = 100. / scene.height;
	double percent_done = 0;

	Log("Raytracing... 0%%\r");

	#pragma omp parallel for schedule(dynamic, 1)
	for (i16 v = 0; v < scene.height; v++) {
		u32 rng_state = 69420 + v;

		for (i16 u = 0; u < scene.width; u++) {
			Vec3 color = {};
			for (i16 i = 0; i < scene.samples; i++) {
				Vec3 ray_origin;
				if (scene.defocus_angle > 0) {

					Vec2 rand_disk = RandomDisk(&rng_state);
					ray_origin = scene.camera + (rand_disk.x * defocus_disk_u) + (rand_disk.y * defocus_disk_v);

				} else {
					ray_origin = scene.camera;
				}

				float u1;
				do {
					u1 = Rand(&rng_state);
				} while (u1 == 0);
				float u2 = Rand(&rng_state);

				float sigma = .5f;

				float jitter_u = sigma * sqrtf(-2 * logf(u1)) * cosf(2 * PI * u2);
				float jitter_v = sigma * sqrtf(-2 * logf(u1)) * sinf(2 * PI * u2);

				Vec3 pixel_center = upper_left + du * (u + .5f + jitter_u) + dv * (v + .5f + jitter_v);
				Vec3 ray_direction = Normalize(pixel_center - ray_origin);
				Ray ray = {ray_origin, ray_direction};

				color += RayTrace(&scene, &ray, &rng_state);
			}

			color /= scene.samples;
			if (color.r < 0 || color.g < 0 || color.b < 0)
				color = {0, 0, 1};
			if (isnan(color.r) || isnan(color.g) || isnan(color.b))
				color = {0, 1, 0};

			int pixel_pos = (v * scene.width + u) * 3;

			byte ir = (byte)(255 * LinearToGamma(color.r, scene.exposure));
			byte ig = (byte)(255 * LinearToGamma(color.g, scene.exposure));
			byte ib = (byte)(255 * LinearToGamma(color.b, scene.exposure));

			image_data[pixel_pos + 2] = ir;
			image_data[pixel_pos + 1] = ig;
			image_data[pixel_pos + 0] = ib;
		}

		#pragma omp critical
		{
			percent_done += percent_row;
			Log("Raytracing... %.0f%%\r", percent_done);
		}
	}
	Log("\n");

	WriteBMP("image.bmp", scene.width, scene.height, image_data);

	return 0;
}
