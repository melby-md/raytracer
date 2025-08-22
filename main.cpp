#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stddef.h>

#include "common.hpp"
#include "mathlib.hpp"
#include "bmp.hpp"
#include "bsdf.hpp"
#include "scene.hpp"
#include "parser.hpp"

struct Hit {
	Vec3 point;
	Vec3 normal;
	i32 obj_idx;
};

struct Ray {
	Vec3 origin, direction;
};

Vec2 RandomDisk(u32 *rng_state) {
	for (;;) {
		float x = Rand(rng_state) * 2 - 1;
		float y = Rand(rng_state) * 2 - 1;

		if (x * y < 1)
			return Vec2{x, y};
	}
}

Vec3 RandomTriangle(u32 *rng_state)
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

Vec3 LinearToGamma(Vec3 color, float exposure)
{
	Vec3 mapped = 1. - Exp((-color) * exposure);

	return Pow(mapped, V3(1 / 2.2f));
}

float HitTriangle(Triangle *tri, Ray ray, float *_u, float *_v)
{
	Vec3 e0 = tri->v0 - tri->v2;
	Vec3 e1 = tri->v1 - tri->v2;
	Vec3 pvec = Cross(ray.direction, e1);
	float det = Dot(e0, pvec);

	if (det > -.0001f && det < .0001f)
		return INFINITY;

	Vec3 tvec = ray.origin - tri->v2;
	float u = Dot(tvec, pvec) / det;
	if (u < 0 || u > 1)
		return INFINITY;

	Vec3 qvec = Cross(tvec, e0);
	float v = Dot(ray.direction, qvec) / det;
	if (v < 0 || u + v > 1)
		return INFINITY;

	float t = Dot(e1, qvec) / det;

	*_u = u;
	*_v = v;

	if (t > .0001f)
		return t;

	return INFINITY;
}

float HitSphere(Sphere *sphere, Ray ray)
{
	Vec3 oc = sphere->center - ray.origin;
	float h = Dot(ray.direction, oc);
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

bool NearestHit(Scene *scene, Ray ray, Hit *hit)
{
	float min_distance = INFINITY;
	Vec3 normal, point;
	i32 obj_idx = -1;

	for (i32 i = 0; i < scene->object_count; i++) {
		Object *obj = &scene->objects[i];

		float t;
		switch (obj->type) {
		case OBJ_SPHERE:
			t = HitSphere(&obj->sphere, ray);
			if (t < min_distance) {
				point = ray.origin + ray.direction * t;
				normal = (point - obj->sphere.center) / obj->sphere.radius;

				obj_idx = i;
				min_distance = t;
			}
			break;

		case OBJ_TRIANGLE:
			float u, v;
			t = HitTriangle(&obj->triangle, ray, &u, &v);
			if (t < min_distance) {
				point = ray.origin + ray.direction * t;

				float w = 1 - u - v;
				normal = Normalize(obj->triangle.n0 * u + obj->triangle.n1 * v + obj->triangle.n2 * w);

				obj_idx = i;
				min_distance = t;
			}
		}
	}

	hit->point = point;
	hit->normal = normal;
	hit->obj_idx = obj_idx;

	return obj_idx >= 0;
}

bool Occluded(Scene *scene, Ray ray, float distance)
{
	for (i32 i = 0; i < scene->object_count; i++) {
		Object *obj = &scene->objects[i];

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

	return false;
}

float PowerHeuristic(float f_pdf, float g_pdf)
{
	return f_pdf*f_pdf / (f_pdf*f_pdf + g_pdf*g_pdf);
}

float TrianglePDF(Triangle *triangle, Vec3 point, Vec3 triangle_point, Vec3 triangle_normal)
{
	Vec3 e0 = triangle->v1 - triangle->v0;
	Vec3 e1 = triangle->v2 - triangle->v0;
	float area = Length(Cross(e0, e1)) / 2;

	Vec3 direction = Normalize(point - triangle_point);
	float length2 = Length2(point - triangle_point);
	return length2 / Dot(triangle_normal, direction) / area;
}

Vec3 RayTrace(Scene *scene, Ray ray, u32 *rng_state)
{
	bool sample_lights = scene->light_count > 0;

	int max_bounces = 10;

	Vec3 color = {};
	Vec3 throughput = {1, 1, 1};
	float bsdf_pdf = 1;

	for (int i = 0; i < max_bounces; i++) {
		Hit hit;
		bool did_hit = NearestHit(scene, ray, &hit);

		if (!did_hit) {
			color += throughput * scene->sun_color;
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
			   !Occluded(scene, shadow_ray, light_distance - .0001f)) {
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
			float roulette_prob = fmaxf(throughput.r, fmaxf(throughput.g, throughput.b));

			if (Rand(rng_state) > roulette_prob)
				break;

			throughput /= roulette_prob;
		}

		ray.origin = hit.point;
		ray.direction = global_basis * sample.l;
		bsdf_pdf = sample.pdf;
	}

	return color;
}

int main()
{
	static Scene scene;

	LoadScene(&scene, "scene.txt");

	float focus_dist = Length(scene.look_at - scene.camera);
	float fov_radians = scene.fov * PI / 180;
	float aspect_ratio = (float)scene.width/scene.height;
	float viewport_height = 2 * tan(fov_radians/2) * focus_dist;
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
	float defocus_radius = focus_dist * tan(defocus_angle_radians / 2);
	Vec3 defocus_disk_u = right * defocus_radius;
	Vec3 defocus_disk_v = up * defocus_radius;

	byte *image_data = (byte *)malloc(scene.width * scene.height * 3);

	double percent_row = 100. / scene.height;
	double percent_done = 0;

	scene.sun_color = V3(0);

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

				color += RayTrace(&scene, ray, &rng_state);
			}

			color /= scene.samples;
			if (color.r < 0 || color.g < 0 || color.b < 0)
				color = {0, 0, 1};
			if (isnan(color.r) || isnan(color.g) || isnan(color.b))
				color = {0, 1, 0};
			Vec3 mapped = LinearToGamma(color, scene.exposure);

			int pixel_pos = (v * scene.width + u) * 3;

			byte ir = (byte)(255 * Clamp(mapped.r, 0, 1));
			byte ig = (byte)(255 * Clamp(mapped.g, 0, 1));
			byte ib = (byte)(255 * Clamp(mapped.b, 0, 1));

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
