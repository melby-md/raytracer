#include <stdlib.h>

#include "common.hpp"
#include "mathlib.hpp"
#include "bmp.hpp"
#include "bsdf.hpp"

enum class ObjectType {
	TRIANGLE,
	SPHERE
};

struct Light {
	Vec3 color;
	i32 object_idx;
};

struct Triangle {
	Vec3 v0, v1, v2;
	Vec3 n0, n1, n2;
};

struct Sphere {
	Vec3 center;
	float radius;
};

struct Object {
	ObjectType type;
	i32 material_idx;
	i32 light_idx;
	union {
		Triangle triangle;
		Sphere sphere;
	};
};

struct Hit {
	Vec3 point;
	Vec3 normal;
	i32 obj_idx;
};

struct Ray {
	Vec3 origin, direction;
};

constexpr i32 MAX_OBJECTS = 64;
constexpr i32 MAX_MATERIALS = 64;
constexpr i32 MAX_LIGHTS = 64;

struct World {
	int object_count;
	Object objects[MAX_OBJECTS];

	int light_count;
	Light lights[MAX_LIGHTS];

	int material_count;
	Material materials[MAX_MATERIALS];

	Vec3 sun_color;
};

Vec2 RandomDisk(u32 *rng_state) {
	float theta = 2.f * PI * Rand(rng_state);
	float r = sqrtf(Rand(rng_state));

	float x = r * cosf(theta);
	float y = r * sinf(theta);

	return Vec2{x, y};
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

bool NearestHit(World *world, Ray ray, Hit *hit)
{
	float min_distance = INFINITY;
	Vec3 normal, point;
	i32 obj_idx = -1;

	for (i32 i = 0; i < world->object_count; i++) {
		Object *obj = &world->objects[i];

		float t;
		switch (obj->type) {
		case ObjectType::SPHERE:
			t = HitSphere(&obj->sphere, ray);
			if (t < min_distance) {
				point = ray.origin + ray.direction * t;
				normal = (point - obj->sphere.center) / obj->sphere.radius;

				obj_idx = i;
				min_distance = t;
			}
			break;

		case ObjectType::TRIANGLE:
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

bool Occluded(World *world, Ray ray, float distance)
{
	for (i32 i = 0; i < world->object_count; i++) {
		Object *obj = &world->objects[i];

		float t;
		switch (obj->type) {
		case ObjectType::SPHERE:
			t = HitSphere(&obj->sphere, ray);
			if (t < distance)
				return true;
			break;

		case ObjectType::TRIANGLE:
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

Vec3 RayTrace(World *world, Ray ray, u32 *rng_state)
{
	bool sample_lights = world->light_count > 0;

	int max_bounces = 10;

	Vec3 color = {};
	Vec3 throughput = {1, 1, 1};
	float bsdf_pdf = 1;

	for (int i = 0; i < max_bounces; i++) {
		Hit hit;
		bool did_hit = NearestHit(world, ray, &hit);

		if (!did_hit) {
			color += throughput * world->sun_color;
			break;
		}

		bool facing_forward = true;
		if (Dot(ray.direction, hit.normal) > 0) {
			hit.normal = -hit.normal;
			facing_forward = false;
		}

		Object *obj = &world->objects[hit.obj_idx];
		Material mat = world->materials[obj->material_idx];
		Light *light = obj->light_idx >= 0 ? &world->lights[obj->light_idx] : nullptr;

		Mat3 global_basis = OrthoNormalBasis(hit.normal);
		Mat3 local_basis = Transpose(global_basis);

		Vec3 v = local_basis * -ray.direction;

		if (facing_forward && light != nullptr) {
			if (i == 0 || !sample_lights) {
				color += throughput * light->color;
			} else {
				float pmf = 1.f / (float)world->light_count;
				float light_pdf = pmf * TrianglePDF(&obj->triangle, ray.origin, hit.point, hit.normal); 
				float mis_weight = PowerHeuristic(bsdf_pdf, light_pdf);

				Assert(light_pdf > 0);

				color += throughput * light->color * mis_weight;
			}
		}

		if (sample_lights) {
			u32 rand_idx = Rand(rng_state, 0, world->light_count-1);
			Light *sampled_light = &world->lights[rand_idx];
			Triangle *tri = &world->objects[sampled_light->object_idx].triangle;

			Vec3 uvw = RandomTriangle(rng_state);
			Vec3 light_point = uvw.x * tri->v0 + uvw.y * tri->v1 + uvw.z * tri->v2;
			Vec3 light_normal = Normalize(uvw.x * tri->n0 + uvw.y * tri->n1 + uvw.z * tri->n2);

			Vec3 light_direction = light_point - hit.point;
			float light_distance = Length(light_direction);
			light_direction /= light_distance;
			Ray shadow_ray = {hit.point, light_direction};
			Vec3 l = local_basis * light_direction;

			if (Dot(light_direction, light_normal) < 0 &&
			   !Occluded(world, shadow_ray, light_distance - .0001f)) {
				float pmf = 1.f / (float)world->light_count;
				float light_pdf = pmf * TrianglePDF(tri, hit.point, light_point, light_normal); 
				float b_pdf = BSDFPDF(v, l, mat);
				float mis_weight = PowerHeuristic(light_pdf, b_pdf);

				Assert(light_pdf > 0);

				color += throughput * sampled_light->color * BSDF(v, l, mat) * mis_weight / light_pdf;
			}
		}

		Sample sample = SampleBSDF(v, mat, rng_state);

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

		Assert(fabsf(Length(ray.direction) - 1) < .0001f);
	}

	return color;
}

i32 AddMaterial(World *world, Vec3 color, float roughness, float ior, float metallic)
{
	if (world->material_count >= MAX_MATERIALS)
		Panic("Too much materials");

	world->materials[world->material_count] = {
		color,
		roughness,
		ior,
		metallic
	};

	return world->material_count++;
}

void AddTriangle(World *world, Vec3 v0, Vec3 v1, Vec3 v2, i32 material_idx, Vec3 light)
{
	if (world->object_count >= MAX_OBJECTS)
		Panic("Too much objects");

	Vec3 e0 = v1 - v0;
	Vec3 e1 = v2 - v0;
	Vec3 normal = Normalize(Cross(e0, e1));

	Triangle tri = {
		v0, v1, v2,
		normal, normal, normal
	};

	Object obj = {ObjectType::TRIANGLE, material_idx, -1, tri};

	i32 obj_idx = world->object_count++;

	if (light.r > 0 || light.g > 0 || light.b > 0) {
		i32 light_idx = world->light_count++;
		Light *l = &world->lights[light_idx];
		l->color = light;
		l->object_idx = obj_idx;
		obj.light_idx = light_idx;
	}

	world->objects[obj_idx] = obj;
}

void LoadCornellBox(World *world)
{
	world->sun_color = {};

	world->object_count = 0;
	world->light_count = 0;
	world->material_count = 0;

	i32 khaki = AddMaterial(world, {0.725f, 0.71f, 0.68f}, 1, 1.5f, 0);
	i32 red = AddMaterial(world, {0.63f, 0.065f, 0.05f}, 1, 1.5f, 0);
	i32 green = AddMaterial(world, {0.14f, 0.45f, 0.091f}, 1, 1.5f, 0);
	i32 light = AddMaterial(world, {0, 0, 0}, 1, 1.5f, 0);

	AddTriangle(world, {-0.99f, -1, 0}, {1.04f, 0.99f, -0}, {-0.99f, 1.01f, 0}, khaki, {});
	AddTriangle(world, {-0.99f, -1, 0}, {1.04f, -1, -0}, {1.04f, 0.99f, -0}, khaki, {});
	AddTriangle(world, {1.04f, 1.02f, 1.99f}, {-0.99f, -1, 1.99f}, {-0.99f, 1.02f, 1.99f}, khaki, {});
	AddTriangle(world, {1.04f, 1.02f, 1.99f}, {1.04f, -1, 1.99f}, {-0.99f, -1, 1.99f}, khaki, {});
	AddTriangle(world, {1.04f, 0.99f, -0}, {1.04f, -1, 1.99f}, {1.04f, 1.02f, 1.99f}, khaki, {});
	AddTriangle(world, {1.04f, 0.99f, -0}, {1.04f, -1, -0}, {1.04f, -1, 1.99f}, khaki, {});
	AddTriangle(world, {-0.99f, -1, 0}, {1.04f, -1, 1.99f}, {1.04f, -1, -0}, green, {});
	AddTriangle(world, {-0.99f, -1, 0}, {-0.99f, -1, 1.99f}, {1.04f, -1, 1.99f}, green, {});
	AddTriangle(world, {-0.99f, 1.01f, 0}, {1.04f, 1.02f, 1.99f}, {-0.99f, 1.02f, 1.99f}, red, {});
	AddTriangle(world, {-0.99f, 1.01f, 0}, {1.04f, 0.99f, -0}, {1.04f, 1.02f, 1.99f}, red, {});
	AddTriangle(world, {-0.75f, -0.53f, 0.6f}, {0, -0.13f, 0.6f}, {-0.57f, 0.05f, 0.6f}, khaki, {});
	AddTriangle(world, {-0.57f, 0.05f, 0.6f}, {-0, -0.13f, 0}, {-0.57f, 0.05f, 0}, khaki, {});
	AddTriangle(world, {-0.75f, -0.53f, 0.6f}, {-0.57f, 0.05f, 0}, {-0.75f, -0.53f, 0}, khaki, {});
	AddTriangle(world, {-0.17f, -0.7f, 0.6f}, {-0.75f, -0.53f, 0}, {-0.17f, -0.7f, 0}, khaki, {});
	AddTriangle(world, {0, -0.13f, 0.6f}, {-0.17f, -0.7f, 0}, {-0, -0.13f, 0}, khaki, {});
	AddTriangle(world, {-0.75f, -0.53f, 0.6f}, {-0.17f, -0.7f, 0.6f}, {0, -0.13f, 0.6f}, khaki, {});
	AddTriangle(world, {-0.57f, 0.05f, 0.6f}, {0, -0.13f, 0.6f}, {-0, -0.13f, 0}, khaki, {});
	AddTriangle(world, {-0.75f, -0.53f, 0.6f}, {-0.57f, 0.05f, 0.6f}, {-0.57f, 0.05f, 0}, khaki, {});
	AddTriangle(world, {-0.17f, -0.7f, 0.6f}, {-0.75f, -0.53f, 0.6f}, {-0.75f, -0.53f, 0}, khaki, {});
	AddTriangle(world, {0, -0.13f, 0.6f}, {-0.17f, -0.7f, 0.6f}, {-0.17f, -0.7f, 0}, khaki, {});
	AddTriangle(world, {0.09f, -0.04f, 1.2f}, {0.49f, 0.71f, 1.2f}, {-0.09f, 0.53f, 1.2f}, khaki, {});
	AddTriangle(world, {-0.09f, 0.53f, 1.2f}, {0.49f, 0.71f, -0}, {-0.09f, 0.53f, 0}, khaki, {});
	AddTriangle(world, {0.49f, 0.71f, 1.2f}, {0.67f, 0.14f, -0}, {0.49f, 0.71f, -0}, khaki, {});
	AddTriangle(world, {0.67f, 0.14f, 1.2f}, {0.09f, -0.04f, -0}, {0.67f, 0.14f, -0}, khaki, {});
	AddTriangle(world, {0.09f, -0.04f, 1.2f}, {-0.09f, 0.53f, 0}, {0.09f, -0.04f, -0}, khaki, {});
	AddTriangle(world, {0.09f, -0.04f, 1.2f}, {0.67f, 0.14f, 1.2f}, {0.49f, 0.71f, 1.2f}, khaki, {});
	AddTriangle(world, {-0.09f, 0.53f, 1.2f}, {0.49f, 0.71f, 1.2f}, {0.49f, 0.71f, -0}, khaki, {});
	AddTriangle(world, {0.49f, 0.71f, 1.2f}, {0.67f, 0.14f, 1.2f}, {0.67f, 0.14f, -0}, khaki, {});
	AddTriangle(world, {0.67f, 0.14f, 1.2f}, {0.09f, -0.04f, 1.2f}, {0.09f, -0.04f, -0}, khaki, {});
	AddTriangle(world, {0.09f, -0.04f, 1.2f}, {-0.09f, 0.53f, 1.2f}, {-0.09f, 0.53f, 0}, khaki, {});
	AddTriangle(world, {0.22f, 0.24f, 1.98f}, {-0.16f, -0.23f, 1.98f}, {-0.16f, 0.24f, 1.98f}, light, {17, 12, 4});
	AddTriangle(world, {0.22f, 0.24f, 1.98f}, {0.22f, -0.23f, 1.98f}, {-0.16f, -0.23f, 1.98f}, light, {17, 12, 4});
}

int main()
{
	static World world;

	LoadCornellBox(&world);

	int width = 400, height = 400;
	int n_samples = 20;
	float exposure = 1;
	float fov = 90;
	Vec3 camera_position = {-1.9f, 0, 1};
	Vec3 looking_at = {0, 0, 1};
	Vec3 vup = {0, 0, 1};
	float defocus_angle = -2;
	float focus_dist = Length(looking_at - camera_position);

	float fov_radians = fov * PI / 180;
	float aspect_ratio = (float)width/(float)height;
	float viewport_height = 2 * tan(fov_radians/2) * focus_dist;
	float viewport_width = viewport_height * aspect_ratio;

	Vec3 forward = Normalize(looking_at - camera_position);
	Vec3 right = Normalize(Cross(forward, vup));
	Vec3 up = Cross(right, forward);

	Vec3 viewport_u = right * viewport_width;
	Vec3 viewport_v = -up * viewport_height;

	Vec3 upper_left = camera_position +
		forward * focus_dist
		- viewport_u/2
		- viewport_v/2;

	Vec3 du = viewport_u / (float)width;
	Vec3 dv = viewport_v / (float)height;

	float defocus_angle_radians = defocus_angle * PI / 180;
	float defocus_radius = focus_dist * tan(defocus_angle_radians / 2);
	Vec3 defocus_disk_u = right * defocus_radius;
	Vec3 defocus_disk_v = up * defocus_radius;

	byte *image_data = (byte *)malloc(sizeof(u8) * width * height * 3);

	double percent_row = 100 / (double)height;
	double percent_done = 0;

	Log("Raytracing... 0%%\r");

	#pragma omp parallel for
	for (int v = 0; v < height; v++) {
		u32 rng_state = 69420 + v;

		for (int u = 0; u < width; u++) {
			Vec3 color = {};
			for (int i = 0; i < n_samples; i++) {
				float random_u = Rand(&rng_state);
				float random_v = Rand(&rng_state);
				
				Vec2 rand_disk = RandomDisk(&rng_state);
				Vec3 rand_defocus = camera_position + (rand_disk.x * defocus_disk_u) + (rand_disk.y * defocus_disk_v);

				Vec3 ray_origin = (defocus_angle <= 0) ? camera_position : rand_defocus;
				Vec3 pixel_center = upper_left + du * ((float)u + random_u) + dv * ((float)v + random_v);
				Vec3 ray_direction = Normalize(pixel_center - ray_origin);
				Ray ray = {ray_origin, ray_direction};

				color += RayTrace(&world, ray, &rng_state);
			}

			color /= (float)n_samples;
			if (color.r < 0 || color.g < 0 || color.b < 0)
				color = {0, 0, 1};
			if (isnan(color.r) || isnan(color.g) || isnan(color.b))
				color = {0, 1, 0};
			Vec3 mapped = LinearToGamma(color, exposure);

			int pixel_pos = (v * width + u) * 3;

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
	putchar('\n');

	WriteBMP("image.bmp", width, height, image_data);

	return 0;
}
