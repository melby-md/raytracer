#include <stdlib.h>

#include "common.hpp"
#include "mathlib.hpp"
#include "bmp.hpp"
#include "bsdf.hpp"

struct Sphere {
	int material_idx;

	Vec3 center;
	float radius;
};

struct Plane {
	int material_idx;
	Vec3 normal;
	float d;
};

struct Triangle {
	int material_idx;
	Vec3 v0, v1, v2;
	Vec3 n0, n1, n2;
};

struct HitInfo {
	Vec3 hit_point;
	Vec3 normal;
	int material_idx;
	int triangle_idx;
};

struct Ray {
	Vec3 origin, direction;
};

constexpr int MAX_OBJECTS = 64;
constexpr int MAX_MATERIALS = 64;

struct World {
	int sphere_count;
	Sphere spheres[MAX_OBJECTS];

	int triangle_count;
	Triangle triangles[MAX_OBJECTS];

	int emissive_triangle_count;
	Triangle emissive_triangles[MAX_OBJECTS];

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

	float g = 1 / 2.2f;
	return Pow(mapped, Vec3{g, g, g});
}

float HitTriangle(Vec3 v0, Vec3 v1, Vec3 v2, Ray ray, float *_u, float *_v)
{
	Vec3 e0 = v0 - v2;
	Vec3 e1 = v1 - v2;
	Vec3 pvec = Cross(ray.direction, e1);
	float det = Dot(e0, pvec);

	if (det > -.0001f && det < .0001f)
		return INFINITY;

	Vec3 tvec = ray.origin - v2;
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
	float a = Length2(ray.direction);
	float h = Dot(ray.direction, oc);
	float c = Length2(oc) - sphere->radius*sphere->radius;
	float delta = h*h - a*c;

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

bool NearestHit(World *world, Ray ray, HitInfo *info)
{
	bool hit = false;
	float min_distance = INFINITY;
	Vec3 normal, hit_point;
	int material_idx = -1;
	int triangle_idx = -1;

	for (int i = 0; i < world->sphere_count; i++) {
		Sphere *sphere = &world->spheres[i];

		float t = HitSphere(sphere, ray);

		if (t < min_distance) {
			hit_point = ray.origin + ray.direction * t;
			normal = (hit_point - sphere->center) / sphere->radius;
			material_idx = sphere->material_idx;
			min_distance = t;
			hit = true;
		}
	}

	for (int i = 0; i < world->triangle_count; i++) {
		Triangle *tri = &world->triangles[i];

		float u, v;
		float t = HitTriangle(tri->v0, tri->v1, tri->v2, ray, &u, &v);

		if (t < min_distance) {
			hit_point = ray.origin + ray.direction * t;
			material_idx = tri->material_idx;
			triangle_idx = i;

			float w = 1 - u - v;
			normal = Normalize(tri->n0 * u + tri->n1 * v + tri->n2 * w);

			min_distance = t;
			hit = true;
		}
	}

	info->hit_point = hit_point;
	info->normal = normal;
	info->material_idx = material_idx;
	info->triangle_idx = triangle_idx;

	return hit;
}

bool Occluded(World *world, Vec3 from, Vec3 to)
{
	Ray ray = {from, to - from};
	float distance = .99f;

	for (int i = 0; i < world->sphere_count; i++) {
		Sphere *sphere = &world->spheres[i];

		float t = HitSphere(sphere, ray);

		if (t < distance) {
			return true;
		}
	}

	for (int i = 0; i < world->triangle_count; i++) {
		Triangle *tri = &world->triangles[i];

		float u, v;
		float t = HitTriangle(tri->v0, tri->v1, tri->v2, ray, &u, &v);

		if (t < distance) {
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
	Vec3 direction = point - triangle_point;

	Vec3 e0 = triangle->v1 - triangle->v0;
	Vec3 e1 = triangle->v2 - triangle->v0;
	float area = Length(Cross(e0, e1)) / 2;

	return Length2(direction) / (Dot(triangle_normal, Normalize(direction))) / area;
}

Vec3 RayTrace(World *world, Ray ray, u32 *rng_state)
{
	bool sample_lights = world->emissive_triangle_count > 0;

	int max_bounces = 10;

	Vec3 color = {};
	Vec3 throughput = {1, 1, 1};
	float bsdf_pdf = 1;

	for (int i = 0; i < max_bounces; i++) {
		HitInfo hit;
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

		Material mat = world->materials[hit.material_idx];

		Mat3 global_basis = OrthoNormalBasis(hit.normal);
		Mat3 local_basis = Transpose(global_basis);

		Vec3 v = local_basis * -ray.direction;

		if (facing_forward && (mat.emission.x > 0 || mat.emission.y > 0 || mat.emission.z > 0)) {
			if (i == 0 || !sample_lights || hit.triangle_idx < 0) {
				color += throughput * mat.emission;
			} else {
				float pmf = 1.f / (float)world->emissive_triangle_count;
				float light_pdf = pmf * TrianglePDF(&world->triangles[hit.triangle_idx], ray.origin, hit.hit_point, hit.normal); 
				float mis_weight = PowerHeuristic(bsdf_pdf, light_pdf);

				Assert(light_pdf > 0);

				color += throughput * mat.emission * mis_weight;
			}
		}

		if (sample_lights) {
			u32 rand_idx = Rand(rng_state, 0, world->emissive_triangle_count-1);
			Triangle *light = &world->emissive_triangles[rand_idx];

			Vec3 uvw = RandomTriangle(rng_state);
			Vec3 light_point = uvw.x * light->v0 + uvw.y * light->v1 + uvw.z * light->v2;
			Vec3 light_normal = Normalize(uvw.x * light->n0 + uvw.y * light->n1 + uvw.z * light->n2);

			Vec3 light_direction = light_point - hit.hit_point;
			Vec3 l = local_basis * Normalize(light_direction);

			if (Dot(light_direction, light_normal) < 0 &&
			    Dot(light_direction, hit.normal) > 0 &&
			    !Occluded(world, light_point, hit.hit_point)) {
				float pmf = 1.f / (float)world->emissive_triangle_count;
				float light_pdf = pmf * TrianglePDF(light, hit.hit_point, light_point, light_normal); 
				float b_pdf = BSDFPDF(v, l, mat);
				float mis_weight = PowerHeuristic(light_pdf, b_pdf);

				Assert(light_pdf > 0);

				Material mat_light = world->materials[light->material_idx];
				color += throughput * mat_light.emission * BSDF(v, l, mat) * mis_weight / light_pdf;
			}
		}

		Sample sample = SampleBSDF(v, mat, rng_state);

		throughput *= sample.bsdf / sample.pdf;

		if (i > 3) {
			float roulette_prob = fmaxf(throughput.x, fmaxf(throughput.y, throughput.z));

			if (Rand(rng_state) > roulette_prob)
				break;

			throughput /= roulette_prob;
		}

		ray.origin = hit.hit_point;
		ray.direction = global_basis * sample.l;

		Assert(fabsf(Length(ray.direction) - 1) < .0001f);
		bsdf_pdf = sample.pdf;
	}

	return color;
}

int AddMaterial(World *world, Vec3 albedo, Vec3 emission, float roughness, float ior, float metallic, bool transparent)
{
	if (world->material_count >= MAX_MATERIALS)
		Panic("Too much materials");

	world->materials[world->material_count] = {
		albedo,
		emission,
		roughness,
		ior,
		metallic,
		transparent
	};

	return world->material_count++;
}

int AddDielectricMaterial(World *world, Vec3 albedo, Vec3 emission, float roughness)
{
	return AddMaterial(world, albedo, emission, roughness, 1.5f, 0, false);
}

int AddMetallicMaterial(World *world, Vec3 albedo, Vec3 emission, float roughness)
{
	return AddMaterial(world, albedo, emission, roughness, 0, 1, false);
}

int AddTransparentMaterial(World *world, Vec3 albedo, Vec3 emission, float ior)
{
	return AddMaterial(world, albedo, emission, 0, ior, false, true);
}

void AddSphere(World *world, Vec3 pos, float radius, int material_idx)
{
	if (world->sphere_count >= MAX_OBJECTS)
		Panic("Too much spheres");

	world->spheres[world->sphere_count++] = {
		material_idx,
		pos,
		radius
	};
}

void AddTriangle(World *world, Vec3 v0, Vec3 v1, Vec3 v2, int material_idx)
{
	if (world->triangle_count >= MAX_OBJECTS)
		Panic("Too much planes");

	Vec3 e0 = v1 - v0;
	Vec3 e1 = v2 - v0;
	Vec3 normal = Normalize(Cross(e0, e1));

	Triangle tri = {
		material_idx,
		v0, v1, v2,
		normal, normal, normal
	};

	world->triangles[world->triangle_count++] = tri;

	Vec3 emission = world->materials[material_idx].emission;

	if (emission.x > 0 || emission.y > 0 || emission.z > 0)
		world->emissive_triangles[world->emissive_triangle_count++] = tri;

}

void LoadCornellBox(World *world)
{
	world->sun_color = {};

	world->sphere_count = 0;
	world->triangle_count = 0;
	world->emissive_triangle_count = 0;
	world->material_count = 0;

	int khaki = AddMaterial(world, {0.725f, 0.71f, 0.68f}, {0, 0, 0}, 1, 1.5f, 0, false);
	int red = AddMaterial(world, {0.63f, 0.065f, 0.05f}, {0, 0, 0}, 1, 1.5f, 0, false);
	int green = AddMaterial(world, {0.14f, 0.45f, 0.091f}, {0, 0, 0}, 1, 1.5f, 0, false);
	int light = AddMaterial(world, {0, 0, 0}, {17, 12, 4}, 1, 1.5f, 0, false);

	AddTriangle(world, {-0.99f, -1, 0}, {1.04f, 0.99f, -0}, {-0.99f, 1.01f, 0}, khaki);
	AddTriangle(world, {-0.99f, -1, 0}, {1.04f, -1, -0}, {1.04f, 0.99f, -0}, khaki);
	AddTriangle(world, {1.04f, 1.02f, 1.99f}, {-0.99f, -1, 1.99f}, {-0.99f, 1.02f, 1.99f}, khaki);
	AddTriangle(world, {1.04f, 1.02f, 1.99f}, {1.04f, -1, 1.99f}, {-0.99f, -1, 1.99f}, khaki);
	AddTriangle(world, {1.04f, 0.99f, -0}, {1.04f, -1, 1.99f}, {1.04f, 1.02f, 1.99f}, khaki);
	AddTriangle(world, {1.04f, 0.99f, -0}, {1.04f, -1, -0}, {1.04f, -1, 1.99f}, khaki);
	AddTriangle(world, {-0.99f, -1, 0}, {1.04f, -1, 1.99f}, {1.04f, -1, -0}, green);
	AddTriangle(world, {-0.99f, -1, 0}, {-0.99f, -1, 1.99f}, {1.04f, -1, 1.99f}, green);
	AddTriangle(world, {-0.99f, 1.01f, 0}, {1.04f, 1.02f, 1.99f}, {-0.99f, 1.02f, 1.99f}, red);
	AddTriangle(world, {-0.99f, 1.01f, 0}, {1.04f, 0.99f, -0}, {1.04f, 1.02f, 1.99f}, red);
	AddTriangle(world, {-0.75f, -0.53f, 0.6f}, {0, -0.13f, 0.6f}, {-0.57f, 0.05f, 0.6f}, khaki);
	AddTriangle(world, {-0.57f, 0.05f, 0.6f}, {-0, -0.13f, 0}, {-0.57f, 0.05f, 0}, khaki);
	AddTriangle(world, {-0.75f, -0.53f, 0.6f}, {-0.57f, 0.05f, 0}, {-0.75f, -0.53f, 0}, khaki);
	AddTriangle(world, {-0.17f, -0.7f, 0.6f}, {-0.75f, -0.53f, 0}, {-0.17f, -0.7f, 0}, khaki);
	AddTriangle(world, {0, -0.13f, 0.6f}, {-0.17f, -0.7f, 0}, {-0, -0.13f, 0}, khaki);
	AddTriangle(world, {-0.75f, -0.53f, 0.6f}, {-0.17f, -0.7f, 0.6f}, {0, -0.13f, 0.6f}, khaki);
	AddTriangle(world, {-0.57f, 0.05f, 0.6f}, {0, -0.13f, 0.6f}, {-0, -0.13f, 0}, khaki);
	AddTriangle(world, {-0.75f, -0.53f, 0.6f}, {-0.57f, 0.05f, 0.6f}, {-0.57f, 0.05f, 0}, khaki);
	AddTriangle(world, {-0.17f, -0.7f, 0.6f}, {-0.75f, -0.53f, 0.6f}, {-0.75f, -0.53f, 0}, khaki);
	AddTriangle(world, {0, -0.13f, 0.6f}, {-0.17f, -0.7f, 0.6f}, {-0.17f, -0.7f, 0}, khaki);
	AddTriangle(world, {0.09f, -0.04f, 1.2f}, {0.49f, 0.71f, 1.2f}, {-0.09f, 0.53f, 1.2f}, khaki);
	AddTriangle(world, {-0.09f, 0.53f, 1.2f}, {0.49f, 0.71f, -0}, {-0.09f, 0.53f, 0}, khaki);
	AddTriangle(world, {0.49f, 0.71f, 1.2f}, {0.67f, 0.14f, -0}, {0.49f, 0.71f, -0}, khaki);
	AddTriangle(world, {0.67f, 0.14f, 1.2f}, {0.09f, -0.04f, -0}, {0.67f, 0.14f, -0}, khaki);
	AddTriangle(world, {0.09f, -0.04f, 1.2f}, {-0.09f, 0.53f, 0}, {0.09f, -0.04f, -0}, khaki);
	AddTriangle(world, {0.09f, -0.04f, 1.2f}, {0.67f, 0.14f, 1.2f}, {0.49f, 0.71f, 1.2f}, khaki);
	AddTriangle(world, {-0.09f, 0.53f, 1.2f}, {0.49f, 0.71f, 1.2f}, {0.49f, 0.71f, -0}, khaki);
	AddTriangle(world, {0.49f, 0.71f, 1.2f}, {0.67f, 0.14f, 1.2f}, {0.67f, 0.14f, -0}, khaki);
	AddTriangle(world, {0.67f, 0.14f, 1.2f}, {0.09f, -0.04f, 1.2f}, {0.09f, -0.04f, -0}, khaki);
	AddTriangle(world, {0.09f, -0.04f, 1.2f}, {-0.09f, 0.53f, 1.2f}, {-0.09f, 0.53f, 0}, khaki);
	AddTriangle(world, {0.22f, 0.24f, 1.98f}, {-0.16f, -0.23f, 1.98f}, {-0.16f, 0.24f, 1.98f}, light);
	AddTriangle(world, {0.22f, 0.24f, 1.98f}, {0.22f, -0.23f, 1.98f}, {-0.16f, -0.23f, 1.98f}, light);
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
			if (color.x < 0 || color.y < 0 || color.z < 0)
				color = {0, 0, 1};
			if (isnan(color.x) || isnan(color.y) || isnan(color.z))
				color = {0, 1, 0};
			Vec3 mapped = LinearToGamma(color, exposure);

			int pixel_pos = (v * width + u) * 3;

			byte ir = (byte)(255 * Clamp(mapped.x, 0, 1));
			byte ig = (byte)(255 * Clamp(mapped.y, 0, 1));
			byte ib = (byte)(255 * Clamp(mapped.z, 0, 1));

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
