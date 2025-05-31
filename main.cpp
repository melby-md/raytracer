#define _USE_MATH_DEFINES
#include <limits>
#include <math.h>
#include <stdlib.h>

#include "bmp.hpp"
#include "bsdf.hpp"
#include "common.hpp"
#include "mathlib.hpp"

constexpr double INF = std::numeric_limits<double>::infinity();

struct Sphere {
	int material_idx;

	Vec3 center;
	double radius;
};

struct Plane {
	int material_idx;
	Vec3 normal;
	double d;
};

struct HitInfo {
	Vec3 hit_point;
	Vec3 normal;
	int material_idx;
};

struct Ray {
	Vec3 origin, direction;
};

constexpr int MAX_OBJECTS = 64;
constexpr int MAX_MATERIALS = 64;

struct World {
	int sphere_count;
	Sphere spheres[MAX_OBJECTS];

	int plane_count;
	Plane planes[MAX_OBJECTS];

	int material_count;
	Material materials[MAX_MATERIALS];

	Vec3 sun_color;
};

Vec2 RandomDisk(u32 *rng_state) {
	double theta = 2.0 * M_PI * Rand(rng_state);
	double r = sqrt(Rand(rng_state));

	double x = r * cos(theta);
	double y = r * sin(theta);

	return Vec2{x, y};
}

Vec3 LinearToGamma(Vec3 color, double exposure)
{
	Vec3 mapped = 1. - Exp((-color) * exposure);

	double g = 1 / 2.2;
	return Pow(mapped, Vec3{g, g, g});
}

double HitSphere(Sphere *sphere, Ray ray)
{
	Vec3 oc = sphere->center - ray.origin;
	double h = Dot(ray.direction, oc);
	double c = Dot(oc, oc) - sphere->radius*sphere->radius;
	double delta = h*h - c;

	if (delta < .001)
		return -1;

	double sqd = sqrt(delta);
	double distance = h - sqd;

	if (distance < .001) {
		distance = h + sqd;
		if (distance < .001)
			return -1;
	}

	return distance;
}

double HitPlane(Plane *plane, Ray ray)
{
	double denominator = Dot(plane->normal, ray.direction);
	if (fabs(denominator) < .001)
		return -1;

	double distance = (plane->d - Dot(ray.origin, plane->normal)) / denominator;
	if (distance < .001)
		return -1;

	return distance;
}

bool NearestHit(World *world, Ray ray, HitInfo *info)
{
	bool hit = false;
	double min_distance = INF;

	for (int i = 0; i < world->plane_count; i++) {
		Plane *plane = &world->planes[i];

		double t = HitPlane(plane, ray);

		if (t > 0 && t < min_distance) {
			info->hit_point = ray.origin + ray.direction * t;
			info->normal = plane->normal;
			info->material_idx = plane->material_idx;

			min_distance = t;
			hit = true;
		}
	}

	for (int i = 0; i < world->sphere_count; i++) {
		Sphere *sphere = &world->spheres[i];

		double t = HitSphere(sphere, ray);

		if (t > 0 && t < min_distance) {
			info->hit_point = ray.origin + ray.direction * t;
			info->normal = (info->hit_point - sphere->center) / sphere->radius;
			info->material_idx = sphere->material_idx;

			min_distance = t;
			hit = true;
		}
	}

	return hit;
}

Vec3 RayTrace(World *world, Ray ray, u32 *rng_state)
{
	int max_bounces = 10;

	Vec3 color = {};
	Vec3 throughput = {1, 1, 1};

	for (int i = 0; i < max_bounces; i++)
	{
		HitInfo hit;
		bool did_hit = NearestHit(world, ray, &hit);

		// hit the sky
		if (!did_hit) {
			color += throughput * world->sun_color;
			break;
		}

		Material mat = world->materials[hit.material_idx];
		color += throughput * mat.emission;

		Mat3 global_basis = OrthoNormalBasis(hit.normal);
		Mat3 local_basis = Transpose(global_basis);

		Vec3 v = local_basis * -ray.direction;
		Sample sample = SampleBSDF(v, mat, rng_state);

		throughput *= sample.bsdf * fabs(sample.l.z) / sample.pdf;

		if (i > 3) {
			double roulette_prob = Max(throughput.x, Max(throughput.y, throughput.z));

			if (Rand(rng_state) > roulette_prob)
				break;

			throughput /= roulette_prob;
		}

		ray.origin = hit.hit_point;
		ray.direction = global_basis * sample.l;
	}

	return color;
}

int AddMaterial(World *world, Vec3 albedo, Vec3 emission, double roughness, double ior, double metallic, bool transparent)
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

int AddDielectricMaterial(World *world, Vec3 albedo, Vec3 emission, double roughness)
{
	return AddMaterial(world, albedo, emission, roughness, 1.5, 0, false);
}

int AddMetallicMaterial(World *world, Vec3 albedo, Vec3 emission, double roughness)
{
	return AddMaterial(world, albedo, emission, roughness, 0, 1, false);
}

int AddTransparentMaterial(World *world, Vec3 albedo, Vec3 emission, double ior)
{
	return AddMaterial(world, albedo, emission, 0, ior, false, true);
}

void AddSphere(World *world, Vec3 pos, double radius, int material_idx)
{
	if (world->sphere_count >= MAX_OBJECTS)
		Panic("Too much spheres");

	world->spheres[world->sphere_count++] = {
		material_idx,
		pos,
		radius
	};
}

void AddPlane(World *world, Vec3 normal, double d, int material_idx)
{
	if (world->plane_count >= MAX_OBJECTS)
		Panic("Too much planes");

	world->planes[world->plane_count++] = {
		material_idx,
		normal,
		d
	};
}

int main()
{
	World world = {};

	int gray = AddDielectricMaterial(&world, Vec3{.7, .7, .7}, Vec3{0, 0, 0}, 1);
	int red = AddDielectricMaterial(&world, Vec3{1, 0, 0}, Vec3{0, 0, 0}, .01);
	int light = AddDielectricMaterial(&world, Vec3{0, 0, 0}, Vec3{.5, .5, 5}, 1);
	int mirror = AddMetallicMaterial(&world, Vec3{.8, .8, .8}, Vec3{0, 0, 0}, .4);

	AddSphere(&world, Vec3{0, 0, 1}, 2, mirror);
	AddSphere(&world, Vec3{0, 3, -1.5}, .5, light);
	AddSphere(&world, Vec3{-1, -3, 0}, 1, red);
	AddPlane(&world, Vec3{0, 0, -1}, 2, gray);

	world.sun_color = Vec3{.5, .5, .8};

	int width = 630, height = 340;
	int n_samples = 200;
	double exposure = 1;
	double fov = 90;
	Vec3 camera_position = {-4, 3, 2};
	Vec3 looking_at = world.spheres[0].center;
	Vec3 vup = {0, 0, 1};
	double defocus_angle = -2;
	double focus_dist = Length(looking_at - camera_position);

	double fov_radians = fov * M_PI / 180;
	double aspect_ratio = (double)width/height;
	double viewport_height = 2 * tan(fov_radians/2) * focus_dist;
	double viewport_width = viewport_height * aspect_ratio;

	Vec3 forward = Normalize(looking_at - camera_position);
	Vec3 right = Normalize(Cross(forward, vup));
	Vec3 up = Cross(right, forward);

	Vec3 viewport_u = right * viewport_width;
	Vec3 viewport_v = -up * viewport_height;

	Vec3 upper_left = camera_position +
		forward * focus_dist
		- viewport_u/2
		- viewport_v/2;

	Vec3 du = viewport_u / width;
	Vec3 dv = viewport_v / height;

	double defocus_angle_radians = defocus_angle * M_PI / 180;
	double defocus_radius = focus_dist * tan(defocus_angle_radians / 2);
	Vec3 defocus_disk_u = right * defocus_radius;
	Vec3 defocus_disk_v = up * defocus_radius;

	byte *image_data = (byte *)malloc(sizeof(u8) * width * height * 3);

	double percent_row = 100 / (double)height;
	double percent_done = 0;

	#pragma omp parallel for
	for (int v = 0; v < height; v++) {
		u32 rng_state = 69420 + v;

		for (int u = 0; u < width; u++) {
			Vec3 color = {};
			for (int i = 0; i < n_samples; i++) {
				double random_u = Rand(&rng_state);
				double random_v = Rand(&rng_state);
				
				Vec2 rand_disk = RandomDisk(&rng_state);
				Vec3 rand_defocus = camera_position + (rand_disk.x * defocus_disk_u) + (rand_disk.y * defocus_disk_v);

				Vec3 ray_origin = (defocus_angle <= 0) ? camera_position : rand_defocus;
				Vec3 pixel_center = upper_left + du * (u + random_u) + dv * (v + random_v);
				Vec3 ray_direction = Normalize(pixel_center - ray_origin);
				Ray ray = {ray_origin, ray_direction};

				color += RayTrace(&world, ray, &rng_state);
			}

			color /= n_samples;
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
