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

u32 rng_state = 69420;

double LinearToGamma(double color, double exposure)
{
	double gamma = 2.2;
	double mapped = 1. - exp(-color * exposure);
	return pow(mapped, 1. / gamma);
}


double Fresnel(double cosine, double refraction_index)
{
	double r0 = (1 - refraction_index) / (1 + refraction_index);
	r0 = r0*r0;
	return r0 + (1-r0)*pow((1 - cosine), 5);
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

Vec3 RayTrace(World *world, Ray ray)
{
	int max_bounces = 10;

	Vec3 color = {}, throughput = {1, 1, 1};

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
		Sample sample = SampleBSDF(v, mat);

		throughput *= sample.bsdf * fabs(sample.l.z) / sample.pdf;

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

	int gray = AddDielectricMaterial(&world, {.7, .7, .7}, {0, 0, 0}, 1);
	int red = AddDielectricMaterial(&world, {1, 0, 0}, {0, 0, 0}, .01);
	int light = AddDielectricMaterial(&world, {0, 0, 0}, {.5, .5, 5}, 1);
	int mirror = AddMetallicMaterial(&world, {.8, .8, .8}, {0, 0, 0}, .4);

	AddSphere(&world, {0, 0, 1}, 2, mirror);
	AddSphere(&world, {0, 3, -1.5}, .5, light);
	AddSphere(&world, {-1, -3, 0}, 1, red);
	AddPlane(&world, {0, 0, -1}, 2, gray);

	world.sun_color = {.5, .5, .8};

	int width = 630, height = 340;
	int n_samples = 500;
	double exposure = 1;
	double fov = 90;
	Vec3 camera_position = {-4, 3, 2};

	double fov_radians = fov * M_PI / 180;
	double aspect_ratio = (double)width/height;
	double viewport_distance = 1;
	double viewport_height = 2 * tan(fov_radians/2) * viewport_distance;
	double viewport_width = viewport_height * aspect_ratio;

	Vec3 looking_at = world.spheres[0].center;

	Vec3 forward = Normalize(looking_at - camera_position);
	Vec3 right = Normalize(Cross(forward, {0, 0, 1}));
	Vec3 down = Normalize(Cross(forward, right));

	Vec3 viewport_u = right * viewport_width;
	Vec3 viewport_v = down * viewport_height;

	Vec3 upper_left = camera_position +
		forward * viewport_distance
		- viewport_u/2
		- viewport_v/2;

	Vec3 du = viewport_u / width;
	Vec3 dv = viewport_v / height;

	byte *image_data = (byte *)malloc(sizeof(u8) * width * height * 3);

	double percent_row = 100 / (double)height;

	for (int v = 0; v < height; v++) {
		Log("Raytracing... %.0f%%\r", percent_row * v);

		for (int u = 0; u < width; u++) {
			Vec3 color = {};
			for (int i = 0; i < n_samples; i++) {
				double random_u = Rand();
				double random_v = Rand();

				Vec3 pixel_center = upper_left + du * (u + random_u) + dv * (v + random_v);
				Vec3 ray_direction = Normalize(pixel_center - camera_position);
				Ray ray = {camera_position, ray_direction};

				color += RayTrace(&world, ray);
			}

			color /= n_samples;

			double r = LinearToGamma(color.x, exposure);
			double g = LinearToGamma(color.y, exposure);
			double b = LinearToGamma(color.z, exposure);

			int pixel_pos = (v * width + u) * 3;

			byte ir = (byte)(255 * Clamp(r, 0, 1));
			byte ig = (byte)(255 * Clamp(g, 0, 1));
			byte ib = (byte)(255 * Clamp(b, 0, 1));

			image_data[pixel_pos + 2] = ir;
			image_data[pixel_pos + 1] = ig;
			image_data[pixel_pos + 0] = ib;
		}
	}
	Log("Raytracing... 100%%\n");

	WriteBMP("image.bmp", width, height, image_data);

	return 0;
}
