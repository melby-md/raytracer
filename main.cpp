#include <limits>
#include <math.h>
#include <stdlib.h>

#include "bmp.hpp"
#include "common.hpp"
#include "vector.hpp"

constexpr double INF = std::numeric_limits<double>::infinity();

struct Material {
	double smoothness;
	Vec3 reflection_color;
	Vec3 emission_color;
};

enum class Geometry {
	SPHERE,
	PLANE
};

struct Object {
	Geometry type;
	int material_idx;
	Vec3 v;
	double s;
};

struct Ray {
	Vec3 origin, direction;
};

constexpr int MAX_OBJECTS = 256;
constexpr int MAX_MATERIALS = 64;

struct World {
	int object_count;
	Object objects[MAX_OBJECTS];

	int material_count;
	Material materials[MAX_MATERIALS];

	Vec3 sun_direction;
	Vec3 sun_color;
};

double LinearToGamma(double linear_component)
{
    if (linear_component > 0)
        return sqrt(linear_component);

    return 0;
}

double Rand(double fMin, double fMax)
{
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

Vec3 RandomUnitVector() {
	for (;;) {
		Vec3 p = {Rand(-1, 1), Rand(-1, 1), Rand(-1, 1)};
		double lensq = Length2(p);
		if (1e-160 < lensq && lensq <= 1)
			return p / sqrt(lensq);
	}
}

double HitSphere(Object *sphere, Ray ray)
{
	Vec3 oc = sphere->v - ray.origin;
	double a = Dot(ray.direction, ray.direction);
	double h = Dot(ray.direction, oc);
	double c = Dot(oc, oc) - sphere->s*sphere->s;
	double delta = h*h - a*c;

	if (delta < .001)
		return -1;

	double sqd = sqrt(delta);
	double distance = (h - sqd)/a;

	if (distance < .001) {
		distance = (h + sqd)/a;
		if (distance < .001)
			return -1;
	}

	return distance;
}

double HitPlane(Object *plane, Ray ray)
{
	double denominator = Dot(plane->v, ray.direction);
	if (denominator < .001)
		return -1;

	double distance = (-plane->s - Dot(ray.origin, plane->v)) / denominator;
	if (distance < .001)
		return -1;

	return distance;
}

Vec3 GetNormal(Object *object, Vec3 hit_point)
{
	switch (object->type) {
	case Geometry::PLANE:
		return -object->v;
	case Geometry::SPHERE:
		return Normalize(hit_point - object->v);
	}

	return {0, 0, 0};
}

Object *NearestObject(World *world, Ray ray, double *t)
{
	Object *nearest = nullptr;
	double min = INF;

	for (int i = 0; i < world->object_count; i++) {
		Object *obj = &world->objects[i];
		double distance = -1;

		switch (obj->type) {
		case Geometry::PLANE:
			distance = HitPlane(obj, ray);
			break;
		case Geometry::SPHERE:
			distance = HitSphere(obj, ray);
		}

		if (distance > 0 && distance < min) {
			min = distance;
			nearest = obj;
		}
	}

	*t = min;
	return nearest;
}

Vec3 RayTrace(World *world, Ray ray)
{
	int max_bounces = 10;

	Vec3 color = {}, attenuation = {1, 1, 1};

	for (int i = 0; i < max_bounces; i++)
	{
		double t;
		Object *obj = NearestObject(world, ray, &t);

		// hit the sky
		if (obj == nullptr) {
			color += attenuation * world->sun_color;
			break;
		}

		Vec3 hit_point = ray.origin + ray.direction*t;
		Vec3 normal = GetNormal(obj, hit_point);

		Vec3 scattered = Normalize(normal + RandomUnitVector());
		Vec3 reflection = Normalize(Reflect(ray.direction, normal));

		Material mat = world->materials[obj->material_idx];
		Vec3 fuzzy_reflection = Lerp(scattered, reflection, mat.smoothness);

		color += attenuation * mat.emission_color;
		attenuation *= mat.reflection_color;

		ray = {hit_point, fuzzy_reflection};
	}

	return color;
}

int AddMaterial(World *world, Vec3 color, Vec3 emission, double smoothness)
{
	world->materials[world->material_count] = {
		smoothness,
		color,
		emission
	};

	return world->material_count++;
}

void AddSphere(World *world, Vec3 pos, double radius, int material_idx)
{
	world->objects[world->object_count++] = {
		Geometry::SPHERE,
		material_idx,
		pos,
		radius
	};
}

void AddPlane(World *world, Vec3 normal, double d, int material_idx)
{
	world->objects[world->object_count++] = {
		Geometry::PLANE,
		material_idx,
		normal,
		d
	};
}

int main()
{
	World world = {};

	int gray = AddMaterial(&world, {.7, .7, .7}, {0, 0, 0}, 0);
	int red = AddMaterial(&world, {.8, 0, 0}, {0, 0, 0}, 0);
	int light = AddMaterial(&world, {0, 0, 0}, {5, 0, 5}, 0);
	int mirror = AddMaterial(&world, {.9, .9, .9}, {0, 0, 0}, .98);

	AddSphere(&world, {0, 0, 1}, 2, mirror);
	AddSphere(&world, {0, 3, 0}, .5, red);
	AddSphere(&world, {-1, -3, 0}, 1, light);
	AddPlane(&world, {0, 0, -1}, -2, gray);

	world.sun_direction = {-5, 0, 0};
	world.sun_color = {1, 1, 1};

	int width = 630, height = 340;
	// 90 degrees
	double fov = M_PI/2;
	Vec3 camera_position = {-4, 3, 2};

	double aspect_ratio = (double)width/height;
	double viewport_distance = 1;
	double viewport_height = 2 * tan(fov/2) * viewport_distance;
	double viewport_width = viewport_height * aspect_ratio;

	Vec3 looking_at = world.objects[0].v;

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

	Vec3 upper_left_pixel = upper_left + (du + dv)/2;

	byte *image_data = (byte *)malloc(sizeof(u8) * width * height * 3);

	int n_samples = 100;

	for (int v = 0; v < height; v++) {
		for (int u = 0; u < width; u++) {
			Vec3 color = {};
			for (int i = 0; i < n_samples; i++) {
				double random_u = Rand(-.5, .5);
				double random_v = Rand(-.5, .5);

				Vec3 pixel_center = upper_left_pixel + du * (u + random_u) + dv * (v + random_v);
				Vec3 ray_direction = pixel_center - camera_position;
				Ray ray = {camera_position, ray_direction};

				color += RayTrace(&world, ray);
			}

			color /= n_samples;

			double r = LinearToGamma(color.x);
			double g = LinearToGamma(color.y);
			double b = LinearToGamma(color.z);

			int pixel_pos = (v * width + u) * 3;

			byte ir = (byte)(255 * Clamp(r, 0, 1));
			byte ig = (byte)(255 * Clamp(g, 0, 1));
			byte ib = (byte)(255 * Clamp(b, 0, 1));

			image_data[pixel_pos + 2] = ir;
			image_data[pixel_pos + 1] = ig;
			image_data[pixel_pos + 0] = ib;
		}
	}

	WriteBMP("image.bmp", width, height, image_data);

	return 0;
}
