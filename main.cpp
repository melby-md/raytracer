#include <limits>
#include <math.h>
#include <stdlib.h>

#include "bmp.hpp"
#include "common.hpp"
#include "mathlib.hpp"

constexpr double INF = std::numeric_limits<double>::infinity();

struct Material {
	Vec3 albedo;
	Vec3 emission;

	double roughness;
	double ior;
	double metallic;

	bool transparent;
};

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

double LinearToGamma(double color, double exposure)
{
	double gamma = 2.2;
	double mapped = 1. - exp(-color * exposure);
	return pow(mapped, 1. / gamma);
}

double Rand(double min, double max)
{
	// Xorshift32
	static u32 state = 69420;

	state ^= state << 13;
	state ^= state >> 17;
	state ^= state << 5;

	double f = (double)state / UINT32_MAX;
	return min + f * (max - min);
}

double Fresnel(double cosine, double refraction_index)
{
	double r0 = (1 - refraction_index) / (1 + refraction_index);
	r0 = r0*r0;
	return r0 + (1-r0)*pow((1 - cosine), 5);
}

Vec3 Fresnel(double cosine, Vec3 f0)
{
	return f0 + (Vec3{1, 1, 1} - f0) * pow(1 - cosine, 5);
}

/*
 * https://jcgt.org/published/0003/02/03/
 * https://registry.khronos.org/glTF/specs/2.0/glTF-2.0.html#appendix-b-brdf-implementation
 */
Vec3 BRDF(Vec3 l, Vec3 n, Vec3 v, Material mat)
{
	Vec3 h = Normalize(v + l);

	double alpha = mat.roughness * mat.roughness;
	double alpha2 = alpha * alpha;

	// GGX normal distribution function
	double n_dot_h = fabs(Dot(n, h));
	double ndf = alpha2 / (M_PI * pow(pow(n_dot_h, 2) * (alpha2 - 1) + 1, 2));

	// Visibility function
	double n_dot_v = fabs(Dot(n, v));
	double n_dot_l = fabs(Dot(n, l));

	double vis_v = n_dot_l * sqrt(n_dot_v * n_dot_v * (1 - alpha2) + alpha2);
	double vis_l = n_dot_v * sqrt(n_dot_l * n_dot_l * (1 - alpha2) + alpha2);

	double vis = .5 / (vis_v + vis_l);

	// Fresnel term
	double r0 = pow((1 - mat.ior) / (1 + mat.ior), 2);
	Vec3 dielectric_fresnel = Vec3{r0, r0, r0};

	Vec3 f0 = Lerp(dielectric_fresnel, mat.albedo, mat.metallic);

	double h_dot_v = fabs(Dot(h, v));
	Vec3 fresnel = Fresnel(h_dot_v, f0);

	Vec3 diffuse = (Vec3{1, 1, 1} - fresnel) * mat.albedo / M_PI * (1 - mat.metallic);

	Vec3 specular = fresnel * (vis * ndf);

	return diffuse + specular;
}

double CosineWeightedPDF(Vec3 wi)
{
	return Max(0, wi.z) / M_PI;
}

Vec3 CosineWeightedSample()
{
	double r1 = Rand(0, 1);
	double r2 = Rand(0, 1);

	double phi = 2 * M_PI * r1;
	double sqrt_r2 = sqrt(r2);

	double x = cos(phi) * sqrt_r2;
	double y = sin(phi) * sqrt_r2;
	double z = sqrt(1 - r2);

	return Vec3{x, y, z};
}

/*
 * https://jcgt.org/published/0007/04/01/
 */
Vec3 GGXVNDFSample(Vec3 wo, double roughness)
{
	double r1 = Rand(0, 1);
	double r2 = Rand(0, 1);

	double alpha = roughness * roughness;

	Vec3 vh = Normalize(Vec3{alpha * wo.x, alpha * wo.y, wo.z});
	double lensq = vh.x * vh.x + vh.y * vh.y;
	Vec3 t1 = lensq > 0 ? Vec3{-vh.y, vh.x, 0} / sqrt(lensq) : Vec3{1, 0, 0};
	Vec3 t2 = Cross(t1, vh);
	double r = sqrt(r1);
	double phi = 2 * M_PI * r2;

	double p1 = r*cos(phi);
	double p2 = r*sin(phi);
	double s = .5 * (1 + vh.z);
	p2 = (1 - s) * sqrt(Max(0.0, 1.0 - p1 * p1)) + s * p2;

	Vec3 n = p1*t1 + p2*t2 + sqrt(Max(0, 1 - p1*p1 - p2*p2))*vh;
	n = Normalize(Vec3{alpha*n.x, alpha*n.y, Max(0, n.z)});

	return 2 * n * Dot(n, wo) - wo;
}

double GGXVNDFPDF(Vec3 wo, Vec3 wi, double roughness)
{
	double alpha = roughness * roughness;
	double alpha2 = alpha * alpha;

	Vec3 h = Normalize(wo + wi);

	double ndf = alpha2 / (M_PI * pow(h.z * h.z * (alpha2 - 1) + 1, 2));

	double geometry_v = 2.0 * fabs(wo.z) / (fabs(wo.z) + sqrt(alpha2 + (1.0 - alpha2) * wo.z * wo.z));

	return ndf * geometry_v / (4 * wo.z);
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

	Vec3 color = {}, attenuation = {1, 1, 1};

	for (int i = 0; i < max_bounces; i++)
	{
		HitInfo hit;
		bool did_hit = NearestHit(world, ray, &hit);

		// hit the sky
		if (!did_hit) {
			color += attenuation * world->sun_color;
			break;
		}

		Material mat = world->materials[hit.material_idx];

		Vec3 next_direction;
		Vec3 attenuation_factor;
		if (mat.transparent) {
			bool front_face = Dot(ray.direction, hit.normal) < 0;

			if (!front_face) {
				hit.normal = -hit.normal;
			}

			double ri = front_face ? (1/mat.ior) : mat.ior;

			double cos_theta = Min(Dot(-ray.direction, hit.normal), 1);
			double sin_theta = sqrt(1 - cos_theta*cos_theta);

			bool cannot_refract = ri * sin_theta > 1.0;

			if (cannot_refract || Fresnel(cos_theta, ri) > Rand(0, 1)) {
				attenuation_factor = {1, 1, 1};
				next_direction = Reflect(ray.direction, hit.normal);
			} else {
				attenuation_factor = mat.albedo;
				next_direction = Refract(ray.direction, hit.normal, ri);
			}

		} else {
			Mat3 global_basis = OrthoNormalBasis(hit.normal);
			Mat3 local_basis = Transpose(global_basis);

			Vec3 wi;
			Vec3 wo = local_basis * -ray.direction;
			Vec3 cosine_sample = CosineWeightedSample();
			Vec3 vndf_sample = GGXVNDFSample(wo, mat.roughness);

			double cosine_weight = 1 - mat.metallic;
			double vndf_weight = 1 - cosine_weight * mat.roughness;

			double sum_weights = vndf_weight + cosine_weight;

			vndf_weight /= sum_weights;
			cosine_weight /= sum_weights;

			if (Rand(0, 1) < cosine_weight) {
				wi = cosine_sample;
			} else {
				wi = vndf_sample;
			}

			double cosine_pdf = CosineWeightedPDF(wi);
			double vndf_pdf = GGXVNDFPDF(wo, wi, mat.roughness);

			next_direction = global_basis * wi;

			double cos_theta = fabs(Dot(next_direction, hit.normal));

			double pdf = cosine_pdf * cosine_weight + vndf_pdf * vndf_weight;

			attenuation_factor = BRDF(next_direction, hit.normal, -ray.direction, mat) * cos_theta / pdf;
		}


		ray = {hit.hit_point, next_direction};

		color += attenuation * mat.emission;
		attenuation *= attenuation_factor;
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
	AddPlane(&world, {0, 0, 1}, 2, gray);

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

	Vec3 upper_left_pixel = upper_left + (du + dv)/2;

	byte *image_data = (byte *)malloc(sizeof(u8) * width * height * 3);

	double percent_row = 100 / (double)height;

	for (int v = 0; v < height; v++) {
		Log("Raytracing... %.0f%%\r", percent_row * v);

		for (int u = 0; u < width; u++) {
			Vec3 color = {};
			for (int i = 0; i < n_samples; i++) {
				double random_u = Rand(-.5, .5);
				double random_v = Rand(-.5, .5);

				Vec3 pixel_center = upper_left_pixel + du * (u + random_u) + dv * (v + random_v);
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
