#ifndef _BSDF_HPP
#define _BSDF_HPP
#include "common.hpp"
#include "mathlib.hpp"

struct Material {
	Vec3 albedo;
	Vec3 emission;

	double roughness;
	double ior;
	double metallic;

	bool transparent;
};


struct Sample {
	Vec3 bsdf;
	double pdf;
	Vec3 l;
};

Vec3 CosineWeightedSample(u32 *rng_state)
{
	double r1 = Rand(rng_state);
	double r2 = Rand(rng_state);

	double phi = 2 * M_PI * r1;
	double sqrt_r2 = sqrt(r2);

	double x = cos(phi) * sqrt_r2;
	double y = sin(phi) * sqrt_r2;
	double z = sqrt(1 - r2);

	return Vec3{x, y, z};
}

double CosineWeightedPDF(Vec3 l)
{
	return l.z / M_PI;
}

/*
 * https://jcgt.org/published/0007/04/01/
 */
Vec3 GGXVNDFSample(Vec3 v, double roughness, u32 *rng_state)
{
	double r1 = Rand(rng_state);
	double r2 = Rand(rng_state);

	double alpha = roughness * roughness;

	Vec3 vh = Normalize(Vec3{alpha * v.x, alpha * v.y, v.z});
	double lensq = vh.x * vh.x + vh.y * vh.y;
	Vec3 t1 = lensq > 0 ? Vec3{-vh.y, vh.x, 0} / sqrt(lensq) : Vec3{1, 0, 0};
	Vec3 t2 = Cross(vh, t1);
	double r = sqrt(r1);
	double phi = 2 * M_PI * r2;

	double p1 = r*cos(phi);
	double p2 = r*sin(phi);
	double s = .5 * (1 + vh.z);
	p2 = (1 - s) * sqrt(Max(0.0, 1.0 - p1 * p1)) + s * p2;

	Vec3 n = p1*t1 + p2*t2 + sqrt(Max(0, 1 - p1*p1 - p2*p2))*vh;
	n = Normalize(Vec3{alpha*n.x, alpha*n.y, Max(0, n.z)});

	return 2 * n * Dot(n, v) - v;
}

double GGXVNDFPDF(Vec3 v, Vec3 l, double roughness)
{
	double alpha = roughness * roughness;
	double alpha2 = alpha * alpha;

	Vec3 h = Normalize(v + l);

	double ndf = alpha2 / (M_PI * pow(h.z * h.z * (alpha2 - 1) + 1, 2));

	double vis_v = 1 / (fabs(v.z) + sqrt(alpha2 + (1.0 - alpha2) * v.z * v.z));

	return ndf * vis_v * .5;
}

Vec3 Fresnel(double cosine, Vec3 f0)
{
	return f0 + (Vec3{1, 1, 1} - f0) * pow(1 - cosine, 5);
}

/*
 * https://jcgt.org/published/0003/02/03/
 * https://google.github.io/filament/Filament.md.html#materialsystem
 */
Vec3 BSDF(Vec3 v, Vec3 l, Material mat)
{
	Vec3 h = Normalize(v + l);

	double alpha = mat.roughness * mat.roughness;
	double alpha2 = alpha * alpha;

	// GGX normal distribution function
	double n_dot_h = fabs(h.z);
	double ndf = alpha2 / (M_PI * pow(pow(n_dot_h, 2) * (alpha2 - 1) + 1, 2));

	// Visibility function
	double n_dot_v = fabs(v.z);
	double n_dot_l = fabs(l.z);

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

Sample SampleBSDF(Vec3 v, Material mat, u32 *rng_state)
{
	double cosine_weight = 1 - mat.metallic;
	double vndf_weight = 1 - cosine_weight * mat.roughness;

	double sum_weights = vndf_weight + cosine_weight;

	vndf_weight /= sum_weights;
	cosine_weight /= sum_weights;

	Vec3 l;
	if (Rand(rng_state) < cosine_weight) {
		l = CosineWeightedSample(rng_state);
	} else {
		l = GGXVNDFSample(v, mat.roughness, rng_state);
	}

	double cosine_pdf = CosineWeightedPDF(l);
	double vndf_pdf = GGXVNDFPDF(v, l, mat.roughness);

	double pdf = cosine_pdf * cosine_weight + vndf_pdf * vndf_weight;

	return Sample{
		BSDF(v, l, mat),
		pdf,
		l
	};
}

#endif
