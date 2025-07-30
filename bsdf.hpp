struct Material {
	Vec3 color;
	float roughness;
	float ior;
	float metallic;
};

struct Sample {
	Vec3 bsdf;
	float pdf;
	Vec3 l;
};

Vec3 CosineWeightedSample(u32 *rng_state)
{
	float r1 = Rand(rng_state);
	float r2 = Rand(rng_state);

	float phi = 2 * PI * r1;
	float sqrt_r2 = sqrtf(r2);

	float x = cosf(phi) * sqrt_r2;
	float y = sinf(phi) * sqrt_r2;
	float z = sqrtf(1 - r2);

	return Vec3{x, y, z};
}

float CosineWeightedPDF(Vec3 l)
{
	return fmaxf(l.z, 0) / PI;
}

/*
 * https://jcgt.org/published/0007/04/01/
 */
Vec3 GGXVNDFSample(Vec3 v, float roughness, u32 *rng_state)
{
	float r1 = Rand(rng_state);
	float r2 = Rand(rng_state);

	float alpha = roughness * roughness;

	Vec3 vh = Normalize(Vec3{alpha * v.x, alpha * v.y, v.z});
	float lensq = vh.x * vh.x + vh.y * vh.y;
	Vec3 t1 = lensq > 0 ? Vec3{-vh.y, vh.x, 0} / sqrtf(lensq) : Vec3{1, 0, 0};
	Vec3 t2 = Cross(vh, t1);
	float r = sqrtf(r1);
	float phi = 2 * PI * r2;

	float p1 = r*cosf(phi);
	float p2 = r*sinf(phi);
	float s = (1 + vh.z) / 2;
	p2 = (1 - s) * sqrtf(fmaxf(0, 1 - p1 * p1)) + s * p2;

	Vec3 n = p1*t1 + p2*t2 + sqrtf(fmaxf(0, 1 - p1*p1 - p2*p2))*vh;
	n = Normalize(Vec3{alpha*n.x, alpha*n.y, fmaxf(0, n.z)});

	return 2 * n * Dot(n, v) - v;
}

float GGXVNDFPDF(Vec3 v, Vec3 l, float roughness)
{
	Vec3 h = Normalize(v + l);

	if (Dot(h, v) <= 0)
		return 0;

	float alpha = roughness * roughness;
	float alpha2 = alpha * alpha;

	float ndf = alpha2 / (PI * powf(h.z * h.z * (alpha2 - 1) + 1, 2));

	float vis_v = 1 / (fabsf(v.z) + sqrtf(alpha2 + (1 - alpha2) * v.z * v.z));

	return ndf * vis_v / 2;
}

Vec3 Fresnel(float cosine, Vec3 f0)
{
	return f0 + (1 - f0) * powf(1 - cosine, 5);
}

/*
 * https://jcgt.org/published/0003/02/03/
 * https://google.github.io/filament/Filament.md.html#materialsystem
 */
Vec3 BSDF(Vec3 v, Vec3 l, Material mat)
{
	Assert(v.z > 0);

	if (l.z <= 0)
		return {};

	Vec3 h = Normalize(v + l);

	float alpha = mat.roughness * mat.roughness;
	float alpha2 = alpha * alpha;

	// GGX normal distribution function
	float ndf = alpha2 / (PI * powf(powf(h.z, 2) * (alpha2 - 1) + 1, 2));

	// Visibility function
	float vis_v = l.z * sqrtf(v.z * v.z * (1 - alpha2) + alpha2);
	float vis_l = v.z * sqrtf(l.z * l.z * (1 - alpha2) + alpha2);

	float vis = .5f / (vis_v + vis_l);

	// Fresnel term
	Vec3 dielectric_f0 = V3(powf((1 - mat.ior) / (1 + mat.ior), 2));

	Vec3 f0 = Lerp(dielectric_f0, mat.color, mat.metallic);

	Vec3 fresnel = Fresnel(Dot(h, v), f0);

	Vec3 diffuse = (1 - fresnel) * mat.color / PI * (1 - mat.metallic);

	Vec3 specular = fresnel * (vis * ndf);

	return (diffuse + specular) * l.z;
}

void GetWeights(Material mat, float *_cosine_weight, float *_vndf_weight)
{
	float cosine_weight = 1 - mat.metallic;
	float vndf_weight = 1 - cosine_weight * mat.roughness;

	float sum_weights = vndf_weight + cosine_weight;

	vndf_weight /= sum_weights;
	cosine_weight /= sum_weights;

	*_cosine_weight = cosine_weight;
	*_vndf_weight = vndf_weight;
}

float BSDFPDF(Vec3 v, Vec3 l, Material mat)
{
	float cosine_weight;
	float vndf_weight;

	GetWeights(mat, &cosine_weight, &vndf_weight);

	float cosine_pdf = CosineWeightedPDF(l);
	float vndf_pdf = GGXVNDFPDF(v, l, mat.roughness);

	return cosine_pdf * cosine_weight + vndf_pdf * vndf_weight;
}

Sample SampleBSDF(Vec3 v, Material mat, u32 *rng_state)
{
	Sample samp = {};
	float cosine_weight;
	float vndf_weight;

	GetWeights(mat, &cosine_weight, &vndf_weight);

	Vec3 l;
	if (Rand(rng_state) < cosine_weight) {
		l = CosineWeightedSample(rng_state);
	} else {
		l = GGXVNDFSample(v, mat.roughness, rng_state);
	}

	float pdf = BSDFPDF(v, l, mat);

	samp.bsdf = BSDF(v, l, mat);
	samp.pdf = pdf;
	samp.l = l;

	return samp;
}
