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
 * https://doi.org/10.1111/cgf.14867
 */
Vec3 GGXVNDFSample(Vec3 v, float alpha, u32 *rng_state)
{
	float r1 = Rand(rng_state);
	float r2 = Rand(rng_state);

	Vec3 vh = Normalize(Vec3{alpha * v.x, alpha * v.y, v.z});

	float phi = 2 * PI * r1;
	float z = fmaf((1 - r2), (1 + vh.z), -vh.z);
	float sin_theta = sqrtf(Clamp(1 - z * z, 0, 1));
	float x = sin_theta * cosf(phi);
	float y = sin_theta * sinf(phi);
	Vec3 cap = {x, y, z};
	Vec3 h = cap + vh;
	Vec3 n = Normalize(Vec3{alpha*h.x, alpha*h.y, h.z});

	return 2 * n * Dot(n, v) - v;

}

float GGXVNDFPDF(Vec3 v, Vec3 l, float alpha)
{
	Vec3 h = Normalize(v + l);

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
Vec3 BSDF(Vec3 v, Vec3 l, Material *mat)
{
	Assert(v.z > 0);

	if (l.z <= 0)
		return {};

	Vec3 h = Normalize(v + l);

	float alpha2 = mat->roughness * mat->roughness;

	// GGX normal distribution function
	float ndf = alpha2 / (PI * powf(powf(h.z, 2) * (alpha2 - 1) + 1, 2));

	// Visibility function
	float vis_v = l.z * sqrtf(v.z * v.z * (1 - alpha2) + alpha2);
	float vis_l = v.z * sqrtf(l.z * l.z * (1 - alpha2) + alpha2);

	float vis = .5f / (vis_v + vis_l);

	// Fresnel term
	Vec3 dielectric_f0 = V3(powf((1 - mat->ior) / (1 + mat->ior), 2));

	Vec3 f0 = Lerp(dielectric_f0, mat->color, mat->metallic);

	Vec3 fresnel = Fresnel(Dot(h, v), f0);

	Vec3 diffuse = (1 - fresnel) * mat->color / PI * (1 - mat->metallic);

	Vec3 specular = fresnel * (vis * ndf);

	return (diffuse + specular) * l.z;
}

void GetWeights(Material *mat, float *_cosine_weight, float *_vndf_weight)
{
	float cosine_weight = 1 - mat->metallic;
	float vndf_weight = 1;

	float sum_weights = vndf_weight + cosine_weight;

	vndf_weight /= sum_weights;
	cosine_weight /= sum_weights;

	*_cosine_weight = cosine_weight;
	*_vndf_weight = vndf_weight;
}

float BSDFPDF(Vec3 v, Vec3 l, Material *mat)
{
	float cosine_weight;
	float vndf_weight;

	GetWeights(mat, &cosine_weight, &vndf_weight);

	float cosine_pdf = CosineWeightedPDF(l);
	float vndf_pdf = GGXVNDFPDF(v, l, mat->roughness);

	return cosine_pdf * cosine_weight + vndf_pdf * vndf_weight;
}

Sample SampleBSDF(Vec3 v, Material *mat, u32 *rng_state)
{
	Sample samp = {};
	float cosine_weight;
	float vndf_weight;

	GetWeights(mat, &cosine_weight, &vndf_weight);

	Vec3 l;
	if (Rand(rng_state) < cosine_weight) {
		l = CosineWeightedSample(rng_state);
	} else {
		l = GGXVNDFSample(v, mat->roughness, rng_state);
	}

	float pdf = BSDFPDF(v, l, mat);

	samp.bsdf = BSDF(v, l, mat);
	samp.pdf = pdf;
	samp.l = l;

	return samp;
}
