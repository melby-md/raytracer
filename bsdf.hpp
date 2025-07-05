struct Material {
	Vec3 albedo;
	Vec3 emission;

	float roughness;
	float ior;
	float metallic;

	bool transmissive;
};

struct Sample {
	Vec3 bsdf;
	float pdf;
	Vec3 l;
};

bool IsSpecular(Material *mat)
{
	return mat->transmissive;
}

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
	return l.z / PI;
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
	float alpha = roughness * roughness;
	float alpha2 = alpha * alpha;

	Vec3 h = Normalize(v + l);

	float ndf = alpha2 / (PI * powf(h.z * h.z * (alpha2 - 1) + 1, 2));

	float vis_v = 1 / (fabsf(v.z) + sqrtf(alpha2 + (1 - alpha2) * v.z * v.z));

	return ndf * vis_v / 2;
}

Vec3 Fresnel(float cosine, Vec3 f0)
{
	return f0 + (Vec3{1, 1, 1} - f0) * powf(1 - cosine, 5);
}

float Fresnel(float cosine, float ri)
{
	cosine = fminf(cosine, 1);
	float sin2theta = (1 - cosine*cosine);

	if (sin2theta * ri*ri >= 1)
		return 1;

	float r0 = (1 - ri) / (1 + ri);
	r0 = r0*r0;
	return r0 + (1 - r0) * powf(1 - cosine, 5);
}

/*
 * https://jcgt.org/published/0003/02/03/
 * https://google.github.io/filament/Filament.md.html#materialsystem
 */
Vec3 BSDF(Vec3 v, Vec3 l, Material mat)
{
	Vec3 h = Normalize(v + l);

	float alpha = mat.roughness * mat.roughness;
	float alpha2 = alpha * alpha;

	// GGX normal distribution function
	float n_dot_h = fabsf(h.z);
	float ndf = alpha2 / (PI * powf(powf(n_dot_h, 2) * (alpha2 - 1) + 1, 2));

	// Visibility function
	float n_dot_v = fabsf(v.z);
	float n_dot_l = fabsf(l.z);

	float vis_v = n_dot_l * sqrtf(n_dot_v * n_dot_v * (1 - alpha2) + alpha2);
	float vis_l = n_dot_v * sqrtf(n_dot_l * n_dot_l * (1 - alpha2) + alpha2);

	float vis = .5f / (vis_v + vis_l);

	// Fresnel term
	float r0 = powf((1 - mat.ior) / (1 + mat.ior), 2);
	Vec3 dielectric_fresnel = Vec3{r0, r0, r0};

	Vec3 f0 = Lerp(dielectric_fresnel, mat.albedo, mat.metallic);

	float h_dot_v = fabsf(Dot(h, v));
	Vec3 fresnel = Fresnel(h_dot_v, f0);

	Vec3 diffuse = (Vec3{1, 1, 1} - fresnel) * mat.albedo / PI * (1 - mat.metallic);

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
	if (IsSpecular(&mat))
		return 0;

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
