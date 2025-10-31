struct BSample {
	Vec3 bsdf;
	float pdf;
	Vec3 l;
};

struct Material;

static Vec3    BSDF(Vec3, Vec3, Material *);
static float   BSDFPDF(Vec3, Vec3, Material *);
static BSample SampleBSDF(Vec3, Material *, u32 *);
