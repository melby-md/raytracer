enum ObjectType {
	OBJ_TRIANGLE,
	OBJ_SPHERE
};

struct Light {
	Vec3 color;
	i32 object_idx;
};

struct Triangle {
	Vec3 v0, v1, v2;
	Vec3 n0, n1, n2;
};

struct Sphere {
	Vec3 center;
	float radius;
};

struct Material {
	Vec3 color;
	float roughness;
	float ior;
	float metallic;
};

struct Object {
	ObjectType type;
	i32 material_idx;
	i32 light_idx;
	union {
		Triangle triangle;
		Sphere sphere;
	};
};

struct BVHNode;

struct BVHTree {
	BVHNode *nodes;
	i32 max_nodes;
	Object *objects;
};

constexpr i32 MAX_OBJECTS = 16384;
constexpr i32 MAX_MATERIALS = 64;
constexpr i32 MAX_LIGHTS = 64;

struct Scene {
	Vec3 camera;
	Vec3 up;
	Vec3 look_at;
	float fov;
	float defocus_angle;
	float exposure;
	i16 width, height;
	i16 samples;

	i32 object_count;
	Object objects[MAX_OBJECTS];

	i32 light_count;
	Light lights[MAX_LIGHTS];

	i32 material_count;
	Material materials[MAX_MATERIALS];

	Vec3 sky_box_color;

	BVHTree bvh;
};

static float Rand(u32 *);
