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

struct Object {
	ObjectType type;
	i32 material_idx;
	i32 light_idx;
	union {
		Triangle triangle;
		Sphere sphere;
	};
};

constexpr i32 MAX_OBJECTS = 64;
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

	int object_count;
	Object objects[MAX_OBJECTS];

	int light_count;
	Light lights[MAX_LIGHTS];

	int material_count;
	Material materials[MAX_MATERIALS];

	Vec3 sun_color;
};
