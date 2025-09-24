enum TokenType {
	TOK_NIL,
	TOK_STRING,
	TOK_L_BRACE, TOK_R_BRACE, TOK_L_BRACKET, TOK_R_BRACKET,
	TOK_END
};

struct Token {
	TokenType type;
	String lexeme;
};

struct Lexer {
	const char *file_name;
	char *input;
	isize pos;
	isize prev_pos;
};

static void ReportError(Lexer *lexer, const char *msg)
{
	Log(
		"ERROR:%s[%ld]: %s\n",
		lexer->file_name,
		lexer->prev_pos+1,
		msg
	);

	exit(1);
}

static bool IsSpace(char c)
{
	return c == ' ' || c == '\t' || c == '\r' || c == '\n';
}

static bool IsAlpha(char c)
{
	switch (c) {
	case '\0': case '{': case '}': case '[': case ']': case '#':
		return false;
	default:
		return !IsSpace(c);
	}
}

static bool StringEquals(String a, String b)
{
	if (a.length != b.length)
		return false;
	return memcmp(a.data, b.data, a.length) == 0;
}

static Token NextToken(Lexer *lexer)
{
	isize start;
	TokenType type = TOK_NIL;
	do {
		start = lexer->pos;
		char c = lexer->input[lexer->pos++];
		switch (c) {
		case ' ': case '\t': case '\r': case '\n':
			while (IsSpace(lexer->input[lexer->pos]))
				lexer->pos++;
			break;
		
		case '#':
			while (lexer->input[lexer->pos++] != '\n');
			break;

		case '\0':
			type = TOK_END;
			break;

		case '{':
			type = TOK_L_BRACE;
			break;

		case '}':
			type = TOK_R_BRACE;
			break;

		case '[':
			type = TOK_L_BRACKET;
			break;

		case ']':
			type = TOK_R_BRACKET;
			break;

		default:
			while (IsAlpha(lexer->input[lexer->pos]))
				lexer->pos++;

			type = TOK_STRING;
		}
	} while (type == TOK_NIL);

	lexer->prev_pos = start;
	return {type, {&lexer->input[start], lexer->pos - start}};
}

static bool ReadCmd(Lexer *lexer, String *cmd)
{
	Token t = NextToken(lexer);

	isize start = lexer->prev_pos;

	if (t.type == TOK_END) {
		return false;
	} else if (t.type == TOK_STRING) {
		if (NextToken(lexer).type != TOK_L_BRACE)
			ReportError(lexer, "Expected '{'");

		*cmd = t.lexeme;
		lexer->prev_pos = start;
		return true;
	}

	ReportError(lexer, "Expected command");
	return false;
}

static bool ReadKey(Lexer *lexer, String *key)
{
	Token t = NextToken(lexer);

	if (t.type == TOK_R_BRACE) {
		return false;
	} else if (t.type == TOK_STRING) {
		*key = t.lexeme;
		return true;
	}

	ReportError(lexer, "Expected key");
	return false;
}

static String ReadString(Lexer *lexer)
{
	Token t = NextToken(lexer);

	if (t.type == TOK_STRING) {
		return t.lexeme;
	}

	ReportError(lexer, "Expected string");
	return {};
}

static int StringToNumber(String s, float *n)
{
	char *end = nullptr;
	errno = 0;
	*n = strtof(s.data, &end);
	if (end - s.data != s.length || errno)
		return -1;
	return 0;
}

static float ReadNumber(Lexer *lexer)
{
	String s = ReadString(lexer);
	float n;
	if (StringToNumber(s, &n) < 0) {
		ReportError(lexer, "Invalid number");
	}
	return n;
}

static i16 ReadI16(Lexer *lexer)
{
	String s = ReadString(lexer);
	long n;
	char *end = nullptr;
	errno = 0;
	n = strtol(s.data, &end, 10);
	if (end - s.data != s.length || errno)
		ReportError(lexer, "Invalid integer");
	if (n > (1 << 16) - 1 || n < 0)
		ReportError(lexer, "Out of bounds integer");
	return (i16)n;
}

static void BeginArray(Lexer *lexer)
{
	Token t = NextToken(lexer);

	if (t.type != TOK_L_BRACKET)
		ReportError(lexer, "Expected array");
}

static bool EndArray(Lexer *lexer)
{
	isize start = lexer->pos;
	Token t = NextToken(lexer);

	if (t.type == TOK_R_BRACKET)
		return true;

	lexer->pos = start;
	return false;
}

static Vec3 ReadVec3(Lexer *lexer)
{
	Vec3 v;

	Token t = NextToken(lexer);
	if (t.type != TOK_L_BRACKET)
		ReportError(lexer, "Expected array");

	v.x = ReadNumber(lexer);
	v.y = ReadNumber(lexer);
	v.z = ReadNumber(lexer);

	t = NextToken(lexer);
	if (t.type != TOK_R_BRACKET)
		ReportError(lexer, "Expected ']'");

	return v;
}

static void NewAreaLight(Scene *scene, i32 obj_idx, Vec3 color)
{
	if (scene->light_count >= MAX_LIGHTS)
		Panic("Too much area lights");

	i32 light_idx = scene->light_count++;
	Light *light = &scene->lights[light_idx];
	light->color = color;
	//light->object_idx = obj_idx;

	Object *obj = &scene->objects[obj_idx];
	obj->light_idx = light_idx;
}

void LoadScene(Scene *scene, const char *file)
{
	scene->camera = {};
	scene->up = {0, 0, 1};
	scene->look_at = {};
	scene->defocus_angle = -1;
	scene->exposure = 1;
	scene->fov = 90;
	scene->width = 512;
	scene->height = 512;
	scene->samples = 20;
	scene->sky_box_color = {};

	scene->object_count = 0;
	scene->material_count = 1;
	scene->light_count = 0;

	scene->materials[0] = {
		{.5f, .5f, .5f},
		1, 1.5f, 0
	};

	FILE *f = fopen(file, "rb");
	if (f == nullptr)
		Panic("Couldn't open '%s'\n", file);

	fseek(f, 0, SEEK_END);
	isize length = ftell(f);
	rewind(f);
	char *data = (char *)malloc(length+1);
	fread(data, 1, length, f);
	fclose(f);

	data[length] = '\0';

	Lexer lexer[1];
	lexer->input = data;
	lexer->pos = 0;
	lexer->prev_pos = 0;
	lexer->file_name = file;

	i32 material_idx = 0;
	bool area_light = false;
	Vec3 area_light_color = {};

	String cmd;
	String key;
	while (ReadCmd(lexer, &cmd)) {
		if (StringEquals(cmd, S("sphere"))) {
			if (scene->object_count >= MAX_OBJECTS)
				Panic("Too much objects\n");

			i32 obj_idx = scene->object_count++;
			Object *obj = &scene->objects[obj_idx];
			obj->type = OBJ_SPHERE;
			obj->material_idx = material_idx;
			obj->light_idx = -1;
			obj->sphere = {};

#if 0
			if (area_light)
				NewAreaLight(scene, obj_idx, area_light_color);
#endif

			while (ReadKey(lexer, &key)) {
				if (StringEquals(key, S("radius"))) {
					obj->sphere.radius = ReadNumber(lexer);
				} else if (StringEquals(key, S("center"))) {
					obj->sphere.center = ReadVec3(lexer);
				} else {
					ReportError(lexer, "Unknown key");
				}
			}
		} else if (StringEquals(cmd, S("triangle_mesh"))) {
			while (ReadKey(lexer, &key)) {
				if (StringEquals(key, S("vertices"))) {

					BeginArray(lexer);

					while (!EndArray(lexer)) {
						if (scene->object_count >= MAX_OBJECTS)
							Panic("Too much objects");

						i32 obj_idx = scene->object_count++;
						Object *obj = &scene->objects[obj_idx];
						obj->type = OBJ_TRIANGLE;
						obj->material_idx = material_idx;
						obj->light_idx = -1;

						obj->triangle.v0.x = ReadNumber(lexer);
						obj->triangle.v0.y = ReadNumber(lexer);
						obj->triangle.v0.z = ReadNumber(lexer);
						obj->triangle.v1.x = ReadNumber(lexer);
						obj->triangle.v1.y = ReadNumber(lexer);
						obj->triangle.v1.z = ReadNumber(lexer);
						obj->triangle.v2.x = ReadNumber(lexer);
						obj->triangle.v2.y = ReadNumber(lexer);
						obj->triangle.v2.z = ReadNumber(lexer);

						Vec3 e0 = obj->triangle.v1 - obj->triangle.v0;
						Vec3 e1 = obj->triangle.v2 - obj->triangle.v0;
						Vec3 normal = Normalize(Cross(e0, e1));

						obj->triangle.n0 = normal;
						obj->triangle.n1 = normal;
						obj->triangle.n2 = normal;

						if (area_light)
							NewAreaLight(scene, obj_idx, area_light_color);
					}

				} else {
					ReportError(lexer, "Unknown key");
				}
			}
		} else if (StringEquals(cmd, S("material"))) {
			if (scene->material_count >= MAX_MATERIALS)
				Panic("Too much materials");

			material_idx = scene->material_count++;
			Material *mat = &scene->materials[material_idx];
			mat->color = V3(.7f);
			mat->roughness = 1;
			mat->ior = 1.5f;
			mat->metallic = 0;

			while (ReadKey(lexer, &key)) {
				if (StringEquals(key, S("color"))) {
					mat->color = ReadVec3(lexer);
				} else if (StringEquals(key, S("roughness"))) {
					float roughness = ReadNumber(lexer);
					mat->roughness = roughness*roughness;
				} else if (StringEquals(key, S("ior"))) {
					mat->ior = ReadNumber(lexer);
				} else if (StringEquals(key, S("metallic"))) {
					mat->metallic = ReadNumber(lexer);
				} else {
					ReportError(lexer, "Unknown key");
				}
			}
		} else if (StringEquals(cmd, S("area_light"))) {
			area_light = true;
			while (ReadKey(lexer, &key)) {
				if (StringEquals(key, S("color"))) {
					area_light_color = ReadVec3(lexer);
				} else {
					ReportError(lexer, "Unknown key");
				}
			}
		} else if (StringEquals(cmd, S("render"))) {
			while (ReadKey(lexer, &key)) {
				if (StringEquals(key, S("camera"))) {
					scene->camera = ReadVec3(lexer);
				} else if (StringEquals(key, S("look_at"))) {
					scene->look_at = ReadVec3(lexer);
				} else if (StringEquals(key, S("up"))) {
					scene->up = ReadVec3(lexer);
				} else if (StringEquals(key, S("fov"))) {
					scene->fov = ReadNumber(lexer);
				} else if (StringEquals(key, S("defocus_angle"))) {
					scene->defocus_angle = ReadNumber(lexer);
				} else if (StringEquals(key, S("exposure"))) {
					scene->exposure = ReadNumber(lexer);
				} else if (StringEquals(key, S("width"))) {
					scene->width = ReadI16(lexer);
				} else if (StringEquals(key, S("height"))) {
					scene->height = ReadI16(lexer);
				} else if (StringEquals(key, S("samples"))) {
					scene->samples = ReadI16(lexer);
				} else if (StringEquals(key, S("sky_box_color"))) {
					scene->sky_box_color = ReadVec3(lexer);
				} else {
					ReportError(lexer, "Unknown key");
				}
			}
		} else {
			ReportError(lexer, "Unknown command");
		}
	}
}
