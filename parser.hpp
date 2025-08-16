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
	char *input;
	isize pos;
};

bool IsSpace(char c)
{
	return c == ' ' || c == '\t' || c == '\r' || c == '\n';
}

bool IsAlpha(char c)
{
	switch (c) {
	case '\0': case '{': case '}': case '[': case ']': case '#':
		return false;
	default:
		return !IsSpace(c);
	}
}

bool StringEquals(String a, String b)
{
	if (a.length != b.length)
		return false;
	return memcmp(a.data, b.data, a.length) == 0;
}

Token NextToken(Lexer *lexer)
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

	return {type, {&lexer->input[start], lexer->pos - start}};
}

bool ReadCmd(Lexer *lexer, String *cmd)
{
	Token t = NextToken(lexer);

	if (t.type == TOK_END) {
		return false;
	} else if (t.type == TOK_STRING) {
		if (NextToken(lexer).type != TOK_L_BRACE)
			Panic("Expected '{'\n");

		*cmd = t.lexeme;
		return true;
	}

	Panic("Expected command\n");
	return false;
}

bool ReadKey(Lexer *lexer, String *key)
{
	Token t = NextToken(lexer);

	if (t.type == TOK_R_BRACE) {
		return false;
	} else if (t.type == TOK_STRING) {
		*key = t.lexeme;
		return true;
	}

	Panic("Expected key\n");
	return false;
}

String ReadString(Lexer *lexer)
{
	Token t = NextToken(lexer);

	if (t.type == TOK_STRING) {
		return t.lexeme;
	}

	Panic("Expected string\n");
	return {};
}

int StringToNumber(String s, float *n)
{
	char *end = nullptr;
	errno = 0;
	*n = strtof(s.data, &end);
	if (end - s.data != s.length || errno)
		return -1;
	return 0;
}

float ReadNumber(Lexer *lexer)
{
	String s = ReadString(lexer);
	float n;
	if (StringToNumber(s, &n) < 0)
		Panic("Invalid number\n");
	return n;
}

i16 ReadI16(Lexer *lexer)
{
	String s = ReadString(lexer);
	long n;
	char *end = nullptr;
	errno = 0;
	n = strtol(s.data, &end, 10);
	if (end - s.data != s.length || errno)
		Panic("Invalid integer\n");
	if (n > (1 << 16) - 1 || n < 0)
		Panic("Out of bounds integer\n");
	return (i16)n;
}

void BeginArray(Lexer *lexer)
{
	Token t = NextToken(lexer);

	if (t.type != TOK_L_BRACKET)
		Panic("Expected array\n");
}

bool EndArray(Lexer *lexer)
{
	isize start = lexer->pos;
	Token t = NextToken(lexer);

	if (t.type == TOK_R_BRACKET)
		return true;

	lexer->pos = start;
	return false;
}

Vec3 ReadVec3(Lexer *lexer)
{
	Vec3 v;

	if (NextToken(lexer).type != TOK_L_BRACKET)
		Panic("Expected array\n");

	v.x = ReadNumber(lexer);
	v.y = ReadNumber(lexer);
	v.z = ReadNumber(lexer);

	if (NextToken(lexer).type != TOK_R_BRACKET)
		Panic("Expected ']'\n");

	return v;
}

void NewAreaLight(Scene *scene, i32 obj_idx, Vec3 color)
{
	if (scene->light_count >= MAX_LIGHTS)
		Panic("Too much area lights");

	i32 light_idx = scene->light_count++;
	Light *light = &scene->lights[light_idx];
	light->color = color;
	light->object_idx = obj_idx;

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
					Panic("Unknown key '%.*s'\n", (int)key.length, key.data);
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
					Panic("Unknown key '%.*s'\n", (int)key.length, key.data);
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
					mat->roughness = ReadNumber(lexer);
				} else if (StringEquals(key, S("ior"))) {
					mat->ior = ReadNumber(lexer);
				} else if (StringEquals(key, S("metallic"))) {
					mat->metallic = ReadNumber(lexer);
				} else {
					Panic("Unknown key '%.*s'\n", (int)key.length, key.data);
				}
			}
		} else if (StringEquals(cmd, S("area_light"))) {
			area_light = true;
			while (ReadKey(lexer, &key)) {
				if (StringEquals(key, S("color"))) {
					area_light_color = ReadVec3(lexer);
				} else {
					Panic("Unknown key '%.*s'\n", (int)key.length, key.data);
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
				} else {
					Panic("Unknown key '%.*s'\n", (int)key.length, key.data);
				}
			}
		} else {
			Panic("Unknown command '%.*s'\n", (int)cmd.length, cmd.data);
		}
	}
}
