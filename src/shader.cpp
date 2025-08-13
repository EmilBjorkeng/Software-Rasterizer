#include "shader.h"
#include "obj.h"
#include "render.h"
#include <algorithm>

#define AMBIENT_LIGHT 0.1f
#define AMBIENT_COLOR float3(1.0f, 1.0f, 1.0f);

// Ambient Light
float LitDiffuseShader::ambientLight = AMBIENT_LIGHT;
float3 LitDiffuseShader::ambientLightColor = AMBIENT_COLOR;

float LitTextureShader::ambientLight = AMBIENT_LIGHT;
float3 LitTextureShader::ambientLightColor = AMBIENT_COLOR;

auto unpack = [](uint32_t c, float &r, float &g, float &b, float &a) {
    a = (c & 0xFF) / 255.0f;
    b = ((c >> 8) & 0xFF) / 255.0f;
    g = ((c >> 16) & 0xFF) / 255.0f;
    r = ((c >> 24) & 0xFF) / 255.0f;
};
auto pack = [](float r, float g, float b, float a) -> uint32_t {
    uint32_t R = (uint32_t)(r * 255.0f) << 24;
    uint32_t G = (uint32_t)(g * 255.0f) << 16;
    uint32_t B = (uint32_t)(b * 255.0f) << 8;
    uint32_t A = (uint32_t)(a * 255.0f);
    return R | G | B | A;
};

// --------------------
// DiffuseShader
// --------------------

// === DiffuseShader ===
uint32_t DiffuseShader::pixelColor(const Material &material, float2, float3, const float3&,
    const std::vector<const Light*>&, uint32_t&, std::vector<Triangle>&) const {
    float r = std::clamp(material.diffuseColor.x(), 0.0f, 1.0f);
    float g = std::clamp(material.diffuseColor.y(), 0.0f, 1.0f);
    float b = std::clamp(material.diffuseColor.z(), 0.0f, 1.0f);
    return pack(r, g, b, 1.0f);
}

bool PointInTriangle(const float3 &P, const float3 &A, const float3 &B, const float3 &C) {
    float3 v0 = C - A;
    float3 v1 = B - A;
    float3 v2 = P - A;

    float d00 = v0.dot(v0);
    float d01 = v0.dot(v1);
    float d11 = v1.dot(v1);
    float d20 = v2.dot(v0);
    float d21 = v2.dot(v1);

    // Compute barycentric coordinates
    float denom = d00 * d11 - d01 * d01;
    if (denom == 0.0f) return false; // Degenerate Triangle

    float v = (d11 * d20 - d01 * d21) / denom;
    float w = (d00 * d21 - d01 * d20) / denom;
    float u = 1.0f - v - w;

    // Check if point is inside the triangle
    return (u >= 0) && (v >= 0) && (w >= 0);
}

// === LitDiffuseShader ===
uint32_t LitDiffuseShader::pixelColor(const Material &material, float2, float3 normal, const float3 &point,
    const std::vector<const Light*> &light, uint32_t&, std::vector<Triangle> &triangles) const {
    float3 diffuse = float3::Zero();

    // Diffuse lighting contribution
    for (const Light* l : light) {
        float3 lightPosition = l->getPosition(point);
        float3 DirectionToLight = l->lightDirection(point);

        float ndotl = normal.dot(DirectionToLight);
        if (ndotl < 1e-6f) continue; // Skip if there is no light hitting the pixel

        float d = l->distance(point);

        // Bounding box between point and light
        float minX = std::min(point.x(), lightPosition.x());
        float maxX = std::max(point.x(), lightPosition.x());
        float minY = std::min(point.y(), lightPosition.y());
        float maxY = std::max(point.y(), lightPosition.y());
        float minZ = std::min(point.z(), lightPosition.z());
        float maxZ = std::max(point.z(), lightPosition.z());

        // Offset origin to avoid self-intersection
        float3 origin = point + normal * 1e-4f;

        // Checks if point is in a shadow
        bool inShadow = false;
        for (int i = 0; i < static_cast<int>(triangles.size()); i++) {
            const float3 &a = triangles[i].a;
            const float3 &b = triangles[i].b;
            const float3 &c = triangles[i].c;
            const float3 &N = triangles[i].normal;

            // Axis-Aligned Bounding Box
            const float3 aabbMin = triangles[i].aabbMin;
            const float3 aabbMax = triangles[i].aabbMax;

            // Skip if all triangles verticies are outside the bounding box between point and light
            if (aabbMax.x() < minX || aabbMin.x() > maxX ||
                aabbMax.y() < minY || aabbMin.y() > maxY ||
                aabbMax.z() < minZ || aabbMin.z() > maxZ) {
                continue;
            }

            // Get the amout of DirectionToLight (vector) it takes to go from point to the plane
            // Skip if it's negative or above 1 (further away then the light)
            float denom = DirectionToLight.dot(N);
            if (denom > -1e-6f) continue; // Parallel and backface culling
            float t = (a - origin).dot(N) / denom;
            if (t <= 0 || t >= d) continue;

            // Check if the intersection point between the vector (DirectionToLight) and the plane is inside the triangle (abc)
            float3 intersectionPoint = origin + (DirectionToLight * t);
            if (PointInTriangle(intersectionPoint, a, b, c)) {
                inShadow = true;
                break;
            }
        }

        // Make the light dimmer the further it is away (inverse square)
        if (d < 1e-6f) d = 1e-6f; // Avoid divide by zero
        float invDistance = 1 / (d * d);

        if (!inShadow)
            diffuse += material.diffuseColor * l->color * l->intensity * ndotl * invDistance;
    }

    // Ambient lighting contribution
    float3 ambient = material.diffuseColor * ambientLightColor * ambientLight;

    float3 color = ambient + diffuse;

    // Clamp each component to [0,1]
    float r = std::clamp(color.x(), 0.0f, 1.0f);
    float g = std::clamp(color.y(), 0.0f, 1.0f);
    float b = std::clamp(color.z(), 0.0f, 1.0f);

    return pack(r, g, b, 1.0f);
}

// === DiffuseAlphaShader ===
uint32_t DiffuseAlphaShader::pixelColor(const Material &material, float2, float3, const float3&,
    const std::vector<const Light*>&, uint32_t &dstColor, std::vector<Triangle>&) const {
    float r = std::clamp(material.diffuseColor.x(), 0.0f, 1.0f);
    float g = std::clamp(material.diffuseColor.y(), 0.0f, 1.0f);
    float b = std::clamp(material.diffuseColor.z(), 0.0f, 1.0f);
    float a = std::clamp(material.opacity, 0.0f, 1.0f);

    // Skip if fully transparent
    if (a < 1e-4f)
    { return 0; }

    if (a < 0.999f) {
        // Partially transparent — blend
        float dr, dg, db, da;
        unpack(dstColor, dr, dg, db, da);

        float outA = a + da * (1.0f - a);
        if (outA < 1e-6f) outA = 1e-6f; // Avoid divide by zero

        float outR = (r * a + dr * da * (1.0f - a)) / outA;
        float outG = (g * a + dg * da * (1.0f - a)) / outA;
        float outB = (b * a + db * da * (1.0f - a)) / outA;

        return pack(outR, outG, outB, outA);
    }

    // Fully opaque
    return pack(r, g, b, 1.0f);
}

// --------------------
// TextureShader
// --------------------

// === TextureShader ===
uint32_t TextureShader::pixelColor(const Material &material, float2 texCoord, float3, const float3&,
    const std::vector<const Light*>&, uint32_t&, std::vector<Triangle>&) const {
    float3 baseColor = material.diffuseColor;

    float3 textureColor = float3::Zero();
    if (material.diffuseTexture) {
        textureColor = material.diffuseTexture->pixelColor(texCoord);
    }
    float3 combinedColor = baseColor * (textureColor + float3::Filled(1e-6f)); // prevent black texture nullifying everything

    return pack(combinedColor.x(), combinedColor.y(), combinedColor.z(), 1.0f);
}

// === LitTextureShader ===
uint32_t LitTextureShader::pixelColor(const Material &material, float2 texCoord, float3 normal, const float3 &point,
    const std::vector<const Light*> &light, uint32_t&, std::vector<Triangle> &triangles) const {
    float3 baseColor = material.diffuseColor;

    float3 textureColor = float3::Zero();
    if (material.diffuseTexture) {
        textureColor = material.diffuseTexture->pixelColor(texCoord);
    }
    float3 combinedColor = baseColor * (textureColor + float3::Filled(1e-6f)); // prevent black texture nullifying everything

    float3 diffuse = float3::Zero();

    // Diffuse lighting contribution
    for (const Light* l : light) {
        float3 lightPosition = l->getPosition(point);
        float3 DirectionToLight = l->lightDirection(point);

        float ndotl = normal.dot(DirectionToLight);
        if (ndotl < 1e-6f) continue; // Skip if there is no light hitting the pixel

        float d = l->distance(point);

        // Bounding box between point and light
        float minX = std::min(point.x(), lightPosition.x());
        float maxX = std::max(point.x(), lightPosition.x());
        float minY = std::min(point.y(), lightPosition.y());
        float maxY = std::max(point.y(), lightPosition.y());
        float minZ = std::min(point.z(), lightPosition.z());
        float maxZ = std::max(point.z(), lightPosition.z());

        // Offset origin to avoid self-intersection
        float3 origin = point + normal * 1e-4f;

        // Checks if point is in a shadow
        bool inShadow = false;
        for (int i = 0; i < static_cast<int>(triangles.size()); i++) {
            const float3 &a = triangles[i].a;
            const float3 &b = triangles[i].b;
            const float3 &c = triangles[i].c;
            const float3 &N = triangles[i].normal;

            // Axis-Aligned Bounding Box
            const float3 aabbMin = triangles[i].aabbMin;
            const float3 aabbMax = triangles[i].aabbMax;

            // Skip if all triangles verticies are outside the bounding box between point and light
            if (aabbMax.x() < minX || aabbMin.x() > maxX ||
                aabbMax.y() < minY || aabbMin.y() > maxY ||
                aabbMax.z() < minZ || aabbMin.z() > maxZ) {
                continue;
            }

            // Get the amout of DirectionToLight (vector) it takes to go from point to the plane
            // Skip if it's negative or above 1 (further away then the light)
            float denom = DirectionToLight.dot(N);
            if (denom > -1e-6f) continue; // Parallel and backface culling
            float t = (a - origin).dot(N) / denom;
            if (t <= 0 || t >= d) continue;

            // Check if the intersection point between the vector (DirectionToLight) and the plane is inside the triangle (abc)
            float3 intersectionPoint = origin + (DirectionToLight * t);
            if (PointInTriangle(intersectionPoint, a, b, c)) {
                inShadow = true;
                break;
            }
        }

        // Make the light dimmer the further it is away (inverse square)
        if (d < 1e-6f) d = 1e-6f; // Avoid divide by zero
        float invDistance = 1 / (d * d);

        if (!inShadow)
            diffuse += combinedColor * l->color * l->intensity * ndotl * invDistance;
    }

    // Ambient lighting contribution
    float3 ambient = combinedColor * ambientLightColor * ambientLight;

    float3 color = ambient + diffuse;

    float r = std::clamp(color.x(), 0.0f, 1.0f);
    float g = std::clamp(color.y(), 0.0f, 1.0f);
    float b = std::clamp(color.z(), 0.0f, 1.0f);
    return pack(r, g, b, 1.0f);
}