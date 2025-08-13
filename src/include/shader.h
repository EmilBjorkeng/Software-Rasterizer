#ifndef __SHADER_H__
#define __SHADER_H__

#include "vector.h"
#include "light.h"
#include <cstdint>
#include <vector>

class Material;
struct Triangle;

class Shader {
protected:
    explicit Shader(bool transparent = false) : hasTransparency(transparent) {}
public:
    bool hasTransparency;

    virtual uint32_t pixelColor(const Material &material, float2 texCoord, float3 normal, const float3 &point,
        const std::vector<const Light*> &light, uint32_t &dstColor, std::vector<Triangle> &triangles) const = 0;
    virtual ~Shader() = default;
};

class DiffuseShader : public Shader {
public:
    uint32_t pixelColor(const Material &material, float2, float3, const float3&,
        const std::vector<const Light*>&, uint32_t&, std::vector<Triangle>&) const override;
};

class LitDiffuseShader : public Shader {
    static float ambientLight;
    static float3 ambientLightColor;
public:
    uint32_t pixelColor(const Material &material, float2, float3 normal, const float3 &point,
        const std::vector<const Light*> &light, uint32_t&, std::vector<Triangle> &triangles) const override;
};

class DiffuseAlphaShader : public Shader {
public:
    DiffuseAlphaShader() : Shader(true) {}
    uint32_t pixelColor(const Material &material, float2, float3, const float3&,
        const std::vector<const Light*>&, uint32_t &dstColor, std::vector<Triangle>&) const override;
};

class TextureShader : public Shader {
public:
    uint32_t pixelColor(const Material &material, float2 texCoord, float3, const float3&,
        const std::vector<const Light*>&, uint32_t&, std::vector<Triangle>&) const override;
};

class LitTextureShader : public Shader {
    static float ambientLight;
    static float3 ambientLightColor;
public:
    uint32_t pixelColor(const Material &material, float2 texCoord, float3 normal, const float3 &point,
        const std::vector<const Light*> &light, uint32_t&, std::vector<Triangle> &triangles) const override;
};

/*class LitSpecularShader : public Shader {
    static float ambientLight;
    static float3 ambientLightColor;
public:
    float4 pixelColor(const Material &material, float2, float3 normal, Light *light, const float3 &point, float3 cameraPosition) override;
};*/

#endif