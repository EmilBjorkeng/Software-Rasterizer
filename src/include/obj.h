#ifndef __OBJ_H__
#define __OBJ_H__

#include "vector.h"
#include "shader.h"
#include "STB/stb_image.h"
#include <vector>
#include <memory>
#include <string>
#include <cstdint>

struct ImageTexture {
    int width, height, channels;
    unsigned char* data;

    ImageTexture(const char* path);
    ~ImageTexture();

    float4 pixelColor(float2 position);
};

struct Material {
    std::string name;

    float3 ambientColor;    // Ka - Ambient color
    float3 diffuseColor;    // Kd - Diffuse color
    float3 specularColor;   // Ks - Specular color
    float3 EmissiveColor;   // Ke - Emissive color
    float shininess;        // Ns - Specular exponent
    float refraction;       // Ni - Optical density
    float opacity;          // d or Tr - Transparency (1.0 = opaque)

    std::string diffuseTexturePath;     // map_Kd - Diffuse texture
    std::string specularTexturePath;    // map_Ks - Specular texture
    std::string normalMapPath;          // bump or map_Bump - Normal map

    ImageTexture *diffuseTexture = nullptr;
    ImageTexture *specularTexture = nullptr;
    ImageTexture *normalMap = nullptr;

    int illumModel; // illum - Illumination model
};

struct Face {
    std::vector<float3> Vertices;               // Vertex Position
    std::vector<float2> VertexTextureCoords;    // Vertex Texture Coordinates
    std::vector<float3> VertexNormals;          // Vertex Normals
    std::shared_ptr<Material> material;
};

struct Object {
    std::vector<Face> faces;
    std::vector<Material> materials;
    float3 position = float3::Zero();
    float3 rotation = float3::Zero();
    float3 scale = float3::One();
    std::unique_ptr<Shader> shader;

    float3 aabbMin = float3::Filled(-INFINITY);
    float3 aabbMax = float3::Filled(INFINITY);


    Object() :
        shader(std::make_unique<DiffuseShader>())
    {}
};

Object LoadObjFile(const char* path);
Material loadMaterial(const char* path, std::string Material);

#endif