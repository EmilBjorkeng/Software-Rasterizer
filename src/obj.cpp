#include "obj.h"
#include "vector.h"
#include <vector>
#include <fstream>
#include <sstream>
#include <algorithm>

ImageTexture::ImageTexture(const char* path) {
    data = stbi_load(path, &width, &height, &channels, 0);
    if (!data) {
        std::cerr << "Failed to load image: " << path << std::endl;
    }
}
ImageTexture::~ImageTexture() {
    if (data) stbi_image_free(data);
}
float4 ImageTexture::pixelColor(float2 position) {
    if (!data) return float4::Filled(-1);

    if (position.x() < 0.0f || position.x() >= 1.0f || position.y() < 0.0f || position.y() >= 1.0f)
        return float4::Filled(-1);

    int px = std::clamp(static_cast<int>(position.x() * width), 0, width - 1);
    int py = std::clamp(static_cast<int>((1.0f - position.y()) * height), 0, height - 1); // flip Y
    int index = (py * width + px) * channels;

    if (index + (channels - 1) >= width * height * channels) return float4::Filled(-1);

    float r = data[index + 0] / 255.0f;
    float g = (channels >= 2) ? data[index + 1] / 255.0f : r;
    float b = (channels >= 3) ? data[index + 2] / 255.0f : r;
    float a = (channels == 4) ? data[index + 3] / 255.0f : 1.0f;

    return float4(r, g, b, a);
}

Object LoadObjFile(const char* path) {
    std::ifstream in(path);
    if (!in.is_open()) {
        std::cerr << "Failed to open OBJ file: " << path << "\n";
        return {};
    }
    std::string line;

    std::vector<float3> vertices;
    std::vector<float2> textureCoordinates;
    std::vector<float3> normals;

    std::string currentMaterial;
    std::string materialLib;

    Object object;

    int lineNumber = 0;
    while (std::getline(in, line)) {
        ++lineNumber;

        std::istringstream ss(line);
        std::string type;
        ss >> type;

        if (type == "v") {
            float3 v;
            ss >> v.x() >> v.y() >> v.z();
            v.z() *= -1; // Flip vertex Z
            vertices.push_back(v);

            // Axis-Aligned Bounding Box
            object.aabbMin.x() = std::min(object.aabbMin.x(), v.x());
            object.aabbMin.y() = std::min(object.aabbMin.y(), v.y());
            object.aabbMin.z() = std::min(object.aabbMin.z(), v.z());

            object.aabbMax.x() = std::max(object.aabbMax.x(), v.x());
            object.aabbMax.y() = std::max(object.aabbMax.y(), v.y());
            object.aabbMax.z() = std::max(object.aabbMax.z(), v.z());
        }
        else if (type == "vt") {
            float2 t;
            ss >> t.x() >> t.y();
            textureCoordinates.push_back(t);
        }
        else if (type == "vn") {
            float3 n;
            ss >> n.x() >> n.y() >> n.z();
            n.z() *= -1; // Flip normal Z
            normals.push_back(n);
        }
        else if (type == "f") {
            std::vector<int> face_vertices;
            std::vector<int> face_texture;
            std::vector<int> face_normal;
            std::string vert;
            while (ss >> vert) {
                int v = 0, vt = 0, vn = 0;

                std::vector<std::string> parts;
                std::stringstream tokenStream(vert);
                std::string token;

                // Split by '/'
                while (std::getline(tokenStream, token, '/')) {
                    parts.push_back(token);
                }

                try {
                    if (parts.size() >= 1 && !parts[0].empty())
                        v = std::stoi(parts[0]);
                    if (parts.size() >= 2 && !parts[1].empty())
                        vt = std::stoi(parts[1]);
                    if (parts.size() >= 3 && !parts[2].empty())
                        vn = std::stoi(parts[2]);
                }
                catch (const std::exception& err) {
                    std::cerr << "Invalid face format: " << vert
                        << " (" << err.what() << ") at line " << lineNumber
                        << " in file: " << path << std::endl;
                    continue;
                }

                face_vertices.push_back(v);
                face_texture.push_back(vt);
                face_normal.push_back(vn);
            }

            int count = static_cast<int>(face_vertices.size());
            Face face;

            for (int i = 0; i < count; i++) {
                int v = face_vertices[i];
                int vt = face_texture[i];
                int vn = face_normal[i];

                if (v > 0 && v <= static_cast<int>(vertices.size())) {
                    face.Vertices.push_back(vertices[v - 1]);
                } else {
                    std::cerr << "Invalid vertex index: " << v
                            << " (vertex array size = " << vertices.size() << ")"
                            << " at line " << lineNumber
                            << " in file: " << path << std::endl;
                    face.Vertices.push_back({0, 0, 0});
                }
                if (vt > 0 && vt <= static_cast<int>(textureCoordinates.size())) {
                    face.VertexTextureCoords.push_back(textureCoordinates[vt - 1]);
                } else {
                    face.VertexTextureCoords.push_back({0, 0});
                }
                if (vn > 0 && vn <= static_cast<int>(normals.size())) {
                    face.VertexNormals.push_back(normals[vn - 1]);
                }
                else {
                    face.VertexNormals.push_back({0, 0, 0});
                }
            }

            // Get the MTL file path
            std::string objDir(path);
            size_t lastSlash = objDir.find_last_of("/\\");
            objDir = (lastSlash != std::string::npos) ? objDir.substr(0, lastSlash + 1) : "";
            std::string mtlPath = objDir + materialLib;

            // Check if material is already loaded
            for (auto& mat : object.materials) {
                if (mat.name == currentMaterial) {
                    face.material = std::make_shared<Material>(mat);
                    break;
                }
            }
            if (face.material == nullptr) {
                // Load material if not found
                Material newMat = loadMaterial(mtlPath.c_str(), currentMaterial);
                object.materials.push_back(newMat);
                face.material = std::make_shared<Material>(object.materials.back());
            }
            object.faces.push_back(face);
        }
        else if (type == "mtllib") {
            ss >> materialLib;
        }
        else if (type == "usemtl") {
            ss >> currentMaterial;
        }
    }

    return object;
}

Material loadMaterial(const char* path, std::string currentMaterial) {
    std::ifstream in(path);
    if (!in.is_open()) {
        std::cerr << "Failed to open MTL file: " << path << "\n";
        return {};
    }

    Material mat;

    // Default values
    mat.name = currentMaterial;
    mat.ambientColor = float3(0.1f, 0.1f, 0.1f);
    mat.diffuseColor = float3(0.8f, 0.8f, 0.8f);
    mat.specularColor = float3(0.0f, 0.0f, 0.0f);
    mat.EmissiveColor = float3(0.0f, 0.0f, 0.0f);
    mat.shininess = 0.0f;
    mat.refraction = 1.0f;
    mat.opacity = 1.0f;
    mat.illumModel = 2;
    mat.diffuseTexture = nullptr;
    mat.specularTexture = nullptr;
    mat.normalMap = nullptr;

    std::string line;
    bool skipping = false;
    while (std::getline(in, line)) {
        std::istringstream ss(line);
        std::string type;
        ss >> type;

        if (type == "newmtl") {
            skipping = false;
            std::string name;
            ss >> name;
            if (name != currentMaterial) {
                // Skip materials until we find the one we want
                skipping = true;
                continue;
            }
            mat.name = name;
        } else if (skipping) {
            continue;
        // Color
        } else if (type == "Ka") {
            ss >> mat.ambientColor.x() >> mat.ambientColor.y() >> mat.ambientColor.z();
        } else if (type == "Kd") {
            ss >> mat.diffuseColor.x() >> mat.diffuseColor.y() >> mat.diffuseColor.z();
        } else if (type == "Ks") {
            ss >> mat.specularColor.x() >> mat.specularColor.y() >> mat.specularColor.z();
        } else if (type == "Ke") {
            ss >> mat.EmissiveColor.x() >> mat.EmissiveColor.y() >> mat.EmissiveColor.z();
        // Properties
        } else if (type == "Ns") {
            ss >> mat.shininess;
        } else if (type == "Ni") {
            ss >> mat.refraction;
        } else if (type == "d") {
            ss >> mat.opacity;
        } else if (type == "Tr") {
            float  transparency;
            ss >> transparency;
            mat.opacity = 1.0f - transparency;
        // Texture
        } else if (type == "map_Kd") {
            ss >> mat.diffuseTexturePath;
            mat.diffuseTexture = new ImageTexture(mat.diffuseTexturePath.c_str());
            if (mat.diffuseTexture->data == nullptr) {
                delete mat.diffuseTexture;
                mat.diffuseTexture = nullptr;
            }
        } else if (type == "map_Ks") {
            ss >> mat.specularTexturePath;
            mat.specularTexture = new ImageTexture(mat.specularTexturePath.c_str());
            if (mat.specularTexture->data == nullptr) {
                delete mat.specularTexture;
                mat.specularTexture = nullptr;
            }
        } else if (type == "map_Bump" || type == "bump") {
            ss >> mat.normalMapPath;
            mat.normalMap = new ImageTexture(mat.normalMapPath.c_str());
            if (mat.normalMap->data == nullptr) {
                delete mat.normalMap;
                mat.normalMap = nullptr;
            }
        // Illumination
        } else if (type == "illum") {
            ss >> mat.illumModel;
        }
    }

    return mat;
}