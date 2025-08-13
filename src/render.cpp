#include "render.h"
#include "matrix.h"
#include "matrix_vector_bridge.h"
#include "shader.h"
#include <algorithm>

// 2D cross product (ap x ab)
float Cross2D(const float2 a, const float2 b, const float2 p) {
    return (p.x() - a.x()) * (b.y() - a.y()) - (p.y() - a.y()) * (b.x() - a.x());
}

void RendersceneObjects(uint32_t *drawBuffer, int pitch, float *depthBuffer, const Camera<float> &camera,
const std::vector<const Object*> &sceneObjects, const std::vector<const Light*> &light) {
    std::vector<RasterizerPoint> rasterizerPoints;
    std::vector<RasterizerPoint> rasterizerPointsAlpha;
    std::vector<Triangle> triangles;

    for (const Object *obj : sceneObjects) {
        // Precompute Matrices
        mat4<float> model_matrix =
            mat4<float>::Translate(obj->position) *
            mat4<float>::RotateZ(obj->rotation.z()) *
            mat4<float>::RotateX(obj->rotation.x()) *
            mat4<float>::RotateY(obj->rotation.y()) *
            mat4<float>::Scale(obj->scale);
        mat4<float> view_matrix = camera.GetViewMatrix();
        mat3<float> normal_matrix = mat3<float>(model_matrix).inverse().transpose();

        bool shaderTransparency = obj->shader->hasTransparency;

        for (int f = 0; f < static_cast<int>(obj->faces.size()); ++f) {
            const Face& face = obj->faces[f];

            float faceOpacity = face.material.get()->opacity;
            bool hasTrancparency = shaderTransparency && (faceOpacity < 1.0f);

            // Skip if face is transparent
            if (faceOpacity < 1e-4f)
            { continue; }

            // Convert face to triangles
            std::vector<int> VertexOrder = TriangulateEarClipping(face.Vertices);

            // Local To World
            std::vector<float3> worldSpaceVertices(face.Vertices.size());
            for (size_t i = 0; i < face.Vertices.size(); ++i) {
                worldSpaceVertices[i] = (float3)(model_matrix * float4(face.Vertices[i], 1.0f));
            }

            // Make a list of all triangles for lit shaders
            for (int i = 0; i < static_cast<int>(VertexOrder.size()); i += 3) {
                const float3 a = worldSpaceVertices[VertexOrder[i]];
                const float3 b = worldSpaceVertices[VertexOrder[i+1]];
                const float3 c = worldSpaceVertices[VertexOrder[i+2]];

                Triangle tri;
                tri.a = a; tri.b = b; tri.c = c;

                float3 ab = b - a;
                float3 ac = c - a;
                tri.normal = ab.cross(ac).normalized();

                // Axis-Aligned Bounding Box
                tri.aabbMin = float3(
                    std::min({a.x(), b.x(), c.x()}),
                    std::min({a.y(), b.y(), c.y()}),
                    std::min({a.z(), b.z(), c.z()})
                );
                tri.aabbMax = float3(
                    std::max({a.x(), b.x(), c.x()}),
                    std::max({a.y(), b.y(), c.y()}),
                    std::max({a.z(), b.z(), c.z()})
                );

                triangles.push_back(tri);
            }

            // World To View
            std::vector<float3> viewSpaceVertices(face.Vertices.size());
            for (size_t i = 0; i < face.Vertices.size(); ++i) {
                viewSpaceVertices[i] = (float3)(view_matrix * float4(worldSpaceVertices[i], 1.0f));
            }

            // Normal to World
            std::vector<float3> worldSpaceNormals(face.VertexNormals.size());
            for (size_t i = 0; i < face.VertexNormals.size(); ++i) {
                worldSpaceNormals[i] = normal_matrix * face.VertexNormals[i];
                worldSpaceNormals[i] = worldSpaceNormals[i].normalized();
            }

            // Slice triangles partially behind the camera
            // And Separate by transparency
            if (!hasTrancparency) {
                triangleSlicer(&rasterizerPoints, worldSpaceVertices, viewSpaceVertices, worldSpaceNormals, VertexOrder, camera.fov, face, obj);
            } else {
                triangleSlicer(&rasterizerPointsAlpha, worldSpaceVertices, viewSpaceVertices, worldSpaceNormals, VertexOrder, camera.fov, face, obj);
            }
        }
    }
    // Add the alpha points to the end
    rasterizerPoints.insert( rasterizerPoints.end(), rasterizerPointsAlpha.begin(), rasterizerPointsAlpha.end() );

    for (int i = 0; i < static_cast<int>(rasterizerPoints.size()); i += 3) {
        const RasterizerPoint &a = rasterizerPoints[i];
        const RasterizerPoint &b = rasterizerPoints[i+1];
        const RasterizerPoint &c = rasterizerPoints[i+2];
        if (a.depth <= 0 && b.depth <= 0 && c.depth <= 0) continue; // Skip if triangle is behind the camera

        // 2D Vertex points
        float2 Vertex0 = a.screenPos;
        float2 Vertex1 = b.screenPos;
        float2 Vertex2 = c.screenPos;

        float area = Cross2D(Vertex0, Vertex1, Vertex2);
        if (area == 0) continue; // Skip if area is 0

        // Precalculate some depth numbers
        float3 depths = float3(a.depth, b.depth, c.depth);
        float3 invDepths = 1.0f / depths;

        // Bound of the triangle
        int minX = std::clamp((int)std::floor( std::min(std::min(Vertex0.x(), Vertex1.x()), Vertex2.x()) ), 0, SCREEN_WIDTH);
        int maxX = std::clamp((int)std::ceil( std::max(std::max(Vertex0.x(), Vertex1.x()), Vertex2.x()) ), 0, SCREEN_WIDTH);
        int minY = std::clamp((int)std::floor( std::min(std::min(Vertex0.y(), Vertex1.y()), Vertex2.y()) ), 0, SCREEN_HEIGHT);
        int maxY = std::clamp((int)std::ceil( std::max(std::max(Vertex0.y(), Vertex1.y()), Vertex2.y()) ), 0, SCREEN_HEIGHT);

        // Top-left rule bias: Ensures shared edges are only drawn by one triangle
        // Applies -1 bias to pixels exactly on the right or bottom edges to exclude them
        constexpr float edgeBias = -1e-4f;
        float Bias0 = (Vertex1.y() > Vertex0.y()) || (Vertex1.y() == Vertex0.y() && Vertex1.x() > Vertex0.x()) ? 0.0f : edgeBias;
        float Bias1 = (Vertex2.y() > Vertex1.y()) || (Vertex2.y() == Vertex1.y() && Vertex2.x() > Vertex1.x()) ? 0.0f : edgeBias;
        float Bias2 = (Vertex0.y() > Vertex2.y()) || (Vertex0.y() == Vertex2.y() && Vertex0.x() > Vertex2.x()) ? 0.0f : edgeBias;

        for (int y = minY; y < maxY; ++y) {
            for (int x = minX; x < maxX; ++x) {
                float2 p(x + 0.5f, y + 0.5f); // Pixel center

                float Weight0 = Cross2D(Vertex1, Vertex2, p) + Bias0;
                float Weight1 = Cross2D(Vertex2, Vertex0, p) + Bias1;
                float Weight2 = Cross2D(Vertex0, Vertex1, p) + Bias2;

                // Check if pixel is in triangle
                if (Weight0 >= 0 && Weight1 >= 0 && Weight2 >= 0) {
                    float weightSum = Weight0 + Weight1 + Weight2;
                    float3 baryCoord(
                        Weight0 / weightSum,
                        Weight1 / weightSum,
                        Weight2 / weightSum
                    );
                    float depth = 1 / (invDepths.dot(baryCoord));

                    // Skip if something nearer has already been drawn
                    if (depth > depthBuffer[(y * SCREEN_WIDTH) + x]) continue;

                    //const Material &material = *face.material;
                    const Material &material = *a.vertex_face->material;

                    // Interpolate Texture
                    float2 texCoord{};
                    texCoord += a.vertexTextureCoords * invDepths[0] * baryCoord[0];
                    texCoord += b.vertexTextureCoords * invDepths[1] * baryCoord[1];
                    texCoord += c.vertexTextureCoords * invDepths[2] * baryCoord[2];
                    texCoord *= depth;

                    // Interpolate Normal
                    float3 normal{};
                    normal += a.vertexNormals * invDepths[0] * baryCoord[0];
                    normal += b.vertexNormals * invDepths[1] * baryCoord[1];
                    normal += c.vertexNormals * invDepths[2] * baryCoord[2];
                    normal *= depth;

                    // Interpolate World Space Point
                    float3 point{};
                    point += a.worldSpaceVertex * invDepths[0] * baryCoord[0];
                    point += b.worldSpaceVertex * invDepths[1] * baryCoord[1];
                    point += c.worldSpaceVertex * invDepths[2] * baryCoord[2];
                    point *= depth;

                    //uint32_t &dstColor = drawBuffer[(y * SCREEN_WIDTH) + x];
                    uint32_t &dstColor = drawBuffer[y * (pitch / sizeof(uint32_t)) + x];
                    //uint32_t srcColor = object.shader->pixelColor(material, texCoord, normal, point, light, dstColor, triangles);
                    uint32_t srcColor = a.object->shader->pixelColor(material, texCoord, normal, point, light, dstColor, triangles);

                    // Only write color & depth if there is anything to draw
                    if (srcColor > 0) {
                        dstColor = srcColor;
                        depthBuffer[(y * SCREEN_WIDTH) + x] = depth;
                    }
                }
            }
        }
    }
}

float2 ViewToScreen(const float3 vertex_view, int fov) {
    float fovRadians = fov * (M_PI / 180.0f);
    float screenHeight_world = tan(fovRadians / 2.0f) * 2.0f;
    float pixelsPerWorldUnit = SCREEN_HEIGHT / screenHeight_world / vertex_view.z();

    float2 pixelOffset = float2(vertex_view.x(), vertex_view.y()) * pixelsPerWorldUnit;
    float2 vertex_screen = float2(SCREEN_WIDTH, SCREEN_HEIGHT) / 2 + pixelOffset;
    return vertex_screen;
}

// Project vertices to 2D plane coordinates
std::vector<float2> ProjectTo2D(const std::vector<float3>& vertices) {
    std::vector<float2> projected;

    // Compute polygon normal (using first 3 vertices)
    float3 edge1 = vertices[1] - vertices[0];
    float3 edge2 = vertices[2] - vertices[0];
    float3 normal = edge1.cross(edge2).normalized();

    // Create orthonormal basis (u,v) on polygon plane
    float3 u = edge1.normalized();
    float3 v = normal.cross(u);  // guaranteed perpendicular to u and normal

    projected.reserve(vertices.size());

    for (const auto& p : vertices) {
        float3 vec = p - vertices[0];
        projected.push_back(
            float2(vec.dot(u), vec.dot(v))
        );
    }
    return projected;
}

bool IsConvex(const float3& prev, const float3& curr, const float3& next) {
    float dx1 = curr.x() - prev.x();
    float dy1 = curr.y() - prev.y();
    float dx2 = next.x() - curr.x();
    float dy2 = next.y() - curr.y();
    float cross = dx1 * dy2 - dy1 * dx2;
    return cross > 0; // Positive means convex
}

// Ear Clipping Algorithm
std::vector<int> TriangulateEarClipping(const std::vector<float3>& vertices) {
    // project to 2D
    std::vector<float2> vertices2D = ProjectTo2D(vertices);

    std::vector<int> VertexOrder;
    int n = static_cast<int>(vertices2D.size());
    if (n < 3) return VertexOrder;

    std::vector<int> indices(n);
    for (int i = 0; i < n; ++i) indices[i] = i;

    while (indices.size() > 3) {
        bool earFound = false;

        for (int i = 0; i < static_cast<int>(indices.size()); ++i) {
            int prev = indices[(i + indices.size() - 1) % indices.size()];
            int curr = indices[i];
            int next = indices[(i + 1) % indices.size()];

            const float3& a = vertices2D[prev];
            const float3& b = vertices2D[curr];
            const float3& c = vertices2D[next];

            if (!IsConvex(a, b, c)) continue;

            // Check no other vertex is inside triangle
            bool hasPointInside = false;
            for (int j = 0; j < static_cast<int>(indices.size()); ++j) {
                int vi = indices[j];
                if (vi == prev || vi == curr || vi == next) continue;

                float3 p = vertices2D[vi];
                bool sameSide =
                    (Cross2D(p, a, b) >= 0.0f) ==
                    (Cross2D(p, b, c) >= 0.0f) &&
                    (Cross2D(p, b, c) >= 0.0f) ==
                    (Cross2D(p, c, a) >= 0.0f);

                if (sameSide) {
                    hasPointInside = true;
                    break;
                }
            }

            if (hasPointInside) continue;

            // Found an ear
            VertexOrder.push_back(prev);
            VertexOrder.push_back(curr);
            VertexOrder.push_back(next);
            indices.erase(indices.begin() + i);
            earFound = true;
            break;
        }

        if (!earFound) {
            break;
        }
    }

    // Final triangle
    if (indices.size() == 3) {
        VertexOrder.push_back(indices[0]);
        VertexOrder.push_back(indices[1]);
        VertexOrder.push_back(indices[2]);
    }

    return VertexOrder;
}

void triangleSlicer(std::vector<RasterizerPoint> *rasterizerPoints,
    const std::vector<float3> &worldSpaceVertices, const std::vector<float3> &viewSpaceVertices,
    const std::vector<float3> &worldSpaceNormals, const std::vector<int> &vertexOrder, int fov,
    const Face &face, const Object *obj) {

    int size = static_cast<int>(vertexOrder.size());

    for (int i = 0; i < size; i += 3) {
        const float3 viewPoints[3] = {
            viewSpaceVertices[vertexOrder[i]],
            viewSpaceVertices[vertexOrder[i+1]],
            viewSpaceVertices[vertexOrder[i+2]]
        };
        const float3 worldPoints[3] = {
            worldSpaceVertices[vertexOrder[i]],
            worldSpaceVertices[vertexOrder[i+1]],
            worldSpaceVertices[vertexOrder[i+2]]
        };
        const float2 vertexTextureCoords[3] = {
            face.VertexTextureCoords[vertexOrder[i]],
            face.VertexTextureCoords[vertexOrder[i+1]],
            face.VertexTextureCoords[vertexOrder[i+2]]
        };
        const float3 vertexNormals[3] = {
            worldSpaceNormals[vertexOrder[i]],
            worldSpaceNormals[vertexOrder[i+1]],
            worldSpaceNormals[vertexOrder[i+2]]
        };

        constexpr float NEAR_CLIP = 0.01f;
        bool clip0 = viewPoints[0].z() <= NEAR_CLIP;
        bool clip1 = viewPoints[1].z() <= NEAR_CLIP;
        bool clip2 = viewPoints[2].z() <= NEAR_CLIP;
        int clipCount = (int)clip0 + (int)clip1 + (int)clip2;

        if (clipCount == 3) continue; // Skip if entire triangle is behind the camera

        switch (clipCount) {
            case 0:
                rasterizerPoints->push_back(RasterizerPoint(viewPoints[0].z(), ViewToScreen(viewPoints[0], fov), worldPoints[0], vertexTextureCoords[0], vertexNormals[0], obj, face));
                rasterizerPoints->push_back(RasterizerPoint(viewPoints[1].z(), ViewToScreen(viewPoints[1], fov), worldPoints[1], vertexTextureCoords[1], vertexNormals[1], obj, face));
                rasterizerPoints->push_back(RasterizerPoint(viewPoints[2].z(), ViewToScreen(viewPoints[2], fov), worldPoints[2], vertexTextureCoords[2], vertexNormals[2], obj, face));
                break;
            case 1: {
                // Find the point that needs to be sliced
                int clipIndex = clip0 ? 0 : clip1 ? 1 : 2;
                int nextIndex = (clipIndex + 1) % 3;
                int prevIndex = (clipIndex - 1 + 3) % 3;

                float3 pointClipped = viewPoints[clipIndex];
                float3 pointA = viewPoints[nextIndex];
                float3 pointB = viewPoints[prevIndex];

                // Find fractional value of the slice
                float fracA = (NEAR_CLIP - pointClipped.z()) / (pointA.z() - pointClipped.z());
                float fracB = (NEAR_CLIP - pointClipped.z()) / (pointB.z() - pointClipped.z());

                // New triangle points
                float3 clipPointAlongEdgeA = pointClipped.lerp(pointA, fracA);
                float3 clipPointAlongEdgeB = pointClipped.lerp(pointB, fracB);

                // New world points
                float3 worldClipPointAlongEdgeA = worldPoints[clipIndex].lerp(worldPoints[nextIndex], fracA);
                float3 worldClipPointAlongEdgeB = worldPoints[clipIndex].lerp(worldPoints[prevIndex], fracB);
                // New texture coord points
                float2 textureClipPointAlongEdgeA = vertexTextureCoords[clipIndex].lerp(vertexTextureCoords[nextIndex], fracA);
                float2 textureClipPointAlongEdgeB = vertexTextureCoords[clipIndex].lerp(vertexTextureCoords[prevIndex], fracB);
                // New normal points
                float3 normalClipPointAlongEdgeA = vertexNormals[clipIndex].lerp(vertexNormals[nextIndex], fracA);
                float3 normalClipPointAlongEdgeB = vertexNormals[clipIndex].lerp(vertexNormals[prevIndex], fracB);

                // Create the new triangles
                rasterizerPoints->push_back(RasterizerPoint(clipPointAlongEdgeB.z(), ViewToScreen(clipPointAlongEdgeB, fov), worldClipPointAlongEdgeB, textureClipPointAlongEdgeB,     normalClipPointAlongEdgeB, obj, face));
                rasterizerPoints->push_back(RasterizerPoint(clipPointAlongEdgeA.z(), ViewToScreen(clipPointAlongEdgeA, fov), worldClipPointAlongEdgeA, textureClipPointAlongEdgeA,     normalClipPointAlongEdgeA, obj, face));
                rasterizerPoints->push_back(RasterizerPoint(pointB.z(),              ViewToScreen(pointB, fov),              worldPoints[prevIndex],   vertexTextureCoords[prevIndex], vertexNormals[prevIndex],  obj, face));

                rasterizerPoints->push_back(RasterizerPoint(clipPointAlongEdgeA.z(), ViewToScreen(clipPointAlongEdgeA, fov), worldClipPointAlongEdgeA, textureClipPointAlongEdgeA,     normalClipPointAlongEdgeA, obj, face));
                rasterizerPoints->push_back(RasterizerPoint(pointA.z(),              ViewToScreen(pointA, fov),              worldPoints[nextIndex],   vertexTextureCoords[nextIndex], vertexNormals[nextIndex],  obj, face));
                rasterizerPoints->push_back(RasterizerPoint(pointB.z(),              ViewToScreen(pointB, fov),              worldPoints[prevIndex],   vertexTextureCoords[prevIndex], vertexNormals[prevIndex],  obj, face));
                break;
            }
            case 2: {
                // Find the point that does not needs to be sliced
                int nonClipIndex = !clip0 ? 0 : !clip1 ? 1 : 2;
                int nextIndex = (nonClipIndex + 1) % 3;
                int prevIndex = (nonClipIndex - 1 + 3) % 3;

                float3 pointNotClipped = viewPoints[nonClipIndex];
                float3 pointA = viewPoints[nextIndex];
                float3 pointB = viewPoints[prevIndex];

                // Find fractional value of the slice
                float fracA = (NEAR_CLIP - pointNotClipped.z()) / (pointA.z() - pointNotClipped.z());
                float fracB = (NEAR_CLIP - pointNotClipped.z()) / (pointB.z() - pointNotClipped.z());

                // New triangle points
                float3 clipPointAlongEdgeA = pointNotClipped.lerp(pointA, fracA);
                float3 clipPointAlongEdgeB = pointNotClipped.lerp(pointB, fracB);

                // New world points
                float3 worldClipPointAlongEdgeA = worldPoints[nonClipIndex].lerp(worldPoints[nextIndex], fracA);
                float3 worldClipPointAlongEdgeB = worldPoints[nonClipIndex].lerp(worldPoints[prevIndex], fracB);
                // New texture coord points
                float3 textureClipPointAlongEdgeA = vertexTextureCoords[nonClipIndex].lerp(vertexTextureCoords[nextIndex], fracA);
                float3 textureClipPointAlongEdgeB = vertexTextureCoords[nonClipIndex].lerp(vertexTextureCoords[prevIndex], fracB);
                // New normal points
                float3 normalClipPointAlongEdgeA = vertexNormals[nonClipIndex].lerp(vertexNormals[nextIndex], fracA);
                float3 normalClipPointAlongEdgeB = vertexNormals[nonClipIndex].lerp(vertexNormals[prevIndex], fracB);

                // Create the new triangle
                rasterizerPoints->push_back(RasterizerPoint(clipPointAlongEdgeB.z(), ViewToScreen(clipPointAlongEdgeB, fov), worldClipPointAlongEdgeB,  textureClipPointAlongEdgeB,        normalClipPointAlongEdgeB,   obj, face));
                rasterizerPoints->push_back(RasterizerPoint(pointNotClipped.z(),     ViewToScreen(pointNotClipped, fov),     worldPoints[nonClipIndex], vertexTextureCoords[nonClipIndex], vertexNormals[nonClipIndex], obj, face));
                rasterizerPoints->push_back(RasterizerPoint(clipPointAlongEdgeA.z(), ViewToScreen(clipPointAlongEdgeA, fov), worldClipPointAlongEdgeA,  textureClipPointAlongEdgeA,        normalClipPointAlongEdgeA,   obj, face));
                break;
            }
        }
    }
}