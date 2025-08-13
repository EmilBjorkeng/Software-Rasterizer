#ifndef __RENDER_H__
#define __RENDER_H__

#include "vector.h"
#include "camera.h"
#include "obj.h"
#include "light.h"
#include <vector>
#include <cstdint>

#define SCREEN_WIDTH 720
#define SCREEN_HEIGHT 480

void RendersceneObjects(uint32_t *drawBuffer, int pitch, float *depthBuffer, const Camera<float> &camera,
    const std::vector<const Object*> &sceneObjects, const std::vector<const Light*> &light);

float2 ViewToScreen(const float3 vertex_view, int fov);
std::vector<int> TriangulateEarClipping(const std::vector<float3>& vertices);

struct Triangle {
    float3 a, b, c;
    float3 normal;
    float3 aabbMin, aabbMax;
};

struct RasterizerPoint {
    float depth;
    float2 screenPos;
    float3 worldSpaceVertex;
    float2 vertexTextureCoords;
    float3 vertexNormals;

    // TODO: DO BETTER  |  FIX
    const Face *vertex_face;
    const Object *object;

    RasterizerPoint(float viewDepth, float3 screenVertex, float3 worldVertex, float3 textureCoords, float3 normal, const Object *obj, const Face &face)
        : depth(viewDepth), screenPos(screenVertex), worldSpaceVertex(worldVertex),
        vertexTextureCoords(textureCoords), vertexNormals(normal),
        vertex_face(&face), object(obj)
    {}
};
void triangleSlicer(std::vector<RasterizerPoint> *rasterizerPoints,
    const std::vector<float3> &worldSpaceVertices, const std::vector<float3> &viewSpaceVertices,
    const std::vector<float3> &worldSpaceNormals, const std::vector<int> &vertexOrder, int fov,
    const Face &face, const Object *obj);

#endif