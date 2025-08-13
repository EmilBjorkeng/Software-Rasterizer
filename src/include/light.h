#ifndef __LIGHT_H__
#define __LIGHT_H__

#include "vector.h"

class Light {
protected:
    explicit Light(float intensity, float3 color = float3(1.0f, 1.0f, 1.0f))
        : intensity(intensity), color(color) {}
public:
    float intensity;
    float3 color;

    virtual float3 lightDirection(const float3& point) const = 0;
    virtual float distance(const float3& point) const = 0;
    virtual float3 getPosition(const float3 &point) const = 0;
    virtual ~Light() = default;
};

class PointLight : public Light {
public:
    float3 position;

    PointLight(float3 position, float intensity)
        : Light(intensity), position(position) {}

    PointLight(float3 position, float intensity, float3 color)
        : Light(intensity, color), position(position) {}

    float3 lightDirection(const float3 &point) const override {
        return (position - point).normalized();
    }

    float distance(const float3 &point) const override {
        return (position - point).length();
    }

    float3 getPosition(const float3&) const override {
        return position;
    }
};

class DirectionalLight : public Light {
public:
    float3 direction;

    DirectionalLight(float3 direction, float intensity)
        : Light(intensity), direction(direction.normalized()) {}

    DirectionalLight(float3 direction, float intensity, float3 color)
        : Light(intensity, color), direction(direction.normalized()) {}

    float3 lightDirection(const float3&) const override {
        return -direction; // constant for all points
    }

    float distance(const float3&) const override {
        return 1; // Directional light has no distance cutoff
    }

    float3 getPosition(const float3 &point) const override {
        return point - (direction * INFINITY); // Set the lights position as far away from the point as possible
    }
};

class Spotlight : public Light {
private:
    bool isInsideCone(const float3 &point) const {
        float3 L = (point - position).normalized();
        return L.dot(direction) > cosCutoff;
    }
public:
    float3 position;
    float3 direction;
    float cutoffAngle;
    float cosCutoff;

    Spotlight(float3 position, float3 direction, float cutoff, float intensity)
        : Light(intensity), position(position), direction(direction.normalized()),
          cutoffAngle(cutoff), cosCutoff(cos(cutoff)) {}

    Spotlight(float3 position, float3 direction, float cutoff, float intensity, float3 color)
        : Light(intensity, color), position(position), direction(direction.normalized()),
          cutoffAngle(cutoff), cosCutoff(cos(cutoff)) {}

    float3 lightDirection(const float3& point) const override {
        if (!isInsideCone(point))
            return float3(0, 0, 0);

        return (position - point).normalized();
    }
    
    float distance(const float3& point) const override {
        return (position - point).length();
    }

    float3 getPosition(const float3&) const override {
        return position;
    }
};

#endif