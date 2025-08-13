#ifndef __CAMERA_H__
#define __CAMERA_H__

#include "vector.h"
#include "matrix.h"
#include "quaternion.h"
#include <type_traits>

template<typename T>
class Camera {
    static_assert(
        std::is_same<T, float>::value || std::is_same<T, double>::value,
        "T must be either float or double"
    );
public:
    vec3<T> position;
    Quaternion<T> rotation;
    int fov = 90;

    // ------------------------------
    // Constructors
    // ------------------------------

    Camera()
    : position(0, 0, 0), rotation(Quaternion<T>::Identity()) {}

    // Constructor with position and quaternion
    template<typename Derived>
    Camera(const vec3_base<Derived, T>& initPosition, const Quaternion<T>& initRotation)
            : position(initPosition), rotation(initRotation){}

    // Constructor with position and Euler angles (pitch, yaw, roll)
    template<typename Derived>
    Camera(const vec3_base<Derived, T>& initPosition, const vec3_base<Derived, T> initRotation)
        : position(initPosition), rotation(Quaternion<T>::EulerToQuaternion(initRotation)) {}


    // ------------------------------
    // Direction vectors
    // ------------------------------

    vec3<T> Forward() const {
        return rotation * vec3<T>(0, 0, 1);
    }

    vec3<T> Right() const {
        return rotation * vec3<T>(1, 0, 0);
    }

    vec3<T> Up() const {
        return rotation * vec3<T>(0, 1, 0);
    }

    // ------------------------------
    // View Matrix
    // ------------------------------

    constexpr mat4<T> GetViewMatrix() const {
        Quaternion<T> invRot = rotation.Inverse();
        mat4<T> FlipY = mat4<T>::Scale(float3(1, -1, 1)); // Flip Y -> +Y points UP
        return mat4<T>::Rotation(invRot) * FlipY * mat4<T>::Translate(-position);
    }
};

#endif