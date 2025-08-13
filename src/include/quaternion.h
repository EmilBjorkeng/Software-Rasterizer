#ifndef __QUATERNION_H__
#define __QUATERNION_H__

#define _USE_MATH_DEFINES
#include <math.h>
#include "vector.h"
#include "matrix.h"

template<typename T>
class Quaternion {
    static_assert(
        std::is_same<T, float>::value || std::is_same<T, double>::value,
        "T must be either float or double"
    );
public:
    T w, x, y, z;

    // ------------------------------
    // Constructors
    // ------------------------------

    Quaternion() : w(1), x(0), y(0), z(0) {}
    Quaternion(T w, T x, T y, T z)
        : w(w), x(x), y(y), z(z) {}

    // Yaw (Y), Pitch (X), Roll (Z)
    // Yaw Roll Pitch
    // Pitch Yaw Roll
    static Quaternion EulerToQuaternion(vec3<T> rotation) {
        T Yaw = rotation.y() * static_cast<T>(0.5);
        T Pitch = rotation.x() * static_cast<T>(0.5);
        T Roll = rotation.z() * static_cast<T>(0.5);

        T cy = cos(Yaw), sy = sin(Yaw);
        T cp = cos(Pitch), sp = sin(Pitch);
        T cr = cos(Roll), sr = sin(Roll);

        return Quaternion(
            cy * cp * cr + sy * sp * sr,  // w
            cy * sp * cr + sy * cp * sr,  // x
            sy * cp * cr - cy * sp * sr,  // y
            cy * cp * sr - sy * sp * cr   // z
        );
    }

    // Axis-angle constructor
    static Quaternion FromAxisAngle(const vec3<T>& axis, T angleRad) {
        T halfAngle = angleRad * static_cast<T>(0.5);
        T s = std::sin(halfAngle);
        vec3<T> normAxis = axis.normalized();
        return Quaternion(std::cos(halfAngle), normAxis.x() * s, normAxis.y() * s, normAxis.z() * s);
    }

    static Quaternion LookRotation(const vec3<T>& forward, const vec3<T>& up) {
        vec3<T> f = forward.normalized();
        vec3<T> r = cross(up, f).normalized();
        vec3<T> u = cross(f, r); // already orthogonalized

        // Create rotation matrix from basis vectors
        // Column-major format:
        // [ r.x  u.x  f.x ]
        // [ r.y  u.y  f.y ]
        // [ r.z  u.z  f.z ]
        T m00 = r.x(), m01 = u.x(), m02 = f.x();
        T m10 = r.y(), m11 = u.y(), m12 = f.y();
        T m20 = r.z(), m21 = u.z(), m22 = f.z();

        T trace = m00 + m11 + m22;
        if (trace > T(0)) {
            T s = std::sqrt(trace + T(1)) * T(0.5);
            T invS = T(0.25) / s;
            return Quaternion(
                s,
                (m21 - m12) * invS,
                (m02 - m20) * invS,
                (m10 - m01) * invS
            );
        } else if (m00 > m11 && m00 > m22) {
            T s = std::sqrt(T(1) + m00 - m11 - m22) * T(0.5);
            T invS = T(0.25) / s;
            return Quaternion(
                (m21 - m12) * invS,
                s,
                (m01 + m10) * invS,
                (m02 + m20) * invS
            );
        } else if (m11 > m22) {
            T s = std::sqrt(T(1) + m11 - m00 - m22) * T(0.5);
            T invS = T(0.25) / s;
            return Quaternion(
                (m02 - m20) * invS,
                (m01 + m10) * invS,
                s,
                (m12 + m21) * invS
            );
        } else {
            T s = std::sqrt(T(1) + m22 - m00 - m11) * T(0.5);
            T invS = T(0.25) / s;
            return Quaternion(
                (m10 - m01) * invS,
                (m02 + m20) * invS,
                (m12 + m21) * invS,
                s
            );
        }
    }

    // ------------------------------
    // Convertion
    // ------------------------------

    // Yaw (Y), Pitch (X), Roll (Z)
    constexpr vec3<T> ToEulerAngles() {
        vec3<T> rotation;

        // sin(pitch)
        T sp = -2 * (y * z - w * x);

        if (std::abs(sp) >= 1) {
            // Gimbal lock: pitch is +-90 degrees
            rotation.x() = std::copysign(static_cast<T>(M_PI) / 2, sp);  // Pitch (X)
            rotation.y() = std::atan2(-x * z + w * y, 0.5 - y * y - z * z); // Yaw (Y)
            rotation.z() = 0; // Roll is zero in gimbal lock
        } else {
            rotation.x() = std::asin(sp); // Pitch (X)
            rotation.y() = std::atan2(x * z + w * y, 0.5 - x * x - y * y); // Yaw (Y)
            rotation.z() = std::atan2(x * y + w * z, 0.5 - x * x - z * z); // Roll (Z)
        }

        return rotation;
    }

    // ------------------------------
    // Operator Quaternion
    // ------------------------------

    // -Quaternion
    Quaternion operator-() const {
        return Quaternion(-w, -x, -y, -z);
    }

    // Quaternion + Quaternion
    Quaternion operator+(const Quaternion& q) const {
        return Quaternion(w + q.w, x + q.x, y + q.y, z + q.z);
    }

    // Quaternion - Quaternion
    Quaternion operator-(const Quaternion& q) const {
        return Quaternion(w - q.w, x - q.x, y - q.y, z - q.z);
    }

    // Quaternion += Quaternion
    Quaternion& operator+=(const Quaternion& q) {
        w += q.w; x += q.x; y += q.y; z += q.z;
        return *this;
    }

    // Quaternion -= Quaternion
    Quaternion& operator-=(const Quaternion& q) {
        w -= q.w; x -= q.x; y -= q.y; z -= q.z;
        return *this;
    }

    // Quaternion == Quaternion
    bool operator==(const Quaternion& q) const {
        return w == q.w && x == q.x && y == q.y && z == q.z;
    }

    // Quaternion != Quaternion
    bool operator!=(const Quaternion& q) const {
        return !(*this == q);
    }

    // Quaternion multiplication (composition of rotations)
    constexpr Quaternion operator*(const Quaternion& q) const {
        return Quaternion(
            w*q.w - x*q.x - y*q.y - z*q.z,
            w*q.x + x*q.w + y*q.z - z*q.y,
            w*q.y - x*q.z + y*q.w + z*q.x,
            w*q.z + x*q.y - y*q.x + z*q.w
        );
    }

    Quaternion& operator*=(const Quaternion& q) {
        *this = *this * q;
        return *this;
    }

    // ------------------------------
    // Operator Scalar
    // ------------------------------

    // Quaternion * Scalar
    template<typename Scalar, typename = std::enable_if_t<std::is_arithmetic_v<Scalar>>>
        Quaternion operator*(Scalar scalar) const {
        return Quaternion(w * scalar, x * scalar, y * scalar, z * scalar);
    }

    // Quaternion / Scalar
    template<typename Scalar, typename = std::enable_if_t<std::is_arithmetic_v<Scalar>>>
    Quaternion operator/(Scalar scalar) const {
        return Quaternion(w / scalar, x / scalar, y / scalar, z / scalar);
    }

    // Quaternion *= Scalar
    template<typename Scalar, typename = std::enable_if_t<std::is_arithmetic_v<Scalar>>>
    Quaternion& operator*=(Scalar scalar) {
        w *= scalar; x *= scalar; y *= scalar; z *= scalar;
        return *this;
    }

    // Quaternion /= Scalar
    template<typename Scalar, typename = std::enable_if_t<std::is_arithmetic_v<Scalar>>>
    Quaternion& operator/=(Scalar scalar) {
        w /= scalar; x /= scalar; y /= scalar; z /= scalar;
        return *this;
    }

    // ------------------------------
    // Operator Vector
    // ------------------------------

    // Quaternion * Vector (Rotate vector)
    constexpr vec3<T> operator*(const vec3<T>& v) const {
        Quaternion<T> p(0, v.x(), v.y(), v.z());
        Quaternion<T> result = (*this) * p * this->Conjugate();
        return vec3<T>(result.x, result.y, result.z);
    }

    // ------------------------------
    // Operator Functions
    // ------------------------------

    // Magnitude
    [[nodiscard]]
    T Magnitude() const {
        return std::sqrt(w * w + x * x + y * y + z * z);
    }

    [[nodiscard]]
    constexpr Quaternion normalized() const noexcept {
        float len = sqrt(w*w + x*x + y*y + z*z);
        return Quaternion(w/len, x/len, y/len, z/len);
    }

    // Conjugate (used to reverse the rotation)
    [[nodiscard]]
    constexpr Quaternion Conjugate() const {
        return Quaternion(w, -x, -y, -z);
    }

    // Inverse of the quaternion
    [[nodiscard]]
    constexpr Quaternion Inverse() const {
        T norm_sq = w*w + x*x + y*y + z*z;
        return Conjugate() / norm_sq;
    }

    [[nodiscard]]
    static constexpr Quaternion identity() {
        return Quaternion(1, 0, 0, 0);
    }

    // Dot product
    T dot(const Quaternion& q) const {
        return w * q.w + x * q.x + y * q.y + z * q.z;
    }

    // Spherical Linear Interpolation
    static Quaternion slerp(const Quaternion& a, const Quaternion& b, T t) {
        T dotProd = a.dot(b);

        // If dot < 0, negate to take shortest path
        Quaternion end = (dotProd < 0) ? -b : b;
        if (dotProd < 0) dotProd = -dotProd;

        const T EPSILON = static_cast<T>(1e-6);
        if (dotProd > 1 - EPSILON) {
            // Linear interpolation for small angles
            return (a * (1 - t) + end * t).normalized();
        }

        T theta = std::acos(dotProd);
        T sinTheta = std::sin(theta);

        T scaleA = std::sin((1 - t) * theta) / sinTheta;
        T scaleB = std::sin(t * theta) / sinTheta;

        return (a * scaleA + end * scaleB).normalized();
    }
};

// Scalar * Quaternion
template<typename T, typename Scalar, typename = std::enable_if_t<std::is_arithmetic_v<Scalar>>>
inline Quaternion<T> operator*(Scalar scalar, const Quaternion<T>& q) {
    return q * scalar;
}

// Stream output operator
template<typename T>
std::ostream& operator<<(std::ostream& os, const Quaternion<T>& q) {
    os << '(' << q.w << ' ' << q.x << "i " << q.y << "j " << q.z << "k)";
    return os;
}

#endif