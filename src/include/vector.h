#ifndef __VECTOR_H__
#define __VECTOR_H__

#define _USE_MATH_DEFINES
#include <math.h>
#include <cstddef>
#include <type_traits>
#include <array>
#include <iostream>
#include <limits>
#include <stdexcept>

// Matrix template for the constructor
template<typename T, std::size_t Rows, std::size_t Cols>
class matrix;

template<typename Derived, typename T, std::size_t Dimensions>
class vec_core {
    static_assert(std::is_arithmetic<T>::value, "T must be an arithmetic type.");
    static_assert(Dimensions > 1, "Dimensions must be over 1.");

    // Convert an array of type U into an array of type T
    template<typename U>
    static constexpr std::array<T, Dimensions> convert_array(const std::array<U, Dimensions>& other) {
        std::array<T, Dimensions> result{};
        for (std::size_t i = 0; i < Dimensions; ++i)
            result[i] = static_cast<T>(other[i]);
        return result;
    }
public:
    std::array<T, Dimensions> data;

    // ------------------------------
    // Constructors
    // ------------------------------

    constexpr vec_core() noexcept : data{} {}
    vec_core(std::initializer_list<T> list) {
        if (list.size() > Dimensions)
            throw std::length_error("Too many initializers for the vector");

        std::size_t i = 0;
        for (auto v : list) {
            if (i < Dimensions)
                data[i++] = v;
        }
        for (; i < Dimensions; ++i)
            data[i] = T(0);
    }
    constexpr vec_core(const std::array<T, Dimensions>& arr) : data(arr) {}

    // vec<T, D> = vec<OtherT, OtherD>
    template<typename OtherDerived, typename OtherT, std::size_t OtherDimensions>
    constexpr vec_core(const vec_core<OtherDerived, OtherT, OtherDimensions>& other) noexcept {
        std::array<T, OtherDimensions> otherData;
        for (std::size_t i = 0; i < OtherDimensions; ++i)
            otherData[i] = static_cast<T>(other.data[i]);

        std::size_t min_dim = (Dimensions < OtherDimensions) ? Dimensions : OtherDimensions;
        for (std::size_t i = 0; i < min_dim; ++i)
            data[i] = otherData[i];
        for (std::size_t i = min_dim; i < Dimensions; ++i)
            data[i] = T(0);
    }

    // vec<T, D> = vec<T, D-1, Value>
    template<typename OtherDerived>
    constexpr vec_core(const vec_core<OtherDerived, T, (Dimensions - 1)>& other, T value) noexcept {
        std::size_t i = 0;
        for (; i < Dimensions - 1; ++i)
            data[i] = other.data[i];
        data[i] = value;
    }

    // vec = matrix (matrix_vector_bridge.h)
    template<std::size_t Rows, std::size_t Cols>
    constexpr vec_core(const matrix<T, Rows, Cols>& mat);

    // ------------------------------
    // Indexing
    // ------------------------------

    auto begin() { return data.begin(); }
    auto end() { return data.end(); }
    auto begin() const { return data.begin(); }
    auto end() const { return data.end(); }

    [[nodiscard]]
    constexpr T& operator[](std::size_t i) { return data[i]; }

    [[nodiscard]]
    constexpr const T& operator[](std::size_t i) const { return data[i]; }

    // ------------------------------
    // Operator Vector
    // ------------------------------

    // -Vector
    constexpr Derived operator-() const noexcept {
        std::array<T, Dimensions> result{};
        for (std::size_t i = 0; i < Dimensions; ++i) {
            result[i] = -data[i];
        }
        return Derived(result);
    }

    // Vector + Vector
    constexpr Derived operator+(const Derived& other) const noexcept {
        std::array<T, Dimensions> result{};
        for (std::size_t i = 0; i < Dimensions; ++i) {
            result[i] = data[i] + other.data[i];
        }
        return Derived(result);
    }

    // Vector - Vector
    constexpr Derived operator-(const Derived& other) const noexcept {
        std::array<T, Dimensions> result{};
        for (std::size_t i = 0; i < Dimensions; ++i) {
            result[i] = data[i] - other.data[i];
        }
        return Derived(result);
    }

    // Vector += Vector
    constexpr Derived& operator+=(const Derived& other) noexcept {
        for (std::size_t i = 0; i < Dimensions; ++i) {
            data[i] += other.data[i];
        }
        return *static_cast<Derived*>(this);
    }

    // Vector -= Vector
    constexpr Derived& operator-=(const Derived& other) noexcept {
        for (std::size_t i = 0; i < Dimensions; ++i) {
            data[i] -= other.data[i];
        }
        return *static_cast<Derived*>(this);
    }

    // Vector == Vector
    template<typename OtherDerived, typename OtherT>
    constexpr bool operator==(const vec_core<OtherDerived, OtherT, Dimensions>& other) const noexcept {
        for (std::size_t i = 0; i < Dimensions; ++i)
            if (data[i] != static_cast<T>(other.data[i]))
                return false;
        return true;
    }

    // Vector != Vector
    constexpr bool operator!=(const Derived& other) const noexcept {
        return !(*this == other);
    }

    // Vector * Vector
    constexpr Derived operator*(const Derived& other) const noexcept {
        std::array<T, Dimensions> result{};
        for (std::size_t i = 0; i < Dimensions; ++i) {
            result[i] = data[i] * other.data[i];
        }
        return Derived(result);
    }

    // Vector / Vector
    constexpr Derived operator/(const Derived& other) const noexcept {
        std::array<T, Dimensions> result{};
        for (std::size_t i = 0; i < Dimensions; ++i) {
            result[i] = data[i] / other.data[i];
        }
        return Derived(result);
    }

    // Vector *= Vector
    constexpr Derived& operator*=(const Derived& other) noexcept {
        for (std::size_t i = 0; i < Dimensions; ++i) {
            data[i] *= other.data[i];
        }
        return *static_cast<Derived*>(this);
    }

    // Vector /= Vector
    constexpr Derived& operator/=(const Derived& other) noexcept {
        for (std::size_t i = 0; i < Dimensions; ++i) {
            data[i] /= other.data[i];
        }
        return *static_cast<Derived*>(this);
    }

    // ------------------------------
    // Operator Scalar
    // ------------------------------

    // Vector * Scalar
    template<typename Scalar, typename = std::enable_if_t<std::is_arithmetic_v<Scalar>>>
    constexpr Derived operator*(Scalar scalar) const noexcept {
        std::array<T, Dimensions> result{};
        for (std::size_t i = 0; i < Dimensions; ++i) {
            result[i] = static_cast<T>(data[i] * scalar);
        }
        return Derived(result);
    }

    // Vector / Scalar
    template<typename Scalar, typename = std::enable_if_t<std::is_arithmetic_v<Scalar>>>
    constexpr Derived operator/(Scalar scalar) const noexcept {
        std::array<T, Dimensions> result{};
        for (std::size_t i = 0; i < Dimensions; ++i) {
            result[i] = static_cast<T>(data[i] / scalar);
        }
        return Derived(result);
    }

    // Vector *= Scalar
    template<typename Scalar, typename = std::enable_if_t<std::is_arithmetic_v<Scalar>>>
    constexpr Derived& operator*=(Scalar scalar) noexcept {
        for (std::size_t i = 0; i < Dimensions; ++i) {
            data[i] = static_cast<T>(data[i] * scalar);
        }
        return *static_cast<Derived*>(this);
    }

    // Vector /= Scalar
    template<typename Scalar, typename = std::enable_if_t<std::is_arithmetic_v<Scalar>>>
    constexpr Derived& operator/=(Scalar scalar) noexcept {
        for (std::size_t i = 0; i < Dimensions; ++i) {
            data[i] = static_cast<T>(data[i] / scalar);
        }
        return *static_cast<Derived*>(this);
    }

    // ------------------------------
    // Operator Functions
    // ------------------------------

    constexpr void fill(T value) noexcept {
        data.fill(value);
    }

    // Create a Vector with all elements set to 0
    static constexpr Derived Zero() noexcept {
        return Derived{};
    }

    // Create a Vector with all elements set to 1
    static constexpr Derived One() noexcept {
        return Filled(T(1));
    }

    // Create a Vector with all elements set to the given value
    static constexpr Derived Filled(T value) noexcept {
        Derived result;
        result.fill(value);
        return result;
    }

    constexpr Derived min(const Derived& other) const noexcept {
        std::array<T, Dimensions> result{};
        for (std::size_t i = 0; i < Dimensions; ++i) {
            result[i] = (T)std::min(data[i], other.data[i]);
        }
        return Derived(result);
    }

    constexpr Derived max(const Derived& other) const noexcept {
        std::array<T, Dimensions> result{};
        for (std::size_t i = 0; i < Dimensions; ++i) {
            result[i] = (T)std::max(data[i], other.data[i]);
        }
        return Derived(result);
    }

    [[nodiscard]]
    constexpr T dot(const Derived& other) const noexcept {
        T result = T(0);
        for (std::size_t i = 0; i < Dimensions; ++i) {
            result += data[i] * other.data[i];
        }
        return result;
    }

    [[nodiscard]]
    constexpr T length_squared() const noexcept {
        T result = T(0);
        for (std::size_t i = 0; i < Dimensions; ++i) {
            result += data[i] * data[i];
        }
        return result;
    }

    [[nodiscard]]
    constexpr auto length() const noexcept {
        if constexpr (std::is_same_v<T, float>)
            return std::sqrt(length_squared());
        else
            return std::sqrt(static_cast<double>(length_squared()));
    }

    [[nodiscard]]
    constexpr Derived normalized() const noexcept {
        static_assert(std::is_floating_point<T>::value, "normalized() requires floating-point type.");
        auto len = length();
        if (std::abs(len) > std::numeric_limits<T>::epsilon()) {
            std::array<T, Dimensions> result{};
            for (std::size_t i = 0; i < Dimensions; ++i) {
                result[i] = static_cast<T>(data[i] / len);
            }
            return Derived(result);
        }
        return Derived::Zero();
    }

    // Linear Interpolation
    [[nodiscard]]
    constexpr Derived lerp(const Derived& other, T t) const noexcept {
        if constexpr (std::is_floating_point<T>::value) {
            if (t < T(0)) t = T(0);
            if (t > T(1)) t = T(1);
        }
        return (*this) * (T(1) - t) + other * t;
    }
};

// Scalar * Vector
template<typename Derived, typename T, std::size_t Dimensions, typename Scalar, typename = std::enable_if_t<std::is_arithmetic_v<Scalar>>>
constexpr Derived operator*(Scalar scalar, const vec_core<Derived, T, Dimensions>& v) noexcept {
    return v * scalar;
}

// Scalar / Vector
template<typename Derived, typename T, std::size_t Dimensions, typename Scalar, typename = std::enable_if_t<std::is_arithmetic_v<Scalar>>>
constexpr Derived operator/(Scalar scalar, const vec_core<Derived, T, Dimensions>& v) noexcept {
        std::array<T, Dimensions> result{};
        for (std::size_t i = 0; i < Dimensions; ++i) {
            result[i] = static_cast<T>(scalar / v[i]);
        }
        return Derived(result);
}

// Stream output operator
template<typename Derived, typename T, std::size_t Dimensions>
std::ostream& operator<<(std::ostream& os, const vec_core<Derived, T, Dimensions>& v) {
    os << '(';
    for (std::size_t i = 0; i < Dimensions; ++i) {
        os << v.data[i];
        if (i + 1 < Dimensions)
            os << ", ";
    }
    os << ')';
    return os;
}

// Vec
template<typename T, std::size_t Dimensions>
class vec : public vec_core<vec<T, Dimensions>, T, Dimensions> {
public:
    using vec_core<vec<T, Dimensions>, T, Dimensions>::vec_core;
};


// ------------------------------
// Vector2
// ------------------------------
template<typename Derived, typename T>
class vec2_base : public vec_core<Derived, T, 2> {
public:
    using vec_core<Derived, T, 2>::vec_core;

    constexpr vec2_base(T x, T y) noexcept {
        this->data[0] = x;
        this->data[1] = y;
    }

    T& x() { return this->data[0]; }
    T& y() { return this->data[1]; }
    const T& x() const { return this->data[0]; }
    const T& y() const { return this->data[1]; }
};

// Vec2
template<typename T>
class vec2 : public vec2_base<vec2<T>, T> {
public:
    using vec2_base<vec2<T>, T>::vec2_base;
};

class float2 : public vec2_base<float2, float> {
public:
    using vec2_base::vec2_base;
};

class double2 : public vec2_base<double2, double> {
public:
    using vec2_base::vec2_base;
};

class int2 : public vec2_base<int2, int> {
public:
    using vec2_base::vec2_base;
};

// ------------------------------
// Vector3
// ------------------------------
template<typename Derived, typename T>
class vec3_base : public vec_core<Derived, T, 3> {
public:
    using vec_core<Derived, T, 3>::vec_core;

        constexpr vec3_base(T x, T y, T z) noexcept {
        this->data[0] = x;
        this->data[1] = y;
        this->data[2] = z;
    }

    T& x() { return this->data[0]; }
    const T& x() const { return this->data[0]; }

    T& y() { return this->data[1]; }
    const T& y() const { return this->data[1]; }

    T& z() { return this->data[2]; }
    const T& z() const { return this->data[2]; }

    constexpr Derived cross(const Derived& other) const noexcept {
        return Derived(
            this->y() * other.z() - this->z() * other.y(),
            this->z() * other.x() - this->x() * other.z(),
            this->x() * other.y() - this->y() * other.x()
        );
    }
};

// Vec3
template<typename T>
class vec3 : public vec3_base<vec3<T>, T> {
public:
    using vec3_base<vec3<T>, T>::vec3_base;
};

class float3 : public vec3_base<float3, float> {
public:
    using vec3_base::vec3_base;
};

class double3 : public vec3_base<double3, double> {
public:
    using vec3_base::vec3_base;
};

class int3 : public vec3_base<int3, int> {
public:
    using vec3_base::vec3_base;
};

// ------------------------------
// Vector4
// ------------------------------
template<typename Derived, typename T>
class vec4_base : public vec_core<Derived, T, 4> {
public:
    using vec_core<Derived, T, 4>::vec_core;

        constexpr vec4_base(T x, T y, T z, T w) noexcept {
        this->data[0] = x;
        this->data[1] = y;
        this->data[2] = z;
        this->data[3] = w;
    }

    T& x() { return this->data[0]; }
    const T& x() const { return this->data[0]; }

    T& y() { return this->data[1]; }
    const T& y() const { return this->data[1]; }

    T& z() { return this->data[2]; }
    const T& z() const { return this->data[2]; }

    T& w() { return this->data[3]; }
    const T& w() const { return this->data[3]; }
};

// Vec4
template<typename T>
class vec4 : public vec4_base<vec4<T>, T> {
public:
    using vec4_base<vec4<T>, T>::vec4_base;
};

class float4 : public vec4_base<float4, float> {
public:
    using vec4_base::vec4_base;
};

class double4 : public vec4_base<double4, double> {
public:
    using vec4_base::vec4_base;
};

class int4 : public vec4_base<int4, int> {
public:
    using vec4_base::vec4_base;
};

#endif