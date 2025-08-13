#ifndef __MATRIX_H__
#define __MATRIX_H__

#include <array>
#include <initializer_list>
#include <stdexcept>
#include <type_traits>
#include <cstddef>
#include <iostream>

// Vector template for the constructor
template<typename Derived, typename T, std::size_t Dimensions>
class vec_core;

// Quaternion template
template<typename T>
class Quaternion;

template<typename T, std::size_t Rows, std::size_t Cols>
class matrix {
    static_assert(std::is_arithmetic<T>::value, "T must be arithmetic.");
    static_assert(Rows > 0, "Rows must be over 0.");
    static_assert(Cols > 0, "Cols must be over 0.");

    // Convert an array of type U into an array of type T
    template<typename U>
    static constexpr std::array<T, Rows * Cols> convert_array(const std::array<U, Rows * Cols>& other) {
        std::array<T, Rows * Cols> result{};
        for (std::size_t i = 0; i < Rows * Cols; ++i)
            result[i] = static_cast<T>(other[i]);
        return result;
    }
public:
    std::array<T, Rows * Cols> data;

    // ------------------------------
    // Constructors
    // ------------------------------

    constexpr matrix() noexcept : data{} {}
    constexpr matrix(std::initializer_list<T> list) {
        if (list.size() > Rows * Cols)
            throw std::length_error("Too many initializers for matrix");
        std::size_t i = 0;
        for (auto v : list) {
            if (i < Rows * Cols)
                data[i++] = v;
        }
        for (; i < Rows * Cols; ++i)
            data[i] = T(0);
    }
    constexpr matrix(std::array<T, Rows * Cols> arr) : data(arr) {}

    // matrix<T, Rows, Cols> = matrix<OtherT, OtherRows, OtherCols>
    template<typename OtherT, std::size_t OtherRows, std::size_t OtherCols>
    constexpr matrix(const matrix<OtherT, OtherRows, OtherCols>& other) noexcept {
        // Initialize all elements to zero
        for (std::size_t i = 0; i < Rows * Cols; ++i)
            data[i] = T(0);

        // Copy overlapping elements from the source matrix
        std::size_t min_Rows = (Rows < OtherRows) ? Rows : OtherRows;
        std::size_t min_Cols = (Cols < OtherCols) ? Cols : OtherCols;
        for (std::size_t r = 0; r < min_Rows; ++r)
            for (std::size_t c = 0; c < min_Cols; ++c)
                (*this)(r, c) = static_cast<T>(other(r, c));
    }

    // matrix = vec (matrix_vector_bridge.h)
    template<typename Derived, std::size_t Dimensions>
    constexpr matrix(const vec_core<Derived, T, Dimensions>& v);

    // ------------------------------
    // Indexing
    // ------------------------------

    auto begin() { return data.begin(); }
    auto end() { return data.end(); }
    auto begin() const { return data.begin(); }
    auto end() const { return data.end(); }

    static constexpr std::size_t row_count() noexcept { return Rows; }
    static constexpr std::size_t col_count() noexcept { return Cols; }

    constexpr std::size_t size() const noexcept {
        return Rows * Cols;
    }

    [[nodiscard]]
    constexpr T& operator()(std::size_t row, std::size_t col) {
        return data[row * Cols + col];
    }

    [[nodiscard]]
    constexpr const T& operator()(std::size_t row, std::size_t col) const {
        return data[row * Cols + col];
    }

    [[nodiscard]]
    constexpr T& operator[](std::size_t i) { return data[i]; }

    [[nodiscard]]
    constexpr const T& operator[](std::size_t i) const { return data[i]; }

    [[nodiscard]]
    constexpr T& at(std::size_t row, std::size_t col) {
        return data[row * Cols + col];
    }

    [[nodiscard]]
    constexpr const T& at(std::size_t row, std::size_t col) const {
        return data[row * Cols + col];
    }

    // ------------------------------
    // Row and Col
    // ------------------------------

    // Return a col
    template<std::size_t N = Cols>
    constexpr std::array<T, N> col(std::size_t r) const {
        static_assert(N <= Cols, "Requested row length exceeds number of columns");
        std::array<T, N> result{};
        for (std::size_t i = 0; i < N; ++i)
            result[i] = (*this)(r, i);
        return result;
    }

    // Return a row
    template<std::size_t N = Rows>
    constexpr std::array<T, N> row(std::size_t c) const {
        static_assert(N <= Rows, "Requested column length exceeds number of rows");
        std::array<T, N> result{};
        for (std::size_t i = 0; i < N; ++i)
            result[i] = (*this)(i, c);
        return result;
    }

    // ------------------------------
    // Operator Matrix
    // ------------------------------

    // -Matrix
    constexpr matrix operator-() const noexcept {
        matrix result;
        for (std::size_t i = 0; i < Rows * Cols; ++i)
            result.data[i] = -data[i];
        return result;
    }

    // Matrix + Matrix
    constexpr matrix operator+(const matrix& other) const noexcept {
        matrix result;
        for (std::size_t i = 0; i < Rows * Cols; ++i)
            result.data[i] = data[i] + other.data[i];
        return result;
    }

    // Matrix - Matrix
    constexpr matrix operator-(const matrix& other) const noexcept {
        matrix result;
        for (std::size_t i = 0; i < Rows * Cols; ++i)
            result.data[i] = data[i] - other.data[i];
        return result;
    }

    // Matrix += Matrix
    constexpr matrix& operator+=(const matrix& other) noexcept {
        for (std::size_t i = 0; i < Rows * Cols; ++i)
            data[i] += other.data[i];
        return *this;
    }

    // Matrix -= Matrix
    constexpr matrix& operator-=(const matrix& other) noexcept {
        for (std::size_t i = 0; i < Rows * Cols; ++i)
            data[i] -= other.data[i];
        return *this;
    }

    // Matrix == Matrix
    template<typename OtherT>
    constexpr bool operator==(const matrix<OtherT, Rows, Cols>& other) const noexcept {
        for (std::size_t i = 0; i < Rows * Cols; ++i)
            if (data[i] != static_cast<T>(other.data[i]))
                return false;
        return true;
    }

    // Matrix != Matrix
    constexpr bool operator!=(const matrix& other) const noexcept {
        return !(*this == other);
    }

    // Matrix *= Matrix
    template<std::size_t N = Cols>
    matrix& operator*=(const matrix<T, N, N>& other) {
        static_assert(Rows == Cols && Cols == N, "Matrix *= only defined for square matrices.");
        *this = *this * other;
        return *this;
    }

    // ------------------------------
    // Operator Scalar
    // ------------------------------

    // Matrix * Scalar
    template<typename Scalar, typename = std::enable_if_t<std::is_arithmetic_v<Scalar>>>
    constexpr matrix operator*(Scalar scalar) const noexcept {
        matrix result;
        for (std::size_t i = 0; i < Rows * Cols; ++i)
            result.data[i] = data[i] * scalar;
        return result;
    }

    // Matrix / Scalar
    template<typename Scalar, typename = std::enable_if_t<std::is_arithmetic_v<Scalar>>>
    constexpr matrix operator/(Scalar scalar) const noexcept {
        matrix result;
        for (std::size_t i = 0; i < Rows * Cols; ++i)
            result.data[i] = data[i] / scalar;
        return result;
    }

    // Matrix *= Scalar
    template<typename Scalar, typename = std::enable_if_t<std::is_arithmetic_v<Scalar>>>
    constexpr matrix& operator*=(Scalar scalar) noexcept {
        for (auto& v : data) v *= scalar;
        return *this;
    }

    // Matrix /= Scalar
    template<typename Scalar, typename = std::enable_if_t<std::is_arithmetic_v<Scalar>>>
    constexpr matrix& operator/=(Scalar scalar) noexcept {
        for (auto& v : data) v /= scalar;
        return *this;
    }

    // ------------------------------
    // Operator Functions
    // ------------------------------

    constexpr void fill(T value) {
        for (auto& v : data) v = value;
    }

    // Transpose
    [[nodiscard]]
    constexpr matrix<T, Cols, Rows> transpose() const noexcept {
        matrix<T, Cols, Rows> result;
        for (std::size_t r = 0; r < Rows; ++r)
            for (std::size_t c = 0; c < Cols; ++c)
                result(c, r) = (*this)(r, c);
        return result;
    }

    // Inverse (square only)
    [[nodiscard]]
    constexpr matrix<T, Cols, Rows> inverse() const {
        static_assert(Rows == Cols, "Inverse only works on a square matrix.");
        static constexpr std::size_t N = Rows;

        matrix<T, N, N> result = matrix::Identity();
        matrix<T, N, N> this_matrix = *this;

        for (std::size_t c = 0; c < N; ++c) {
            // Find pivot in column c
            std::size_t pivot = c;
            for (std::size_t row = c + 1; row < N; ++row) {
                if ((T)std::abs(this_matrix(row, c)) > (T)std::abs(this_matrix(pivot, c))) {
                    pivot = row;
                }
            }

            // Check for singularity
            if ((T)std::abs(this_matrix(pivot, c)) < std::numeric_limits<T>::epsilon()) {
                // Singular matrix: no inverse
                return {};
            }

            // Swap rows to move pivot to diagonal
            if (pivot != c) {
                for (std::size_t j = 0; j < N; ++j) {
                    std::swap(this_matrix(c, j), this_matrix(pivot, j));
                    std::swap(result(c, j), result(pivot, j));
                }
            }

            // Normalize pivot row
            T divisor = this_matrix(c, c);
            for (std::size_t j = 0; j < N; ++j) {
                this_matrix(c, j) /= divisor;
                result(c, j) /= divisor;
            }

            // Eliminate other rows
            for (std::size_t r = 0; r < N; ++r) {
                if (r != c) {
                    T factor = this_matrix(r, c);
                    for (std::size_t j = 0; j < N; ++j) {
                        this_matrix(r, j) -= factor * this_matrix(c, j);
                        result(r, j) -= factor * result(c, j);
                    }
                }
            }
        }

        return result;
    }

    // ------------------------------
    // Matrix types
    // ------------------------------

    // Create a Matrix with all elements set to 0
    static constexpr matrix Zero() noexcept {
        return matrix{};
    }

    // Create a Matrix with all elements set to 1
    static constexpr matrix One() noexcept {
        return Filled(T(1));
    }

    // Create a Matrix with all elements set to the given value
    static constexpr matrix Filled(T value) noexcept {
        matrix result;
        result.fill(value);
        return result;
    }

    // Identity Matrix (square only)
    static constexpr matrix Identity() noexcept {
        static_assert(Rows == Cols, "Identity is only for square matrices.");
        matrix result;
        for (std::size_t i = 0; i < Rows; ++i)
            result(i, i) = T(1);
        return result;
    }

    // Translate Matrix (square only)
    static constexpr matrix Translate(float3 v) noexcept {
        static_assert(Rows == Cols, "Translate matrix has to be a square matrix.");
        static_assert(Rows >= 4, "Translate matrix has to be a 4x4 matrix or bigger.");

        matrix result = matrix::Identity();
        result.at(0, 3) = v.x();
        result.at(1, 3) = v.y();
        result.at(2, 3) = v.z();
        return result;
    }

    // Scale Matrix (square only)
    static constexpr matrix Scale(float3 v) noexcept {
        static_assert(Rows == Cols, "Scale matrix has to be a square matrix.");
        static_assert(Rows >= 3, "Scale matrix has to be a 3x3 matrix or bigger.");

        matrix result = matrix::Identity();
        result.at(0, 0) = v.x();
        result.at(1, 1) = v.y();
        result.at(2, 2) = v.z();
        return result;
    }

    // RotateX Matrix (square only)
    static matrix RotateX(float angle) noexcept {
        static_assert(Rows == Cols, "RotateX matrix has to be a square matrix.");
        static_assert(Rows >= 3, "RotateX matrix has to be a 3x3 matrix or bigger.");

        matrix result = matrix::Identity();
        result.at(1, 1) = (T)cos(angle);
        result.at(1, 2) = (T)-sin(angle);
        result.at(2, 1) = (T)sin(angle);
        result.at(2, 2) = (T)cos(angle);
        return result;
    }

    // RotateY Matrix (square only)
    static matrix RotateY(float angle) noexcept {
        static_assert(Rows == Cols, "RotateY matrix has to be a square matrix.");
        static_assert(Rows >= 3, "RotateX matrix has to be a 3x3 matrix or bigger.");

        matrix result = matrix::Identity();
        result.at(0, 0) = (T)cos(angle);
        result.at(0, 2) = (T)sin(angle);
        result.at(2, 0) = (T)-sin(angle);
        result.at(2, 2) = (T)cos(angle);
        return result;
    }

    // RotateZ Matrix (square only)
    static matrix RotateZ(float angle) noexcept {
        static_assert(Rows == Cols, "RotateZ matrix has to be a square matrix.");
        static_assert(Rows >= 3, "RotateX matrix has to be a 3x3 matrix or bigger.");

        matrix result = matrix::Identity();
        result.at(0, 0) = (T)cos(angle);
        result.at(0, 1) = (T)-sin(angle);
        result.at(1, 0) = (T)sin(angle);
        result.at(1, 1) = (T)cos(angle);
        return result;
    }

    // Rotation Matrix (square only)
    static matrix Rotation(Quaternion<T> q) noexcept {
        static_assert(Rows == Cols, "Rotation matrix has to be a square matrix.");
        static_assert(Rows >= 3, "Rotation matrix has to be a 3x3 matrix or bigger.");

        q = q.normalized();

        float xx = q.x * q.x, yy = q.y * q.y, zz = q.z * q.z;
        float xy = q.x * q.y, xz = q.x * q.z, yz = q.y * q.z;
        float wx = q.w * q.x, wy = q.w * q.y, wz = q.w * q.z;

        matrix result = matrix::Identity();
        result.at(0, 0) = 1 - 2 * (yy + zz);
        result.at(0, 1) = 2 * (xy - wz);
        result.at(0, 2) = 2 * (xz + wy);
        result.at(1, 0) = 2 * (xy + wz);
        result.at(1, 1) = 1 - 2 * (xx + zz);
        result.at(1, 2) = 2 * (yz - wx);
        result.at(2, 0) = 2 * (xz - wy);
        result.at(2, 1) = 2 * (yz + wx);
        result.at(2, 2) = 1 - 2 * (xx + yy);
        return result;
    }
};

// Matrix Multiplication
template<typename T, std::size_t R, std::size_t K, std::size_t C>
constexpr matrix<T, R, C> operator*(const matrix<T, R, K>& a, const matrix<T, K, C>& b) {
    matrix<T, R, C> result;
    for (std::size_t i = 0; i < R; ++i) {
        for (std::size_t j = 0; j < C; ++j) {
            auto sum = T{};
            for (std::size_t k = 0; k < K; ++k) {
                sum += a(i, k) * b(k, j);
            }
            result(i, j) = sum;
        }
    }
    return result;
}

// Scalar * Matrix
template<typename T, std::size_t R, std::size_t C>
constexpr matrix<T, R, C> operator*(T scalar, const matrix<T, R, C>& mat) {
    return mat * scalar;
}

// Stream output operator
template<typename T, std::size_t R, std::size_t C>
std::ostream& operator<<(std::ostream& os, const matrix<T, R, C>& mat) {
    size_t maxLength = 0;
    for (const auto& num : mat.data) {
        size_t len = std::to_string(num).size();
        if (len > maxLength) maxLength = len;
    }

    for (std::size_t r = 0; r < R; ++r) {
        os << '[';
        for (std::size_t c = 0; c < C; ++c) {
            std::string str = std::to_string(mat(r, c));
            os << str << std::string(maxLength - str.size(), ' ');
            /*if (c < C-1)
                os << str << ',' << std::string(maxLength - str.size(), ' ');
            else
                os << str << std::string(maxLength - str.size(), ' ');*/
            if (c < C-1) os << ' ';
        }
        os << ']';
        if (r < R-1) os << '\n';
    }
    return os;
}

// Aliases
template<typename T> using mat2 = matrix<T, 2, 2>;
template<typename T> using mat3 = matrix<T, 3, 3>;
template<typename T> using mat4 = matrix<T, 4, 4>;
template<typename T> using mat5 = matrix<T, 5, 5>;
template<typename T> using mat6 = matrix<T, 6, 6>;
template<typename T> using mat7 = matrix<T, 7, 7>;
template<typename T> using mat8 = matrix<T, 8, 8>;
template<typename T> using mat9 = matrix<T, 9, 9>;

#endif