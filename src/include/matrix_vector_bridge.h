#ifndef __MATRIX_VECTOR_BRIDGE_H__
#define __MATRIX_VECTOR_BRIDGE_H__

#include "vector.h"
#include "matrix.h"

// ------------------------------
// Constructors
// ------------------------------

// Vector = Matrix<N, 1> || Vector = Matrix<1, N>
template<typename Derived, typename T, std::size_t Dimensions>
template<std::size_t Rows, std::size_t Cols>
constexpr vec_core<Derived, T, Dimensions>::vec_core(const matrix<T, Rows, Cols>& mat) {
    static_assert(Rows == 1 || Cols == 1, "Either Rows or Cols has to be 1.");
    constexpr std::size_t matrix_size = (Rows > Cols) ? Rows : Cols;
    constexpr std::size_t min_size = (matrix_size < Dimensions) ? matrix_size : Dimensions;

    data.fill(0);
    for (std::size_t i = 0; i < min_size; ++i) {
        data[i] = mat[i];
    }
}

// Matrix = Vector
template<typename T, std::size_t Rows, std::size_t Cols>
template<typename Derived, std::size_t Dimensions>
constexpr matrix<T, Rows, Cols>::matrix(const vec_core<Derived, T, Dimensions>& v) {
    static_assert(Rows == 1 || Cols == 1, "Either Rows or Cols has to be 1.");
    constexpr std::size_t matrix_size = (Rows > Cols) ? Rows : Cols;
    constexpr std::size_t min_size = (matrix_size < Dimensions) ? matrix_size : Dimensions;

    data.fill(0);
    for (std::size_t i = 0; i < min_size; ++i) {
        data[i] = v[i];
    }
}


// ------------------------------
// Operator Functions
// ------------------------------

// Matrix * Vector
template<typename T, std::size_t Rows, std::size_t Cols, typename VecDerived>
constexpr matrix<T, Rows, 1> operator*(const matrix<T, Rows, Cols>& mat, const vec_core<VecDerived, T, Cols>& vector) noexcept {
    matrix<T, Rows, 1>result{};
    for (std::size_t r = 0; r < Rows; ++r) {
        for (std::size_t c = 0; c < Cols; ++c) {
            result[r] += mat(r, c) * vector[c];
        }
    }
    return result;
}

// Vector * Matrix
template<typename T, std::size_t Rows, std::size_t Cols, typename VecDerived>
constexpr matrix<T, 1, Cols> operator*(const vec_core<VecDerived, T, Rows>& vector, const matrix<T, Rows, Cols> mat) noexcept {
    matrix<T, 1, Cols> result{};
    for (std::size_t c = 0; c < Cols; ++c) {
        for (std::size_t r = 0; r < Rows; ++r) {
            result[c] += vector[r] * mat(r, c);
        }
    }
    return result;
}

// Vector *= Matrix
template<typename T, std::size_t Rows, std::size_t Cols, typename VecDerived>
constexpr vec_core<VecDerived, T, Cols>& operator*=(vec_core<VecDerived, T, Rows>& vector, const matrix<T, Rows, Cols>& mat) noexcept {
    static_assert(Rows == Cols, "In-place multiplication requires square matrix");
    vec_core<VecDerived, T, Cols> result{};
    for (std::size_t c = 0; c < Cols; ++c) {
        for (std::size_t r = 0; r < Rows; ++r) {
            result[c] += vector[r] * mat(r, c);
        }
    }
    for (std::size_t i = 0; i < Cols; ++i) {
        vector[i] = result[i];
    }
    return vector;
}

#endif