#include <initializer_list>
#include <cassert>
#include "defines.h"
#include "math.h"
#include "vec.h"
#include "mat.h"

template<size_t cols, size_t rows>
inline constexpr Matrix<cols, rows>::Matrix() : Matrix<cols, rows>::Matrix(frac()) {
}

template<size_t cols, size_t rows>
inline constexpr Matrix<cols, rows>::Matrix(const frac& value) {
	for (size_t i = 0; i < cols * rows; ++i)
		data[i] = value;
}

template<size_t cols, size_t rows>
inline constexpr Matrix<cols, rows>::Matrix(const Matrix& another) {
	for (size_t i = 0; i < cols * rows; ++i)
		data[i] = another.data[i];
}

template<size_t cols, size_t rows>
template<size_t N>
inline constexpr Matrix<cols, rows>::Matrix(const frac(&init)[N]) {
	static_assert(N == cols * rows, "Wrong init list size!");
	for (size_t i = 0; i < N; ++i)
		data[i] = init[i];
}

template<size_t cols, size_t rows>
inline constexpr Matrix<cols, rows>& Matrix<cols, rows>::operator+=(const Vector<cols>& vec) noexcept {
	for (size_t i = 0; i < cols; ++i) {
		const frac& v = vec[i];
		frac* d = data + i;
		for (size_t j = 0; j < rows; ++j)
			data[j * cols] += v;
	}
	return *this;
}

template<size_t cols, size_t rows>
inline constexpr Matrix<cols, rows>& Matrix<cols, rows>::operator-=(const Vector<cols>& vec) noexcept {
	for (size_t i = 0; i < cols; ++i) {
		const frac& v = vec[i];
		frac* d = data + i;
		for (size_t j = 0; j < rows; ++j)
			data[j * cols] -= v;
	}
	return *this;
}

template<size_t cols, size_t rows>
inline constexpr Matrix<cols, rows>& Matrix<cols, rows>::operator*=(const Vector<cols>& vec) noexcept {
	for (size_t i = 0; i < cols; ++i) {
		const frac& v = vec[i];
		frac* d = data + i;
		for (size_t j = 0; j < rows; ++j)
			data[j * cols] *= v;
	}
	return *this;
}

template<size_t cols, size_t rows>
inline constexpr Matrix<cols, rows>& Matrix<cols, rows>::operator/=(const Vector<cols>& vec) noexcept {
	for (size_t i = 0; i < cols; ++i) {
		const frac& v = vec[i];
		frac* d = data + i;
		for (size_t j = 0; j < rows; ++j)
			data[j * cols] /= v;
	}
	return *this;
}

template<size_t cols, size_t rows>
inline constexpr Vector<rows> Matrix<cols, rows>::operator+(const Vector<cols>& another) noexcept {
	Vector<rows> result(0.0);

	for (size_t i = 0; i < rows; ++i) {
		frac& r = result[i];
		const frac* d = data + i * cols;
		for (size_t j = 0; j < cols; ++j) {
			r += d[j] + another.data[j];
		}
	}

	return result;
}

template<size_t cols, size_t rows>
inline constexpr Vector<rows> Matrix<cols, rows>::operator-(const Vector<cols>& another) noexcept {
	Vector<rows> result(0.0);

	for (size_t i = 0; i < rows; ++i) {
		frac& r = result[i];
		const frac* d = data + i * cols;
		for (size_t j = 0; j < cols; ++j) {
			r += d[j] - another.data[j];
		}
	}

	return result;
}

template<size_t cols, size_t rows>
inline constexpr Vector<rows> Matrix<cols, rows>::operator*(const Vector<cols>& another) noexcept {
	Vector<rows> result(0.0);

	for (size_t i = 0; i < rows; ++i) {
		frac& r = result[i];
		const frac* d = data + i * cols;
		for (size_t j = 0; j < cols; ++j) {
			r += d[j] * another.data[j];
		}
	}

	return result;
}

template<size_t cols, size_t rows>
inline constexpr Vector<rows> Matrix<cols, rows>::operator/(const Vector<cols>& another) noexcept {
	Vector<rows> result(0.0);

	for (size_t i = 0; i < rows; ++i) {
		frac& r = result[i];
		const frac* d = data + i * cols;
		for (size_t j = 0; j < cols; ++j) {
			r += d[j] / another.data[j];
		}
	}

	return result;
}

template<size_t cols, size_t rows>
inline constexpr Matrix<cols, rows>& Matrix<cols, rows>::operator=(const Matrix<cols, rows>& another) noexcept {
	for (size_t i = 0; i < rows; ++i) {
		frac* a1 = data + i * cols;
		const frac* a2 = another.data + i * cols;
		for (size_t j = 0; j < cols; ++j)
			a1[j] = a2[j];
	}
	return *this;
}

template<size_t cols, size_t rows>
inline constexpr Matrix<cols, rows>& Matrix<cols, rows>::operator+=(const Matrix<cols, rows>& another) noexcept {
	for (size_t i = 0; i < rows; ++i) {
		frac* d = data + i * cols;
		const frac* a = another.data + i * cols;
		for (size_t j = 0; j < cols; ++j) {
			d[j] += a[j];
		}
	}
	return *this;
}

template<size_t cols, size_t rows>
inline constexpr Matrix<cols, rows>& Matrix<cols, rows>::operator-=(const Matrix<cols, rows>& another) noexcept {
	for (size_t i = 0; i < rows; ++i) {
		frac* d = data + i * cols;
		const frac* a = another.data + i * cols;
		for (size_t j = 0; j < cols; ++j) {
			d[j] -= a[j];
		}
	}
	return *this;
}

template<size_t cols, size_t rows>
inline constexpr Matrix<cols, rows>& Matrix<cols, rows>::operator*=(const Matrix<cols, rows>& another) noexcept {
	for (size_t i = 0; i < rows; ++i) {
		frac* d = data + i * cols;
		const frac* a = another.data + i * cols;
		for (size_t j = 0; j < cols; ++j) {
			d[j] *= a[j];
		}
	}
	return *this;
}

template<size_t cols, size_t rows>
inline constexpr Matrix<cols, rows>& Matrix<cols, rows>::operator/=(const Matrix<cols, rows>& another) noexcept {
	for (size_t i = 0; i < rows; ++i) {
		frac* d = data + i * cols;
		const frac* a = another.data + i * cols;
		for (size_t j = 0; j < cols; ++j) {
			d[j] /= a[j];
		}
	}
	return *this;
}

template<size_t cols, size_t rows>
inline frac& Matrix<cols, rows>::operator()(const size_t& col, const size_t& row) noexcept {
	assert(col < cols&& row < rows);
	return data[col + row * cols];
}

//template<size_t cols, size_t rows>
//inline constexpr Vector<cols>& Matrix<cols, rows>::operator[](const size_t& row) noexcept {
//	assert(row < rows);
//	return data_vec[row];
//}

template<size_t cols, size_t rows>
inline constexpr Matrix<cols, rows> Matrix<cols, rows>::operator+(const Matrix<cols, rows>& another) noexcept {
	Matrix<cols, rows> result;
	for (size_t i = 0; i < rows; ++i) {
		const frac* d = data + i * cols;
		const frac* a = another.data + i * cols;
		frac* r = result.data + i * cols;
		for (size_t j = 0; j < cols; ++j) {
			r[j] = d[j] + a[j];
		}
	}
	return result;
}

template<size_t cols, size_t rows>
inline constexpr Matrix<cols, rows> Matrix<cols, rows>::operator-(const Matrix<cols, rows>& another) noexcept {
	Matrix<cols, rows> result;
	for (size_t i = 0; i < rows; ++i) {
		const frac* d = data + i * cols;
		const frac* a = another.data + i * cols;
		frac* r = result.data + i * cols;
		for (size_t j = 0; j < cols; ++j) {
			r[j] = d[j] - a[j];
		}
	}
	return result;
}

template<size_t cols, size_t rows>
template<size_t p>
inline constexpr Matrix<cols, p> Matrix<cols, rows>::operator*(const Matrix<rows, p>& another) noexcept {
	Matrix<cols, p> matrix(0.0);

	for (size_t i = 0; i < p; ++i) {
		const frac* a_row = another.data + i * rows;
		frac* r_row = matrix.data + i * cols;

		for (size_t j = 0; j < rows; ++j) {
			const frac& a = a_row[j];
			const frac* d_row = data + j * cols;

			for (int k = 0; k < cols; ++k)
				r_row[k] += d_row[k] * a;
		}
	}
	
	return matrix;
}

template<size_t cols, size_t rows>
inline constexpr void Matrix<cols, rows>::sum(const Matrix<cols, rows>& first, const Matrix<cols, rows>& second) noexcept {
	for (size_t i = 0; i < rows; ++i) {
		frac* d = data + i * cols;
		const frac* f = first.data + i * cols;
		const frac* s = second.data + i * cols;
		for (size_t j = 0; j < cols; ++j)
			d[j] = f[j] + s[j];
	}
}

template<size_t cols, size_t rows>
inline constexpr void Matrix<cols, rows>::sub(const Matrix<cols, rows>& first, const Matrix<cols, rows>& second) noexcept {
	for (size_t i = 0; i < rows; ++i) {
		frac* d = data + i * cols;
		const frac* f = first.data + i * cols;
		const frac* s = second.data + i * cols;
		for (size_t j = 0; j < cols; ++j)
			d[j] = f[j] - s[j];
	}
}

template<size_t cols, size_t rows>
template<size_t p>
inline constexpr void Matrix<cols, rows>::mult(const Matrix<cols, p>& first, const Matrix<p, rows>& second) noexcept {
	this->set_null();
	for (size_t i = 0; i < rows; ++i) {
		const frac* a_row = second.data + i * p;
		frac* r_row = data + i * cols;

		for (size_t j = 0; j < p; ++j) {
			const frac& a = a_row[j];
			const frac* f_row = first.data + j * cols;

			for (int k = 0; k < cols; ++k)
				r_row[k] += f_row[k] * a;
		}
	}
}

template<size_t cols, size_t rows>
constexpr Matrix<rows, cols> Matrix<cols, rows>::transpose() noexcept {
	Matrix<rows, cols> result;
	for (size_t i = 0; i < cols; ++i) {
		frac* r = result + i * rows;
		for (size_t j = 0; j < rows; ++j) {
			r[j] = data[i + j * cols];
		}
	}
	return result;
}

template<size_t cols, size_t rows>
constexpr size_t Matrix<cols, rows>::columns() const noexcept {
	return cols;
}

template<size_t cols, size_t rows>
constexpr size_t Matrix<cols, rows>::lines() const noexcept {
	return rows;
}

template<size_t cols, size_t rows>
inline constexpr void Matrix<cols, rows>::set_null() noexcept {
	for (size_t i = 0; i < cols * rows; ++i)
		data[i] = frac();
}

template<size_t cols, size_t rows>
template<size_t m, size_t n>
constexpr Matrix<cols - m + 1, rows - n + 1> Matrix<cols, rows>::operator&(const Matrix<m, n> kernel) {
	constexpr int64_t m_size = cols - m + 1;
	constexpr int64_t n_size = rows - n + 1;
	static_assert(m_size > 0 && n_size > 0, "Wrong kernel size!");

	Matrix<m_size, n_size> result(0.0);

	for (size_t j = 0; j < n_size; ++j) {
		frac* r_col = result.data + j * m_size;
		const frac* d_col = data + j * cols;

		for (size_t i = 0; i < m_size; ++i) {
			frac* r = r_col + i;
			const frac* d = d_col + i;

			for (size_t jj = 0; jj < n; ++jj) {
				frac* r_col2 = r + jj * m_size;
				const frac* d_col2 = d + jj * cols;
				const frac* k_col = kernel.data + jj * m;
				
				for (size_t ii = 0; ii < m; ++ii)
					r_col2[ii] += d_col2[ii] * k_col[ii];
			}
		}
	}

	return result;
}

template<size_t cols, size_t rows>
Matrix<cols, rows>::~Matrix() {}