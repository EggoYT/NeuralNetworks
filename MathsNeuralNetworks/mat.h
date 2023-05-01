#ifndef MAT_H
#define MAT_H

template<size_t cols, size_t rows>
class Matrix {
private:
	union {
		frac data[cols * rows];
		//Vector<cols> data_vec[rows];
	};
public:
	template<size_t, size_t> friend class Matrix;
	template<size_t> friend class Vector;
	constexpr Matrix();
	constexpr Matrix(const frac& value);
	constexpr Matrix(const Matrix& another);
	template<size_t N> constexpr Matrix(const frac(&init)[N]);

	constexpr Matrix<cols, rows>& operator+=(const Vector<cols>& vec) noexcept;
	constexpr Matrix<cols, rows>& operator-=(const Vector<cols>& vec) noexcept;
	constexpr Matrix<cols, rows>& operator*=(const Vector<cols>& vec) noexcept;
	constexpr Matrix<cols, rows>& operator/=(const Vector<cols>& vec) noexcept;

	constexpr Vector<rows> operator+(const Vector<cols>& another) noexcept;
	constexpr Vector<rows> operator-(const Vector<cols>& another) noexcept;
	constexpr Vector<rows> operator*(const Vector<cols>& another) noexcept;
	constexpr Vector<rows> operator/(const Vector<cols>& another) noexcept;

	constexpr Matrix<cols, rows>& operator=(const Matrix& another) noexcept;
	constexpr Matrix<cols, rows>& operator+=(const Matrix& another) noexcept;
	constexpr Matrix<cols, rows>& operator-=(const Matrix& another) noexcept;
	constexpr Matrix<cols, rows>& operator*=(const Matrix& another) noexcept;
	constexpr Matrix<cols, rows>& operator/=(const Matrix& another) noexcept;

	constexpr Matrix<cols, rows> operator+(const Matrix<cols, rows>& another) noexcept;
	constexpr Matrix<cols, rows> operator-(const Matrix<cols, rows>& another) noexcept;
	template<size_t p> constexpr Matrix<cols, p> operator*(const Matrix<rows, p>& another) noexcept;
	
	frac& operator()(const size_t& col, const size_t& row) noexcept;
	//constexpr Vector<cols>& operator[](const size_t& row) noexcept;

	constexpr void sum(const Matrix<cols, rows>& first, const Matrix<cols, rows>& second) noexcept;
	constexpr void sub(const Matrix<cols, rows>& first, const Matrix<cols, rows>& second) noexcept;
	template<size_t p> constexpr void mult(const Matrix<cols, p>& first, const Matrix<p, rows>& second) noexcept;

	constexpr Matrix<rows, cols> transpose() noexcept;

	constexpr size_t columns() const noexcept;
	constexpr size_t lines() const noexcept;

	constexpr void set_null() noexcept;

	template<size_t m, size_t n>
	constexpr Matrix<cols - m + 1, rows - n + 1> operator&(const Matrix<m, n> kernel);

	~Matrix();
};
#endif