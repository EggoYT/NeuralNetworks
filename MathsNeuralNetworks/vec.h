#ifndef VEC_H
#define VEC_H

template <size_t size>
class Vector {
private:
	frac data[size]{};
public:
	template<size_t, size_t> friend class Matrix;
	template<size_t> friend class Vector;

	constexpr Vector();
	constexpr Vector(const frac& filler);
	constexpr Vector(const Vector& another);
	template<size_t p> constexpr Vector(const frac (&filler)[p]);

	constexpr Vector<size>& operator=(const Vector&) noexcept;
	constexpr Vector<size>& operator+=(const Vector&) noexcept;
	constexpr Vector<size>& operator-=(const Vector&) noexcept;
	constexpr Vector<size>& operator*=(const Vector&) noexcept;
	constexpr Vector<size>& operator/=(const Vector&) noexcept;

	constexpr Vector<size> operator+(const Vector& another) noexcept;
	constexpr Vector<size> operator-(const Vector& another) noexcept;
	constexpr Vector<size> operator*(const Vector& another) noexcept;
	constexpr Vector<size> operator/(const Vector& another) noexcept;

	constexpr frac& operator[](const size_t& index) noexcept;
	constexpr frac& operator()(const size_t& index) noexcept;

	constexpr void sum(const Vector& first, const Vector& second);
	constexpr void sub(const Vector& first, const Vector& second);
	constexpr void mult(const Vector& first, const Vector& second);
	constexpr void div(const Vector& first, const Vector& second);

	template<size_t p> constexpr void mult(const Matrix<p, size>& first, const Vector<p>& second);

	~Vector();
};
#endif