#include <initializer_list>
#include "defines.h"
#include "math.h"
#include "mat.h"
#include "vec.h"

template<size_t size>
inline constexpr Vector<size>::Vector() : Vector<size>::Vector(frac()) {}

template<size_t size>
inline constexpr Vector<size>::Vector(const frac& filler) {
	for (frac& i : this->data)
		i = filler; 
}

template<size_t size>
inline constexpr Vector<size>::Vector(const Vector& another) {
	for (size_t i = 0; i < size; ++i)
		data[i] = another.data[i];
}

template<size_t size>
template<size_t p>
inline constexpr Vector<size>::Vector(const frac(&filler)[p]) {
	static_assert(p == size, "Wrong init list size!");
	for (size_t i = 0; i < size; ++i) {
		data[i] = filler[i];
	}
}

template<size_t size>
inline constexpr Vector<size>& Vector<size>::operator=(const Vector& another) noexcept {
	for (size_t i = 0; i < size; ++i)
		this->data[i] = another.data[i];
	return *this;
}

template<size_t size>
inline constexpr Vector<size>& Vector<size>::operator+=(const Vector& another) noexcept {
	for (size_t i = 0; i < size; ++i)
		this->data[i] += another.data[i];
	return *this;
}

template<size_t size>
inline constexpr Vector<size>& Vector<size>::operator-=(const Vector& another) noexcept {
	for (size_t i = 0; i < size; ++i)
		this->data[i] -= another.data[i];
	return *this;
}

template<size_t size>
inline constexpr Vector<size>& Vector<size>::operator*=(const Vector& another) noexcept {
	for (size_t i = 0; i < size; ++i)
		this->data[i] *= another.data[i];
	return *this;
}

template<size_t size>
inline constexpr Vector<size>& Vector<size>::operator/=(const Vector& another) noexcept {
	for (size_t i = 0; i < size; ++i)
		this->data[i] /= another.data[i];
	return *this;
}

template<size_t size>
inline constexpr Vector<size> Vector<size>::operator+(const Vector<size>& another) noexcept {
	Vector<size> result;
	for (size_t i = 0; i < size; ++i)
		result.data[i] = data[i] + another.data[i];
	return result;
}

template<size_t size>
inline constexpr Vector<size> Vector<size>::operator-(const Vector& another) noexcept {
	Vector<size> result;
	for (size_t i = 0; i < size; ++i)
		result.data[i] = data[i] - another.data[i];
	return result;
}

template<size_t size>
inline constexpr Vector<size> Vector<size>::operator*(const Vector& another) noexcept {
	Vector<size> result;
	for (size_t i = 0; i < size; ++i)
		result.data[i] = data[i] * another.data[i];
	return result;
}

template<size_t size>
inline constexpr Vector<size> Vector<size>::operator/(const Vector& another) noexcept {
	Vector<size> result;
	for (size_t i = 0; i < size; ++i)
		result.data[i] = data[i] / another.data[i];
	return result;
}

template<size_t size>
inline constexpr frac& Vector<size>::operator[](const size_t& index) noexcept {
	return data[index];
}

template<size_t size>
inline constexpr frac& Vector<size>::operator()(const size_t& index) noexcept {
	return data[index];
}

template<size_t size>
inline constexpr void Vector<size>::sum(const Vector& first, const Vector& second) {
	for (size_t i = 0; i < size; ++i)
		data[i] = first.data[i] + second.data[i];
}

template<size_t size>
inline constexpr void Vector<size>::sub(const Vector& first, const Vector& second) {
	for (size_t i = 0; i < size; ++i)
		data[i] = first.data[i] - second.data[i];
}

template<size_t size>
inline constexpr void Vector<size>::mult(const Vector& first, const Vector& second) {
	for (size_t i = 0; i < size; ++i)
		data[i] = first.data[i] * second.data[i];
}

template<size_t size>
inline constexpr void Vector<size>::div(const Vector& first, const Vector& second) {
	for (size_t i = 0; i < size; ++i)
		data[i] = first.data[i] / second.data[i];
}

template<size_t size>
template<size_t p>
inline constexpr void Vector<size>::mult(const Matrix<p, size>& first, const Vector<p>& second) {
	for (size_t i = 0; i < size; ++i) {
		const frac* f = first.data + i * p;
		frac& d = data[i];
		d = 0.0;
		for (size_t j =  0; j < p; ++j)
			d += f[j] * second.data[j];
	}
}

template<size_t size>
Vector<size>::~Vector() {}

