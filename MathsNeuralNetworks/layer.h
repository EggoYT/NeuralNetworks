#ifndef LAYER_H
#define LAYER_H

template<size_t In, size_t Out, typename Type, typename = std::enable_if_t<In != 0>>
class Layer {
private:
	Matrix<In, Out> data;
public:
	Layer(){}

	//constexpr void calculate(const Vector<In>& in, Vector<Out>& out);
};
#endif