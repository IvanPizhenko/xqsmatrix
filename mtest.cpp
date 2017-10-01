#include <iostream>
#include <random>
#include "xqsmatrix.h"

int main()
{
	std::mt19937 rng(111);
	std::uniform_real_distribution<double> d(0, 1.0);
	
	constexpr std::size_t N = 10;
	XQSMatrix<double> m(N, N);
	for (std::size_t i = 0; i < N; ++i) {
		for (std::size_t j = 0; j < N; ++j) {
			m(i, j) = d(rng);
		}
	} 

	const auto mm1 = m.inverse_v1();
	const auto pr1 = m * mm1;

	const auto mm2 = m.inverse_v2();
	const auto pr2 = m * mm2;

	std::cout << pr1 << std::endl << std::endl;
	std::cout << pr2 << std::endl << std::endl;

	return 0;
}