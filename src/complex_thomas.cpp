#include "complex_thomas.hpp"

void Thomas(ThomasVectorType & a,
		ThomasVectorType & b,		
		ThomasVectorType & c,		
		ThomasVectorType & d) {
	/*
	| b c . . . . | |x| = |d|
	| a b c . . . | |x| = |d|
	| . a b c . . | |x| = |d|
	.                      .
	.                      .
	| . . . . a b | |x| = |d|
	*/

    auto n = d.size();

    for (size_t i = 1; i < n; ++i) {
		std::complex<double> w = a[i-1] / b[i-1];
        b[i] -= w * c[i-1];
        d[i] -= w * d[i-1];
    }

    b[n-1] = d[n-1] / b[n-1];
    for (int i = n - 2; i >= 0; --i) {
        b[i] = (d[i] - c[i] * b[i+1]) / b[i];
    }
}

void Thomas(TridiagonalSystem & s) { Thomas(s.a, s.b, s.c, s.d); }
