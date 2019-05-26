#pragma once

#include<complex>
#include<vector>

using ThomasVectorType = std::vector<std::complex<double>>;

struct TridiagonalSystem {
	TridiagonalSystem(size_t N) :
		a(N-1), b(N), c(N-1), d(N)
	{}
	ThomasVectorType a;
	ThomasVectorType b;
	ThomasVectorType c;
	ThomasVectorType d;
};

void Thomas(ThomasVectorType & a,
		ThomasVectorType & b,		
		ThomasVectorType & c,		
		ThomasVectorType & d);
void Thomas(TridiagonalSystem &);
