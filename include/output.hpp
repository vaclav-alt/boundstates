#pragma once

#include <comptools/grid.hpp>

class BasicPrinter {
protected:
	void Print(std::FILE * fp, double x, double y) {
		fprintf(fp, "%.10f %.10f\n", x, y);
	}

	void Print(std::FILE * fp, double x, std::complex<double> y) {
		fprintf(fp, "%.10f %.10f %.10f\n", x, y.real(), y.imag());
	}
};

template<typename PrintPolicy>
class OutputHandler : private PrintPolicy {
public:
	using PrintPolicy::Print;

	template<class FunctionType>
	void PrintGridFunction(comptools::Grid g, const FunctionType & f, std::string filename) {
		std::FILE * fp = std::fopen(filename.c_str(), "w");
		for (auto x : g) {
			Print(fp, x, f(x));
		}
		std::fclose(fp);
	}
};
