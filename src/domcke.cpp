#include <iostream>
#include <cstdio>
#include <functional>
#include <string>

#include <comptools/grid.hpp>

#include "domcke-model.hpp"

using namespace comptools;

void PrintRow(std::FILE * fp, double x, double y) {
		fprintf(fp, "%.10f %.10f\n", x, y);
}

void PrintRow(std::FILE * fp, double x, std::complex<double> y) {
		fprintf(fp, "%.10f %.10f %.10f\n", x, y.real(), y.imag());
}

template<class FunctionType>
void PrintGridFunction(Grid g, const FunctionType & f, std::string filename) {
	std::FILE * fp = std::fopen(filename.c_str(), "w");
	for (auto x : g) {
		PrintRow(fp, x, f(x));
	}
	std::fclose(fp);
}

int main() {
	DomckeModel::Parameters p;
	p.al = 1.25;
	p.A = 1.757;
	p.B = 1.667;
	p.C = 0.98;
	p.mu = 1.0;
	p.a = 1.96;
	p.D0 = 5.0;
	p.D1 = 3.2;
	p.R0 = 2.5;
	p.t = 0.2;
	p.Qd = 1.2;
	p.bisa = -5;
	p.bisb = 15;

	auto RGrid = Grid::FromMinMaxN(1.5, 10, 1000);
	auto RGrid2 = Grid::FromMinMaxN(2.1, 5, 1000);
	auto EGrid = Grid::FromMinMaxN(-10, 15, 1500);


	DomckeModel model(p);
	PrintGridFunction(RGrid, [&model](double R) { return model.V0(R); }, "V0.dat");
	PrintGridFunction(RGrid, [&model](double R) { return model.Vd(R); }, "Vd.dat");
	
	PrintGridFunction(EGrid, [&model](double E) { return model.Gamma(E, 2.5);}, "gamma.dat");
	PrintGridFunction(EGrid, [&model](double E) { return model.Delta(E, 2.5);}, "delta.dat");
	PrintGridFunction(EGrid, [&model](double E) { return model.EresEquation(E, 2.5);}, "eresEquation.dat");
	PrintGridFunction(RGrid2, [&model](double R) { return model.Vres(R); }, "Vres.dat");
	PrintGridFunction(RGrid2, [&model](double R) { return model.Vloc(R); }, "Vloc.dat");
	return 0;
}
