#include <iostream>
#include <cstdio>
#include <functional>
#include <string>

#include <comptools/grid.hpp>

#include "domcke-model.hpp"
#include "output.hpp"

using namespace comptools;


int main() {
	DomckeModelUnitWrapper::UnitConversionFactors u;
	u.DToNew = 1.8897259886;
	u.DToOld = 1.0 / u.DToNew;
	u.EToOld = 27.211386245988;
	u.EToNew = 1.0 / u.EToOld;

	DomckeModel::Parameters p;
	p.al = 1.25;
	p.A = 1.757;
	p.B = 1.667;
	p.C = 0.98;
	p.mu = 1.0;
	p.a = 1.96;
	p.D0 = 5.0;
	p.D1 = 3.2;
	p.R0 = 2.5 * u.DToOld;
	p.t = 0.2;
	p.Qd = 1.2;
	p.bisa = -50;
	p.bisb = 20;


	auto RGrid = Grid::FromMinMaxN(1.5, 10, 1000);
	auto RGrid2 = Grid::FromMinMaxN(1.5, 10, 1000);
	auto EGrid = Grid::FromMinMaxN(-10, 15, 1500);


	DomckeModelUnitWrapper model(u, p);
	OutputHandler<BasicPrinter> output;
	output.PrintGridFunction(RGrid, [&model](double R) { return model.V0(R); }, "V0.dat");
	output.PrintGridFunction(RGrid, [&model](double R) { return model.Vd(R); }, "Vd.dat");
	output.PrintGridFunction(RGrid2, [&model](double R) { return model.Vres(R); }, "Vres.dat");
	output.PrintGridFunction(RGrid2, [&model](double R) { return model.Vloc(R); }, "Vloc.dat");
	
	output.PrintGridFunction(EGrid, [&model](double E) { return model.Gamma(E, 2.5);}, "gamma.dat");
	output.PrintGridFunction(EGrid, [&model](double E) { return model.Delta(E, 2.5);}, "delta.dat");
	output.PrintGridFunction(EGrid, [&model](double E) { return model.EresEquation(E, 2.5);}, "eresEquation.dat");
	return 0;
}
