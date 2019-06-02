#include "lcp-app.hpp"

void LcpApp::SaveGridFtorToFile(comptools::Grid g, const std::function<double(double)> & f, std::string filename) const {
		std::FILE * fp = std::fopen(filename.c_str(), "w");
		for (auto x : g) {
			fprintf(fp, "%.10f %.10f\n", x, f(x));
		}
		std::fclose(fp);
}

void LcpApp::SaveFunctionToFile(ComplexFunctionType & f, std::string filename) const {
		std::FILE * fp = std::fopen(filename.c_str(), "w");
		for (auto p : f) {
			fprintf(fp, "%.10f %.10f %.10f\n", p.x, p.y.real(), p.y.imag());
		}
		std::fclose(fp);
}

void LcpApp::SaveFunctionToFile(RealFunctionType & f, std::string filename) const {
		std::FILE * fp = std::fopen(filename.c_str(), "w");
		for (auto p : f) {
			fprintf(fp, "%.10f %.10f\n", p.x, p.y);
		}
		std::fclose(fp);
}

void LcpApp::GetVibrationalState(size_t nu) {
	auto V = comptools::LoadFunctionFromFile(settings_.V0file);
	auto Vinterpolated = comptools::interpol::NaturalSplineInterpol(V);
	
	Parameters p;
	p.a = settings_.a;
	p.b = settings_.b;
	p.mu = settings_.mu;
	p.basisSize = settings_.basisSize;
	p.V = std::bind(Vinterpolated, std::placeholders::_1);
	
	SchrodingerBox box(p);
	
	rr_.Enu = box.Energies()(nu);

	for (auto x : xgrid_) {
		ComplexFunctionType::Point p = { x, box.wf(0, x) };
		functions_.chi.push_back(p);
	}
}

void LcpApp::PrecalculateFunctions() {
	GetVibrationalState(settings_.nu);

	RealFunctionType vres;
	RealFunctionType gamma;
	io::CSVReader<3, io::trim_chars<' ', '\t'>, io::no_quote_escape<' '>, io::single_line_comment<'#'> > in(settings_.Vlocfile);
	double x;
	double v;
	double g;
	while(in.read_row(x, v, g)) {
		comptools::Function<double,double>::Point vp = {x, v};
		comptools::Function<double,double>::Point gp = {x, g};
		vres.push_back(vp);
		gamma.push_back(gp);
	}

	auto gammaInterpolated = comptools::interpol::NaturalSplineInterpol(gamma);
	auto vresInterpolated = comptools::interpol::NaturalSplineInterpol(vres);

	for (auto x : xgrid_) {
		RealFunctionType::Point p = { x, vresInterpolated(x) - rr_.Enu - settings_.ea };
		functions_.vres.push_back(p);
		p.y = gammaInterpolated(x);
		if (p.y < 0.0)
			p.y = 0.0;
		functions_.gamma.push_back(p);
	}
	rr_.VresAs = functions_.vres[xgrid_.size() - 1].y;

}

TridiagonalSystem LcpApp::GenerateTridiagonalSystem() {
	size_t N = ecsgrid_.size();
	size_t ic = ecsgrid_.ic();
	double h = ecsgrid_.h;

	std::complex<double> phase(cos(ecsgrid_.eta), sin(ecsgrid_.eta));
	std::complex<double> ii(0.0, 1.0);
	auto OneOverPhaseSquared = 1.0 / phase / phase;
	auto OneOverOnePlusPhase = 1.0 / (phase + 1.0);
	std::complex<double> kinfac = 0.5 / h / h / settings_.mu;

	TridiagonalSystem s(N);

	s.b[0] = kinfac - functions_.vres.y(0) + functions_.gamma.y(0) * ii;
	s.c[0] = kinfac;
	for (size_t i = 1; i < ic; ++i) {
		s.a[i-1] = kinfac;
		s.b[i] = -2.0 *kinfac - functions_.vres.y(i) + 0.5 * functions_.gamma.y(i) * ii;
		s.c[i] = kinfac;
	}

	s.a[ic-1] = 2.0 * kinfac * OneOverOnePlusPhase;
	s.b[ic] = -2.0 * kinfac / phase - functions_.vres.y(ic) + 0.5 * functions_.gamma.y(ic) * ii;
	s.c[ic] = 2.0 * kinfac * OneOverOnePlusPhase / phase;

	kinfac *= OneOverPhaseSquared;
	for (size_t i = ic + 1; i < N-1; ++i) {
		s.a[i-1] = kinfac;
		s.b[i] = -2.0 * kinfac - functions_.vres.y(i) + 0.5 * functions_.gamma.y(i) * ii;
		s.c[i] = kinfac;
	}

	s.a[N-2] = kinfac;
	s.b[N-1] =-2.0 * kinfac - functions_.vres.y(N-1) + 0.5 * functions_.gamma.y(N-1) * ii;

	for (size_t i = 0; i < N; ++i) {
		s.d[i] = sqrt(0.5*functions_.gamma.y(i) / M_PI) * functions_.chi.y(i);
	}
	return s;
}

void LcpApp::UpdateEnergyTerm(double E, TridiagonalSystem & s) {
	for (auto & x : s.b)
		x += E;
}
