#include "schrodinger-box.hpp"

void SchrodingerBox::Run() {
	hmat = MatrixType(p.basisSize, p.basisSize);
	HamMatPotential(hmat);
	HamMatKinetic(hmat);
	energies = Diagonalize(hmat);
}

void SchrodingerBox::HamMatPotential(MatrixType & m) {
	FillWithXMean(m);
	auto D = Diagonalize(m);
	EvaluatePotential(D);
	TransformBack(m, D);
}

void SchrodingerBox::printMatrix(MatrixType & m) {
	for (size_t i = 0; i < m.rows(); ++i) {
		for (size_t j = 0; j < m.cols(); ++j) {
			printf("%12.8f + %12.8fi", m(i,j).real(), m(i,j).imag());
		}
		printf("\n");
	}
}

void SchrodingerBox::FillWithXMean(MatrixType & m) {
	for (size_t i = 0; i < m.rows(); ++i) {
		for (size_t j = 0; j < i; ++j) {
			if ((i+j+1)%2 == 0) {
				m(i, j) = PontentialMatrixElement(i+1 , j+1);
			} else {
				m(i,j) = 0;
			}
		}
		m(i,i) = 0;
	}
}

SchrodingerBox::NumberType SchrodingerBox::PontentialMatrixElement(size_t i, size_t j) {
	size_t dif = i - j;
	size_t sum = i + j;
	double Vdiag = -8.0 * i * j * (p.b - p.a) / (M_PI*M_PI*sum*sum*dif*dif);
	return Vdiag;
}

SchrodingerBox::RealVectorType SchrodingerBox::Diagonalize(MatrixType & A) {
	char JOBZ = 'V';
	char UPLO = 'U';
	int N = A.rows();
	RealVectorType W(N);
	RealVectorType RWORK(3*N-2);
	int LWORK=2*N-1;
	VectorType WORK(LWORK);
	int INFO;

	zheev_(&JOBZ,
			&UPLO,
			&N,
			A.data(),
			&N,
			W.data(),
			WORK.data(),
			&LWORK,
			RWORK.data(),
			&INFO);

	return W;
}

void SchrodingerBox::EvaluatePotential(SchrodingerBox::RealVectorType & v) {
	for (auto & x : v) {
		// printf("%10f", x);
		x = p.V(x);
		// printf(" %10f\n", x);
	}
}

void SchrodingerBox::TransformBack(MatrixType & m, RealVectorType & D) {
	auto A = m;
	for (size_t i = 0; i < m.rows(); ++i) {
		for (size_t j = 0; j < m.cols(); ++j) {
			NumberType temp = 0;
			for (size_t k = 0; k < m.cols(); ++k) {
				temp += A(k,i) * A(k,j) * D(k);
			}
			m(i,j) = temp;
		}
	}
}

void SchrodingerBox::HamMatKinetic(MatrixType & m) {
	double f1 = M_PI / (p.b - p.a);
	for (size_t i = 0; i < m.rows(); ++i) {
		double f = (i+1) * f1;
		m(i, i) += f * f / 2 / p.mu;
	}
}

SchrodingerBox::NumberType SchrodingerBox::wf(size_t i, double x) {
	NumberType y  = 0;
	for (size_t j = 0; j < p.basisSize; ++j) {
		y += hmat(i, j) * basis(j, x);
	}
	return y;
}
