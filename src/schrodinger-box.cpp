#include "schrodinger-box.hpp"

void SchrodingerBox::Run() {
	MatrixType hmat(p.basisSize, p.basisSize);
	HamMatPotential(hmat);
	HamMatKinetic(hmat);
	auto D = Diagonalize(hmat);
	for (auto & x : D) 
		printf("%12.8f\n", x);
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
			printf("%12.8f ", m(i,j));
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

double SchrodingerBox::PontentialMatrixElement(size_t i, size_t j) {
	size_t dif = i - j;
	size_t sum = i + j;
	double Vdiag = -8.0 * i * j * (p.b - p.a) / (M_PI*M_PI*sum*sum*dif*dif);
	return Vdiag;
}

SchrodingerBox::VectorType SchrodingerBox::Diagonalize(MatrixType & A) {
	char JOBZ = 'V';
	char UPLO = 'U';
	int N = A.rows();
	VectorType W(N);
	int LWORK=3*N-1;
	std::vector<double> WORK(LWORK);
	int INFO;

	dsyev_(&JOBZ,
			&UPLO,
			&N,
			A.data(),
			&N,
			W.data(),
			WORK.data(),
			&LWORK,
			&INFO);

	return W;
}

void SchrodingerBox::EvaluatePotential(VectorType & v) {
	for (auto & x : v)
		x = p.V(x);
}

void SchrodingerBox::TransformBack(MatrixType & m, VectorType & D) {
	auto A = m;
	for (size_t i = 0; i < m.rows(); ++i) {
		for (size_t j = 0; j < m.cols(); ++j) {
			double temp = 0;
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
