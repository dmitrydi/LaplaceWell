/*
 * gwell.cpp
 *
 *  Created on: 25 ����. 2021 �.
 *      Author: Dmitry_Di
 */

#include "gwell.h"

using namespace std;

LaplWell::LaplWell(): stehf_coefs(CalcStehf(NCOEF)) {};
LaplWell::~LaplWell() {};
double LaplWell::pwd(const double td) const {
	return InverseLaplace(pwd_lapl, td);
}

namespace Rectangular {
Well::Well(): LaplWell(), bess(false), dx(1./NSEG) {};
Well::~Well() {};

double Well::SEXP(const double y, const double e) const {
	double b = exp(-2.*y*e);
	return b/(1.-b);
};



void Well::fill_if1(const double u,
		const double ywd, const double yed,
		const double alpha,
		Eigen::MatrixXd& matrix) const {
	double squ = sqrt(u+alpha*alpha);
	double ans = 0.5*dx/squ;
	double dy = std::abs(ywd-ywd);
	double sumy = ywd+ywd;
	ans *= exp(-squ*(2.*yed-sumy))+exp(-squ*sumy)+exp(-squ*(2.*yed-dy))+exp(-squ*dy);
	ans *= (1+SEXP(yed, squ));
	matrix = Eigen::MatrixXd::Ones(2*NSEG, 2*NSEG)*ans;
}

void Well::vect_if1_yd(const double u,
		const double yd, const double ywd, const double yed,
		const double alpha,
		Eigen::VectorXd& buf) {
	double squ = sqrt(u+alpha*alpha);
	double ans = 0.5*dx/squ;
	double dy = std::abs(yd-ywd);
	double sumy = yd+ywd;
	ans *= exp(-squ*(2.*yed-sumy))+exp(-squ*sumy)+exp(-squ*(2.*yed-dy))+exp(-squ*dy);
	ans *= (1+SEXP(yed, squ));
	buf = Eigen::VectorXd::Ones(2*NSEG)*ans;
}

void Well::fill_if2e(const double u,
		const double xwd,
		const double xed, const double xede,
		const double ywd, const double yed,
		const double alpha,
		Eigen::MatrixXd& matrix, Eigen::VectorXd& buf) const {
	matrix = Eigen::MatrixXd::Zero(2*NSEG, 2*NSEG);
	const double term1 = PI/xed*(2*yed-(ywd+ywd));
	const double term2 = PI/xed*(2*yed-abs(ywd-ywd));
	const double term3 = PI/xed*(ywd+ywd);
	const double term4 = PI/xed*(2*yed+abs(ywd-ywd));
	const double dterm1 = 1.-exp(-term1);
	const double dterm2 = 1.-exp(-term2);
	const double dterm3 = 1.-exp(-term3);
	const double dterm4 = 1.-exp(-term4);
	double ek_term, ek_, sexp_, aydywd, kpiOxed, ydPywd, mmult, A, d, max_mat=0.;
	for (int k = 1; k <= KMAX; ++k) {
		ek_term = k*PI/xede;
		ek_ = sqrt(u + ek_term*ek_term + alpha*alpha);
		sexp_ = SEXP(yed, ek_);
		aydywd = abs(ywd-ywd); //!
		kpiOxed = k*PI/xed;
		ydPywd = ywd + ywd; //!
		mmult = 2./kpiOxed/ek_*((exp(-ek_*(2.*yed - ydPywd)) + exp(-ek_*ydPywd) + exp(-ek_*(2.*yed)))*(1 + sexp_)
				+ sexp_);
		for (int j = 0; j < 2*NSEG; ++j) {
			double x1 = -1.+j*dx;
			double x2 = x1 + dx;
			buf(j)  = mmult*sin(0.5*kpiOxed*(x2 - x1))*std::cos(0.5*kpiOxed*(2.*xwd + x1 + x2));
			if (abs(buf(j)) > max_mat) max_mat = abs(buf(j));
		}
		for (int i = 0; i < 2*NSEG; ++i) {
			double xd = (xwd-1.+(i+0.5)*dx);
			double row_mult = cos(kpiOxed*xd);
			for (int j = 0; j < 2*NSEG; ++j) {
				matrix(i,j) += row_mult*buf(j);
			}
		}
		A = 2*xed/PI/(1-exp(-2*ek_*yed));
		d = A*(exp(-k*term1)/dterm1+exp(-k*term2)/dterm2+exp(-k*term3)/dterm3 + exp(-k*term4)/dterm4);
		if (isnan(d)) d = 0.;
		if (k > KMIN && (abs(d) <= TINY || abs(max_mat) <= TINY || abs(d/max_mat) < SUM_EPS)) break;
	}
}

void Well::vect_if2e_yd(const double u, const double xd, const double xwd,
		const double xed, const double xede,
		const double yd, const double ywd, const double yed,
		const double alpha, Eigen::VectorXd& buf) const {
	buf = Eigen::VectorXd::Zero(2*NSEG);
	const double term1 = PI/xed*(2*yed-(yd+ywd));
	const double term2 = PI/xed*(2*yed-abs(yd-ywd));
	const double term3 = PI/xed*(yd+ywd);
	const double term4 = PI/xed*(2*yed+abs(yd-ywd));
	const double dterm1 = 1.-exp(-term1);
	const double dterm2 = 1.-exp(-term2);
	const double dterm3 = 1.-exp(-term3);
	const double dterm4 = 1.-exp(-term4);
	double ek_term, ek_, sexp_, aydywd, kpiOxed,ydPywd, mmult, A, d, max_mat=0., row_mult;
	for (int k = 1; k <= KMAX; ++k) {
		ek_term = k*PI/xede;
		ek_ = sqrt(u + ek_term*ek_term + alpha*alpha);
		sexp_ = SEXP(yed, ek_);
		aydywd = abs(yd-ywd); //!
		kpiOxed = k*PI/xed;
		ydPywd = yd + ywd; //!
		mmult =  2./kpiOxed/ek_*((std::exp(-ek_*(2.*yed - ydPywd)) + std::exp(-ek_*ydPywd) + std::exp(-ek_*(2.*yed-aydywd)))*(1 + sexp_)
				+ std::exp(-ek_*aydywd)*sexp_);
		row_mult = cos(kpiOxed*xd);
		for (int j = 0; j < 2*NSEG; ++j) {
			double x1 = -1.+j*dx;
			double x2 = x1 + dx;
			buf(j)  += row_mult*mmult*sin(0.5*kpiOxed*(x2 - x1))*std::cos(0.5*kpiOxed*(2.*xwd + x1 + x2));
			if (abs(buf(j)) > max_mat) max_mat = abs(buf(j));
		}
		A = 2*xed/PI/(1-exp(-2*ek_*yed));
		d = A*(exp(-k*term1)/dterm1+exp(-k*term2)/dterm2+exp(-k*term3)/dterm3 + exp(-k*term4)/dterm4);
		if (isnan(d)) d = 0.;
		if (k > KMIN && (abs(d) <= TINY || abs(max_mat) <= TINY || abs(d/max_mat) < SUM_EPS)) break;
	}
}

void Well::fill_i1f2h(const double u,
		const double xwd, const double xed, const double xede, const double alpha,
		Eigen::MatrixXd& matrix, Eigen::VectorXd& buf) const {
	for (int i: {0, 2*NSEG-1}) {
		double xd = xwd-1.+(i+0.5)*dx;
		vect_i1f2h(u, xd, xwd, xed, xede, alpha, buf);
		for (int j = 0; j < 2*NSEG; ++j)
			matrix(i,j) = buf(j);
	}
	for (int i = 1; i < 2*NSEG - 1; ++i) {
		for (int j = 0; j < 2*NSEG; ++j) {
			if (j >= i) {
				matrix(i,j) = matrix(0, j-i);
			}
			else {
				matrix(i,j) = matrix(2*NSEG-1, 2*NSEG-1-i+j);
			}
		}
	}
}

void Well::vect_i1f2h(const double u,
		const double xd, const double xwd, const double xed, const double xede,
		const double alpha,
		Eigen::VectorXd& buf) const {
	buf = Eigen::VectorXd::Zero(2*NSEG);
	const double squ = sqrt(u+alpha*alpha);
	int kmin = static_cast<int>(abs(0.5*(2./squ/xede-(-1.+(2*NSEG-1)*dx)/xed+(xwd-1+0.5*dx)/xed+xwd/xed)));
	double dmult = 8.*dx*exp(squ*xede/xed*abs((xwd-1+0.5*dx) + xwd - (-1.+(2*NSEG-1)*dx)))/1.-exp(-squ*2.*xede); //?
	double mult = xed/xede/squ*0.5*xede/PI;
	double t1, t2, x1, x2, elem, d;
	for (int k = 0; k <=KMAX; ++k) {
		for (double beta: {-1., 1.}) {
			for (int j = 0; j < 2*NSEG; ++j) {
				x1 = -1.+j*dx;
				x2 = x1 + dx;
				t1 = squ*xede*((xd)/xed+beta*xwd/xed-x2/xed-2.*k);
				t2 = squ*xede*((xd)/xed+beta*xwd/xed-x1/xed-2.*k);
				elem = mult*bess.abs_ik0ab(t1, t2);
				if (k > 0) {
					t1 = squ*xede*((xd)/xed+beta*xwd/xed-x2/xed+2.*k);
					t2 = squ*xede*((xd)/xed+beta*xwd/xed-x1/xed+2.*k);
					elem += mult*bess.abs_ik0ab(t1, t2);
				}
				buf(j) += elem;
			}
		}
		d = dmult*exp(-squ*2*k*xede);
		if (isnan(d)) d = 0.;
		if (k>kmin && (d <= TINY || buf(2*NSEG-1) <= TINY || abs(d/buf(2*NSEG-1)) < SUM_EPS)) break;
	}
};

void Well::vect_i1f2h_yd(const double u,
		const double xd, const double xwd, const double xed, const double xede,
		const double yd, const double ywd,
		const double alpha, Eigen::VectorXd& buf) const {
	double adyd = abs(ywd-yd);
	if (adyd < 1e-16)
	{
		vect_i1f2h(u, xd, xwd, xed, xede, alpha, buf);
	} else {
		buf = Eigen::VectorXd::Zero(2*NSEG);
		const double squ = sqrt(u+alpha*alpha);
		auto func = [this, adyd, squ](double x) {return this->bess.k0(squ*std::sqrt(x*x+adyd*adyd));};
		double t1, t2, x1, x2, elem;
		double mult = xed/xede*0.5*xede/PI;
		for (int k = 0; k <=KMAX; ++k) {
			for (double beta: {-1., 1.}) {
				for (int j = 0; j < 2*NSEG; ++j) {
					x1 = -1.+j*dx;
					x2 = x1 + dx;
					t1 = xede*((xd)/xed+beta*xwd/xed-x2/xed-2.*k);
					t2 = xede*((xd)/xed+beta*xwd/xed-x1/xed-2.*k);
					elem = mult*qromb(func, t1, t2, INT_EPS);
					if (k > 0) {
						t1 = xede*((xd)/xed+beta*xwd/xed-x2/xed+2.*k);
						t2 = xede*((xd)/xed+beta*xwd/xed-x1/xed+2.*k);
						elem += mult*qromb(func, t1, t2, INT_EPS);
					}
					buf(j) += elem;
				}
			}
			if ( elem <= TINY || buf(2*NSEG-1) <= TINY || abs(elem/buf(2*NSEG-1)) < SUM_EPS) break;
		}
	}
}

void Well::fill_i2f2h(const double u, const double ywd, const double alpha, Eigen::MatrixXd& matrix) const {
	double squ = sqrt(u+alpha*alpha);
	matrix = -0.5*exp(-squ*abs(ywd-ywd))/squ*(dx)*Eigen::MatrixXd::Ones(2*NSEG, 2*NSEG);
}

void Well::vect_i2f2h_yd(const double u, const double yd, const double ywd, const double alpha, Eigen::VectorXd& buf) const {
	double squ = sqrt(u+alpha*alpha);
	buf = -0.5*exp(-squ*abs(yd-ywd))/squ*(dx)*Eigen::VectorXd::Ones(2*NSEG);
}


Fracture::Fracture(const Boundary boundary, const double xwd, const double xed,
		const double ywd, const double yed,
		const double Fcd, const double alpha): Well(),
				xwd(xwd), xed(xed), xede(xed), ywd(ywd), yed(yed), Fcd(Fcd), alpha(alpha), boundary(boundary),
				_src_matrix(MakeSrcMatrix()){};

double Fracture::pd_lapl(const double u, const double xd, const double yd, const double zd) const {
	return 0;
};

double Fracture::pwd_lapl(const double u) const {
	return MakeMatrix(u).colPivHouseholderQr().solve(MakeRhs(u))(0);
};

double Fracture::qwd_lapl(const double u) const {
	return 1./u/u/pwd_lapl(u);
}


Eigen::MatrixXd Fracture::MakeMatrix(const double u) const {
	Eigen::MatrixXd source_matrix_ = Eigen::MatrixXd::Zero(2*NSEG+1, 2*NSEG+1);
	Eigen::MatrixXd if1_matrix_(2*NSEG, 2*NSEG),if2e_matrix_(2*NSEG, 2*NSEG), i1f2h_matrix_(2*NSEG, 2*NSEG), i2f2h_matrix_(2*NSEG, 2*NSEG);
	Eigen::VectorXd if2e_buf_(2*NSEG), i1f2h_buf_(2*NSEG);
	double mult = -1.*PI/xed;
	for (int i = 0; i < 2*NSEG; ++i) {
		source_matrix_(i,0) = 1.;
		source_matrix_(2*NSEG, i+1) = 1.;
	}
	fill_if1(u, ywd, yed, alpha, if1_matrix_);
	fill_if2e(u, xwd, xed, xede, ywd, yed, alpha, if2e_matrix_, if2e_buf_);
	fill_i1f2h(u, xwd, xed, xede, alpha, i1f2h_matrix_, i1f2h_buf_);
	fill_i2f2h(u,  ywd, alpha, i2f2h_matrix_);
	for (int i = 0; i < 2*NSEG; ++i) {
		for (int j = 0; j < 2*NSEG; ++j) {
			source_matrix_(i, j+1) = mult*(if1_matrix_(i,j)+ if2e_matrix_(i,j)+i1f2h_matrix_(i,j)+ i2f2h_matrix_(i,j))+_src_matrix(i,j);
		}
	}
	return source_matrix_;
}

Eigen::VectorXd Fracture::MakeRhs(const double u) const {
	Eigen::VectorXd rhs(2*NSEG+1);
	double coef = PI/Fcd/NSEG/u;
	for (int i = 0; i < NSEG; ++i) {
		rhs(NSEG+i) = coef*(i+0.5);
		rhs(NSEG-i-1) = rhs[NSEG+i];
	}
	rhs(2*NSEG) = 2.*NSEG/u;
	return rhs;
}

Eigen::MatrixXd Fracture::MakeSrcMatrix() const {
	Eigen::MatrixXd ans(2*NSEG, 2*NSEG);
	for (int i = 0; i < 2*NSEG; ++i) {
		for (int j = 0; j < 2*NSEG; ++j) {
			ans(i,j) = 0.;
		}
	}
	double dx = 1./NSEG;
	double dx2_8 = 0.125*dx*dx;
	double dx2_2 = 0.5*dx*dx;
	double coef = PI/Fcd;
	for (int j = 0; j < NSEG; ++j) {
		ans(j+NSEG,j+NSEG) = coef*dx2_8;
		ans(NSEG-j-1, NSEG-j-1) = ans(j+NSEG,j+NSEG);
		double xj = dx*(j+0.5);
		for (int i = 0; i < j; ++i) {
			ans(j+NSEG,i+NSEG) = coef*(dx2_2 + dx*(xj - (i+1)*dx));
			ans(NSEG-j-1, NSEG-i-1) = ans(j+NSEG,i+NSEG);
		}
	}
	return ans;
}

Eigen::VectorXd Fracture::MakeGreenVector(const double u, const double xd, const double yd, const double zd) const {
	Eigen::VectorXd ans(2*NSEG), buf(2*NSEG);
	return ans;
}

}

std::vector<double> CalcStehf(const int n) { //OK
	std::vector<double> v(n+1);
	std::vector<double> g(161);
	std::vector<double> h(81);
	g[1] = 1.;
	int NH = n/2;
	for (int i = 2; i < n+1; ++i) {
		g[i] = g[i-1]*i;
	}
	h[1] = 2./g[NH-1];
	for (int i = 2; i < NH+1; i++){
		double fi = i;
		if (i != NH) {
			h[i] = std::pow(fi, NH)*g[2*i]/(g[NH - i]*g[i]*g[i - 1]);
		} else {
			h[i] = std::pow(fi, NH)*g[2*i]/(g[i]*g[i - 1]);
		}
	}
	int SN = 2 * (NH - (NH/2)*2) - 1;
	for (int i = 1; i < n+1; ++i) {
		v[i] = 0.;
		int K1 = (i+1)/2;
		int K2 = i;
		if(K2>NH) K2=NH;
		for (int k = K1; k < K2+1; ++k) {
			if (2*k - i == 0) {
				v[i] += h[k]/g[i-k];
			} else if (i==k) {
				v[i] += h[k]/g[2*k - i];
			} else {
				v[i] += h[k]/(g[i - k]*g[2*k - i]);
			}
		}
		v[i] *= SN;
		SN = -1*SN;
	}
	return v;
}


