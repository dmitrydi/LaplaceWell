/*
 * gwell.h
 *
 *  Created on: 25 февр. 2021 г.
 *      Author: Dmitry_Di
 */

#pragma once

#include <Eigen/Dense>
#include <cmath>
#include <limits>
#include <exception>
#include <iostream>
#include <sstream>
#include <iomanip>
#include "chbessel.h"
#include "profile.h"
#include "quadrature.h"
#include <mutex>

static const int NCOEF = 12;

std::vector<double> CalcStehf(const int n);

class LaplWell {
public:
	LaplWell();
	double pwd(const double td) const;
	double qwd(const double td) const;
	double pd(const double u, const double xd, const double yd, const double zd = 0.) const;
	template <typename Func>
	double InverseLaplace(Func func, const double td) const {
		double s_mult = std::log(2.)/td;
		double ans = 0.;
		for (int i = 1; i <= NCOEF; ++i) {
			double s = i*s_mult;
			{
			ans += (this->*func)(s)*s*stehf_coefs[i]/i;
			}
		}
		return ans;
	}
	virtual double pd_lapl(const double u, const double xd, const double yd, const double zd = 0.) const = 0;
	virtual double pwd_lapl(const double u) const = 0;
	virtual double qwd_lapl(const double u) const = 0;
	virtual ~LaplWell();
protected:
	const std::vector<double> stehf_coefs;
};

namespace Rectangular {

static const int NSEG = 20;
static const double SUM_EPS = 1e-8;
static const double INT_EPS = 1e-9;
static const int KMAX = 10000;
static const int KMIN = 10;
static const double PI = 3.141592653589793;
static const double TINY = std::numeric_limits<double>::min();


enum class Boundary {
	NNNN,
	CCCC
};

class Well: public LaplWell {
public:
	Well();
	virtual ~Well();

protected:
	const FastBessel::Bess bess;
	const double dx;
	virtual Eigen::MatrixXd MakeMatrix(const double u) const = 0;
	virtual Eigen::VectorXd MakeRhs(const double u) const = 0;
	virtual Eigen::MatrixXd MakeSrcMatrix() const = 0;
	virtual Eigen::VectorXd MakeGreenVector(const double u, const double xd, const double yd, const double zd = 0.) const = 0;
	double SEXP(const double y, const double e) const;
	void fill_if1(const double u,
			const double ywd, const double yed,
			const double alpha,
			Eigen::MatrixXd& matrix) const;
	void vect_if1_yd(const double u,
			const double yd, const double ywd, const double yed,
			const double alpha,
			Eigen::VectorXd& buf);

	void fill_if2e(const double u,
			const double xwd,
			const double xed, const double xede,
			const double ywd, const double yed,
			const double alpha,
			Eigen::MatrixXd& matrix, Eigen::VectorXd& buf) const; // OK
	void vect_if2e_yd(const double u, const double xd, const double xwd,
			const double xed, const double xede,
			const double yd, const double ywd, const double yed,
			const double alpha, Eigen::VectorXd& buf) const;

	void fill_i1f2h(const double u,
			const double xwd, const double xed, const double xede, const double alpha,
			Eigen::MatrixXd& matrix, Eigen::VectorXd& buf) const;
	void vect_i1f2h(const double u,
			const double xd, const double xwd, const double xed, const double xede,
			const double alpha,
			Eigen::VectorXd& buf) const;
	void vect_i1f2h_yd(const double u,
			const double xd, const double xwd, const double xed, const double xede,
			const double yd, const double ywd,
			const double alpha, Eigen::VectorXd& buf) const;

	void fill_i2f2h(const double u, const double ywd, const double alpha, Eigen::MatrixXd& matrix) const;
	void vect_i2f2h_yd(const double u, const double yd, const double ywd, const double alpha, Eigen::VectorXd& buf) const;

};

class Fracture: public Well {
public:
	Fracture(const Boundary boundary, const double xwd, const double xed,
			const double ywd, const double yed,
			const double Fcd, const double alpha = 0.);
	double pd_lapl(const double u, const double xd, const double yd, const double zd = 0.) const override;
	double pwd_lapl(const double u) const override;
	double qwd_lapl(const double u) const override;
private:
	const double xwd, xed, xede, ywd, yed, Fcd, alpha;
	const Boundary boundary;
	const Eigen::MatrixXd _src_matrix;
	Eigen::MatrixXd MakeMatrix(const double u) const override;
	Eigen::VectorXd MakeRhs(const double u) const override;
	Eigen::MatrixXd MakeSrcMatrix() const override;
	Eigen::VectorXd MakeGreenVector(const double u, const double xd, const double yd, const double zd = 0.) const override;

};

}
