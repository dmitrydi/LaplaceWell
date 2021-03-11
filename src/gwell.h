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
#include <future>
#include <iterator>
#include <exception>
#include <cassert>
#include "chbessel.h"
#include "profile.h"
#include "quadrature.h"
#include "auxillary.h"


static const int NCOEF = 12;

std::vector<double> CalcStehf(const int n);



class LaplWell {
public:
	LaplWell();
	virtual double pd_lapl(const double u, const double xd, const double yd, const double zd = 0.) const = 0;
	virtual double pwd_lapl(const double u) const = 0;
	virtual double qwd_lapl(const double u) const = 0;
	virtual ~LaplWell();

	double pwd(const double td) const;
	double qwd(const double td) const;
	double pd(const double td, const double xd, const double yd, const double zd = 0.) const;
	void pwd_parallel(const std::vector<double>& tds, std::vector<double>& pwds, int nthreads) const;
	void qwd_parallel(const std::vector<double>& tds, std::vector<double>& qwds, int nthreads) const;

	Matrix3DV pd_m_parallel(const double td, int nthreads, const std::vector<double>& xs,
			const std::vector<double>& ys,
			const std::vector<double>& zs = {0.}) const;
	Matrix3DV pd_lapl_m(const double u, const Matrix3DV& grid, int nthread) const;
	void pd_lapl_m(const double u, const Matrix3DV& grid, Matrix3DV& buf, int nthread) const;


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

	template <typename Func>
	void InverseLaplaceParallel(Func func, const std::vector<double>& tds, std::vector<double>& props, int nthreads) const {
		assert (tds.size() == props.size());
		std::vector<std::future<void>> fut;
		auto xbegin = tds.begin();
		auto ybegin = props.begin();
		for (int i = 0; i < nthreads; ++i) {
			fut.push_back(async(std::launch::async, [&]{return TSingleThread(func, xbegin++, tds.end(), ybegin++, nthreads);}));
		}
	}

	template <typename Func>
	double InverseLaplaceXYZ(Func func, const double td, const double x, const double y, const double z = 0.) const {
		double s_mult = std::log(2.)/td;
		double ans = 0.;
		for (int i = 1; i <= NCOEF; ++i) {
			double s = i*s_mult;
			{
			ans += (this->*func)(s, x, y, z)*s*stehf_coefs[i]/i;
			}
		}
		return ans;
	}
	template <typename Func, typename Iterator1, typename Iterator2>
	void TSingleThread(Func func, Iterator1 x_begin, Iterator1 x_end, Iterator2 y_begin, int step) const {
		do {
			*y_begin = (this->*func)(*x_begin);
			y_begin = next(y_begin, step);
		} while(SafeNext(x_begin, x_end, step));
	}

protected:
	const std::vector<double> stehf_coefs;
};

namespace Rectangular {

static const int NSEG = 20;
static const double SUM_EPS = 1e-8;
static const double INT_EPS = 1e-12;
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
			Eigen::VectorXd& buf) const;

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

	void vect_i1f2h_yd_verbose(const double u,
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
