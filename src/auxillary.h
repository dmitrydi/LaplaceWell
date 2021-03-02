/*
 * auxillary.h
 *
 *  Created on: 7 дек. 2020 г.
 *      Author: Dmitry_Di
 */

#pragma once

#include <cmath>
#include <vector>

template <typename T>
std::vector<T> LogSpaced(const T min, const T max, const int steps) {
	std::vector<double> ans(steps);
	double ln_min = std::log(min);
	double ln_max = std::log(max);
	double d = (ln_max - ln_min)/(steps - 1);
	for (int i = 0; i < steps; ++i) {
		ans[i] = std::exp(ln_min + i*d);
	}
	return ans;
}

template <typename T>
std::vector<T> LinSpaced(const T min, const T max, const int steps) {
	std::vector<T> ans(steps);
	double d = (max - min)/(steps - 1);
	for (int i = 0; i < steps; ++i) {
		ans[i] = min + i*d;
	}
	return ans;
}

template <typename It>
bool SafeNext(It& it, It it_end, int step) {
	for (int i = 0; i < step; ++i) {
		if (++it == it_end) return false;
	}
	return true;
}

template <typename It> class IteratorRange {
public:
	IteratorRange(It first_, It last_) {first = first_, last = last_ ; };
	It begin() const { return first; }
	It end() const { return last; }
	size_t size() const { return distance(first, last); };
private:
	It first, last;
};

template <typename Iterator>
class Paginator {
public:
	Paginator(Iterator begin, Iterator end, size_t page_size): p_size(page_size) {
		Iterator cbegin = begin;
		while (cbegin != end) {
			auto cend = next(cbegin, std::min(p_size, static_cast<size_t>(distance(cbegin, end))));
			pages.push_back( {cbegin, cend} );
			cbegin = cend;
		}
	}
	size_t size() const {
		return pages.size();
	}
	auto begin() const { return pages.begin(); }
	auto end() const { return pages.end(); }
private:
	size_t p_size;
	std::vector<IteratorRange<Iterator>> pages;
};

template <typename C>
auto Paginate(C& c, size_t page_size) {
	return Paginator(c.begin(), c.end(), page_size);
}

template <typename C>
auto NPaginate(C& c, size_t n) {
	size_t n_el = c.size();
	size_t p_size;
	if (n_el % n == 0) {
		p_size = n_el/n;
	} else {
		p_size = n_el/n + 1;
	}
	return Paginate(c, p_size);
}


struct PointXYZV {
	PointXYZV(): x(0), y(0), z(0), val(0) {};
	PointXYZV(double x, double y, double z): x(x), y(y), z(z), val(0) {};
	double x, y, z, val;
};

typedef std::vector<PointXYZV> VectorXYZV;
typedef std::vector<VectorXYZV> MatrixXYZV;


MatrixXYZV MakeGrid(const std::vector<double>& xs, const std::vector<double>& ys, const std::vector<double>& zs = {0.});
