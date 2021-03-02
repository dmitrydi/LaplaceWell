#pragma once

#include <sstream>
#include <stdexcept>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>
#include <unordered_map>
#include <limits>
#include <cmath>

static const double SMALL = std::numeric_limits<double>::min()*10.0;

template <class T>
std::ostream& operator << (std::ostream& os, const std::vector<T>& s) {
  os << "{";
  bool first = true;
  for (const auto& x : s) {
    if (!first) {
      os << ", ";
    }
    first = false;
    os << x;
  }
  return os << "}";
}

template <class T>
std::ostream& operator << (std::ostream& os, const std::set<T>& s) {
  os << "{";
  bool first = true;
  for (const auto& x : s) {
    if (!first) {
      os << ", ";
    }
    first = false;
    os << x;
  }
  return os << "}";
}

template <class K, class V>
std::ostream& operator << (std::ostream& os, const std::map<K, V>& m) {
  os << "{";
  bool first = true;
  for (const auto& kv : m) {
    if (!first) {
      os << ", ";
    }
    first = false;
    os << kv.first << ": " << kv.second;
  }
  return os << "}";
}

template <class K, class V>
std::ostream& operator << (std::ostream& os, const std::unordered_map<K, V>& m) {
	os << "{";
	  bool first = true;
	  for (const auto& kv : m) {
	    if (!first) {
	      os << ", ";
	    }
	    first = false;
	    os << kv.first << ": " << kv.second;
	  }
	  return os << "}";
}

template<class T, class U>
void AssertEqual(const T& t, const U& u, const std::string& hint = {}) {
  if (!(t == u)) {
    std::ostringstream os;
    os << "Assertion failed: " << t << " != " << u;
    if (!hint.empty()) {
       os << " hint: " << hint;
    }
    throw std::runtime_error(os.str());
  }
}

template <typename T>
bool IsClose(const T x, const T y, const T eps) {
	return (std::abs(x-y) < eps*std::abs(y) || std::abs(x-y) < std::numeric_limits<T>::min());
}


template <class T>
bool AlmostEqual(T x, T y) {
	return std::abs(x - y) <= std::numeric_limits<T>::epsilon()*std::abs(x+y) || std::abs(x-y) < std::numeric_limits<T>::min();
}

template<class T>
void AssertClose(const T& t, const T& u, const T& eps, const std::string& hint = {}) {
	if (!IsClose(t, u, eps)) {
		AssertEqual(t, u, hint);
	}
}

inline void Assert(bool b, const std::string& hint) {
  AssertEqual(b, true, hint);
}

class TestRunner {
public:
  template <class TestFunc>
  void RunTest(TestFunc func, const std::string& test_name) {
    try {
      func();
      std::cerr << test_name << " OK" << std::endl;
    } catch (std::exception& e) {
      ++fail_count;
      std::cerr << test_name << " fail: " << e.what() << std::endl;
    } catch (...) {
      ++fail_count;
      std::cerr << "Unknown exception caught" << std::endl;
    }
  }

  ~TestRunner() {
    if (fail_count > 0) {
      std::cerr << fail_count << " unit tests failed. Terminate" << std::endl;
      exit(1);
    }
  }

private:
  int fail_count = 0;
};

#define ASSERT_EQUAL(x, y) {            \
  std::ostringstream os;                     \
  os << #x << " != " << #y << ", "      \
    << __FILE__ << ":" << __LINE__;     \
  AssertEqual(x, y, os.str());          \
}

#define ASSERT_CLOSE(x, y, eps) { 		\
	std::ostringstream os;					\
	os << std::abs(x-y)/std::abs(y) << ", "	\
	<<__FILE__ << ":" << __LINE__;		\
	AssertClose(x, y, eps, os.str());				\
}

#define ASSERT(x) {                     \
  std::ostringstream os;                     \
  os << #x << " is false, "             \
    << __FILE__ << ":" << __LINE__;     \
  Assert(x, os.str());                  \
}

#define RUN_TEST(tr, func) \
  tr.RunTest(func, #func)
