# LaplaceWell

LaplaceWell is a C++ library providing classes for pressure distribution and productivity calculation of oil wells located in various drainage areas with various boundary conditions.
All parameters (pressure and productivity) calculations are based on unsteady-state solution of piezoconductivity equation with sources in all time domain. Thus no artifical splices are used for parameters calculation during transient states.

## Compile

Uses [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page) library and ```C ``` quadmath library

Requirements:
1. Eigen library installed
2. C++17 or higher compiler

Required compile options:
```bash
-fext-numeric-literals
```
Required linker options:
```
-lquadmath
```

## Usage

There are several namespaces denoting type of drainage area:
1. Rectangular
2. Circle (not implemented)
3. Infinite (not implemented)

For every namespace several classes denoting a type of a well provided:
1. Fracture
2. Multifractured (not implemented)
2. Horizontal (not implemented)
3. Vertical (not implemented)

Thus, for example, instantiation of a fractured well located in a rectangular drainage area will look like this:
```C++
#include "gwell.h"

int main() {
    // parameters definition
    Rectangular::Fracture frac(boundary,
        xwd, xed,
        ywd, yed,
        Fcd, alpha);
    // usage of frac
}
```

For every class a bunch of methods provided.

### For calculations in Laplace space:
1. ```pd_lapl(const double u, const double xd, const double yd, const double zd = 0.) const``` - returns dimentionless pressure at space point ```(xd, yd, zd)``` in Laplace space at Laplace' parameter ```u```.
2. ```double pwd_lapl(const double u) const``` - returns dimentionless wellbore pressure in Laplace space at Laplace' parameter ```u```.
3. ```double qwd_lapl(const double u) const```- returns dimentionless wellbore flow in Laplace space at Laplace' parameter ```u```.

### For calculations in real space:
1. ``` double pd(const double td, const double xd, const double yd, const double zd = 0.) const``` - returns dimentionless pressure at space point ```(xd, yd, zd)``` at time ```td```.
2. ```double pwd(const double td) const``` - returns dimentionless wellbore pressure at time ```td```.
3. ```double qwd(const double td) const``` - returns dimentionless wellbore flow at time ```td```.
4. ```void pwd_parallel(const std::vector<double>& tds, std::vector<double>& pwds, int nthreads) const```;
5. ```void qwd_parallel(const std::vector<double>& tds, std::vector<double>& qwds, int nthreads) const```;
6. ```Matrix3DV pd_m_parallel(const double td, int nthreads, const std::vector<double>& xs, const std::vector<double>& ys, const std::vector<double>& zs = {0.}) const``` - returns an ```Matrix3DV``` object containing grid ```xs x ys x zs``` and pressure values at time ```td```.

## License
[MIT](https://choosealicense.com/licenses/mit/)