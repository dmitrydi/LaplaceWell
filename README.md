# LaplaceWell

LaplaceWell is a C++ classes for pressure distribution and productivity calculation of a finite-conductivity fractured oil well located in a closed rectangular drainage area.
All parameters (pressure and productivity) calculations are based on unsteady-state solution of piezoconductivity equation with sources in all time domain. Thus no artifical splices are used for calculation during transient states.

## Compile

Uses [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page) library and C-quadmath library

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

Instantiation of a fractured well located in a rectangular drainage area will look like this:
```C++
#include "gwell.h"

int main() {
    // parameters definition
	...
    //
    Rectangular::Fracture frac(Boundary::NNNN,
        xwd, xed,
        ywd, yed,
        Fcd, alpha);
    // usage of frac
    ...
    //
}
```

### For calculations in Laplace space:
1. ```double pd_lapl(const double u, const double xd, const double yd, const double zd = 0.) const``` - returns dimentionless pressure at space point ```(xd, yd, zd)``` in Laplace space at Laplace' parameter ```u```.
2. ```double pwd_lapl(const double u) const``` - returns dimentionless wellbore pressure in Laplace space at Laplace' parameter ```u```.
3. ```double qwd_lapl(const double u) const```- returns dimentionless wellbore flow in Laplace space at Laplace' parameter ```u```.

### For calculations in real space:
1. ``` double pd(const double td, const double xd, const double yd, const double zd = 0.) const``` - returns dimentionless pressure at space point ```(xd, yd, zd)``` at time ```td```.
2. ```double pwd(const double td) const``` - returns dimentionless wellbore pressure at time ```td```.
3. ```double qwd(const double td) const``` - returns dimentionless wellbore flow at time ```td```.
4. ```void pwd_parallel(const std::vector<double>& tds, std::vector<double>& pwds, int nthreads) const``` - same as (2) done in parallel with ```nthreads``` threads;
5. ```void qwd_parallel(const std::vector<double>& tds, std::vector<double>& qwds, int nthreads) const``` - same as (3) done in parallel with ```nthreads``` threads;
6. ```Matrix3DV pd_m_parallel(const double td, int nthreads, const std::vector<double>& xs, const std::vector<double>& ys, const std::vector<double>& zs = {0.}) const``` - returns an ```Matrix3DV``` object containing grid ```xs x ys x zs``` and pressure values at time ```td```. Done in parallel with ```nthreads``` threads.

## License
[MIT](https://choosealicense.com/licenses/mit/)