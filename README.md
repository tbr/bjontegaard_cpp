# Bjontegaard Metric implementation for C++17 (or later).

## Description:

This project contains an implementation of the Bjontegaard metric to calculate Delta SNR and Delta Rate, using arbitrary number of data points. The calculated results comply with VCEG-M33 when using 4 data points. This code relies on the Eigen3 library to perform the polynomial fitting of curves.

The library provides two functions (see `bdcpp.h` for type details):

```cpp
value_type bdsnr(const curve_data_type& curveA, const curve_data_type& curveB, const int polyOrder = 3);

value_type bdbr(const curve_data_type& curveA, const curve_data_type& curveB, const int polyOrder = 3);
```

Normally, curve fitting is done using a 3rd order polynomial, however this code can also do higher order approximations to improve precision.

Note: This code is a port of the original VBA implementation (see [Bjontegaard ETRO project](https://github.com/tbr/bjontegaard_etro)) to C++.

## Building:

No special steps are required.

The example build uses Conan and CMake. Create a build folder and from there issue:

```
conan install <path_to_source_root> -s build_type=Release
conan install <path_to_source_root> -s build_type=Debug
cmake <path_to_source_root>
```

## Author:

Tim Bruylants

## License:

MIT license, see included LICENSE.txt file.

## References:

	[1] G. Bjontegaard, Calculation of average PSNR differences between RD-curves (VCEG-M33)

_Copyright (C) 2022 Tim Bruylants._
