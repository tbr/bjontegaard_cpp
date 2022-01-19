/*++
    BDcpp -- Simple Bjontegaard Delta metric implementation for C++.

    MIT License

    Copyright (c) 2022 Tim Bruylants

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.
--*/

#pragma once

#include <vector>

namespace bdcpp
{
typedef double value_type;  // base value type
typedef std::pair<value_type, value_type> curve_data_point_type;  // one (x, y) pair as (rate, psnr)
typedef std::vector<curve_data_point_type> curve_data_type;  // (X, Y) values of an RD curve (X=rate, Y=psnr)

value_type bdsnr(const curve_data_type& curveA, const curve_data_type& curveB, const int polyOrder = 3);
value_type bdbr(const curve_data_type& curveA, const curve_data_type& curveB, const int polyOrder = 3);
}