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

#include "bdcpp.h"

#include <Eigen/Dense>
#include <cmath>

namespace bdcpp
{
namespace details
{
// Fit a polynomial of given order on the curve (minimizing squared distance).
std::vector<value_type> polyFit(const curve_data_type& curve, const size_t order)
{
    const size_t numCoefficients = order + 1;
    const size_t nCount = curve.size();

    Eigen::MatrixX<value_type> X(nCount, numCoefficients);
    Eigen::MatrixX<value_type> Y(nCount, 1);

    // fill X and Y matrices (X is a Vandermonde matrix)
    for (size_t row = 0; row < nCount; ++row)
    {
        Y(row, 0) = curve[row].second;
        value_type v = (value_type)1;
        for (size_t col = 0; col < numCoefficients; ++col)
        {
            X(row, col) = v;
            v *= curve[row].first;
        }
    }

    // Solve for the polynomial coefficients (one column) and return.
    const Eigen::VectorX<value_type> coefficients = X.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(Y);
    return std::vector<value_type>(coefficients.data(), coefficients.data() + numCoefficients);
}

// Calculates Y(x), where the polynomial coefficients are given.
value_type polyVal(const std::vector<value_type>& coefficients, const value_type x)
{
    assert(!coefficients.empty());
    size_t c = coefficients.size();
    value_type r = coefficients[--c];
    while (c != 0)
    {
        r *= x;
        r += coefficients[--c];
    }
    return r;
}

std::vector<value_type> polyIntegrate(const std::vector<value_type>& coefficients, const value_type constant = 0)
{
    const size_t numCoefficients = coefficients.size();
    std::vector<value_type> ic(numCoefficients + 1);
    ic[0] = constant;
    for (size_t c = 0; c < numCoefficients; ++c)
    {
        ic[c + 1] = coefficients[c] / (c + 1);
    }
    return ic;
}

// The main Bjontegaard calculation to get the area surface difference between two curves.
value_type bdDiff(const curve_data_type& curveA, const curve_data_type& curveB, const int polyOrder)
{
    assert(polyOrder >= 3);
    // Take lowest and highest X values (assumes sorted curves and relevant range overlap).
    const auto lowX = std::max(curveA.front().first, curveB.front().first);
    const auto highX = std::min(curveA.back().first, curveB.back().first);

    // Fit curves as polynomials and integrate them.
    const auto iCoefficientsA = polyIntegrate(polyFit(curveA, (size_t)polyOrder));
    const auto iCoefficientsB = polyIntegrate(polyFit(curveB, (size_t)polyOrder));

    // Calculate the definite integrals.
    const auto intA = polyVal(iCoefficientsA, highX) - polyVal(iCoefficientsA, lowX);
    const auto intB = polyVal(iCoefficientsB, highX) - polyVal(iCoefficientsB, lowX);

    // Return the BD diff (as the area over range).
    return (intB - intA) / (highX - lowX);
}

// Prepare curve points for BD calculations.
template<bool TRANSPOSE>
curve_data_type prepareCurve(const curve_data_type& curve)
{
    assert(curve.size() >= 4);
    auto newCurve(curve);
    std::sort(newCurve.begin(), newCurve.end(), [](const curve_data_point_type& vA, const curve_data_point_type& vB)
        {
            return vA.first < vB.first;
        });
    std::for_each(newCurve.begin(), newCurve.end(), [](curve_data_point_type& v)
        {
            v.first = std::log(v.first);
            if constexpr (TRANSPOSE)
            {
                std::swap(v.first, v.second);
            }
        });
    return newCurve;
}
}

// Calculate the BD-SNR for the two given curves (returns a distortion improvement in dB).
value_type bdsnr(const curve_data_type& curveA, const curve_data_type& curveB, const int polyOrder)
{
    return details::bdDiff(details::prepareCurve<false>(curveA), details::prepareCurve<false>(curveB), polyOrder);
}

// Calculate the BD-BR for the two given curves (returns a rate improvement in %).
value_type bdbr(const curve_data_type& curveA, const curve_data_type& curveB, const int polyOrder)
{
    return (std::exp(details::bdDiff(details::prepareCurve<true>(curveA), details::prepareCurve<true>(curveB), polyOrder)) - 1) * 100;
}
}
