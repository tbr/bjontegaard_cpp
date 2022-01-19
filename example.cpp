/*++
    BDcpp -- Simple Bjontegaard Delta metric implementation for C++.

    Example code to demonstrate usage.

    Copyright (c) Tim bruylants. All rights reserved.
--*/

#include "bdcpp.h"

#include <cassert>
#include <iostream>

bdcpp::value_type round(const bdcpp::value_type& v, const int digits)
{
    const double p = std::pow(10., digits);    
    return (bdcpp::value_type)((int64_t)(v * p + (v > 0 ? 0.5 : -0.5)) / p);
}

int main(int argc, char* argv[])
{
    (void)argc;
    (void)argv;

    // VCEG-M33 numbers.
    bdcpp::curve_data_type m33_curveA = {
        {90.33, 27.95},
        {181.03, 31.17},
        {332.99, 34.44},
        {547.79, 38.11}
    };
    bdcpp::curve_data_type m33_curveB = {
        {127.1719, 29.9286},
        {186.14, 32.4165},
        {307.2924, 34.7013},
        {481.5588, 37.4561}
    };

    // ETRO barb512 numbers.
    bdcpp::curve_data_type etro_curveA = {
        {2.99899, 48.6681},
        {1.99884, 43.8357},
        {1.49673, 41.3982},
        {0.99707, 38.0124},
        {0.745941, 35.5122},
        {0.596436, 33.8673},
        {0.495148, 32.6658},
        {0.29425, 29.4505},
        {0.244049, 28.484},
    };
    bdcpp::curve_data_type etro_curveB = {
        {2.997253, 48.5637},
        {1.968292, 43.7318},
        {1.472565, 41.2727},
        {0.994965, 38.0001},
        {0.748871, 35.6375},
        {0.58902, 33.9057},
        {0.499176, 32.8318},
        {0.296753, 29.6851},
        {0.249634, 28.8244},
    };

    {
        const auto bdsnr = bdcpp::bdsnr(m33_curveA, m33_curveB);
        const auto bdbr = bdcpp::bdbr(m33_curveA, m33_curveB);
        std::cout << bdsnr << " -- " << bdbr << "\n";
        assert(round(bdsnr, 6) == 0.800713);
        assert(round(bdbr, 6) == -13.261773);
    }

    {
        const auto bdsnr = bdcpp::bdsnr(etro_curveA, etro_curveB);
        const auto bdbr = bdcpp::bdbr(etro_curveA, etro_curveB);
        std::cout << bdsnr << " -- " << bdbr << "\n";
        assert(round(bdsnr, 6) == 0.072297);
        assert(round(bdbr, 6) == -0.903450);
    }

	return 0;
}