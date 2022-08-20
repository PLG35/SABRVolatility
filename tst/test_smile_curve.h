#ifndef TEST_SMILE_CURVE_H
#define TEST_SMILE_CURVE_H

#include "lib_smile_curve.h"

using namespace std;

class Test_Sabr : public Sabr
{
public:
    // Constructor & destructor
    Test_Sabr();
    ~Test_Sabr();

    // Tests Z
    bool test_dzdAlpha(double money, double bumpRatio, double precision);
    bool test_dzdBeta(double money, double bumpRatio, double precision);
    bool test_dzdVovol(double money, double bumpRatio, double precision);

    // Tests XZ
    bool test_dxzdz(double money, double bumpRatio, double precision);
    bool test_dxzdAlpha(double money, double bumpRatio, double precision);
    bool test_dxzdBeta(double money, double bumpRatio, double precision);
    bool test_dxzdRho(double money, double bumpRatio, double precision);
    bool test_dxzdVovol(double money, double bumpRatio, double precision);

    // Tests ZXZ
    bool test_dzxzdAlpha(double money, double bumpRatio, double precision);
    bool test_dzxzdBeta(double money, double bumpRatio, double precision);
    bool test_dzxzdRho(double money, double bumpRatio, double precision);
    bool test_dzxzdVovol(double money, double bumpRatio, double precision);

    // Tests W
    bool test_dWdBeta(double money, double bumpRatio, double precision);
    bool test_dLfdBeta(double money, double bumpRatio, double precision);

    // Tests Rf
    bool test_dRfdAlpha(double money, double bumpRatio, double precision);
    bool test_dRfdBeta(double money, double bumpRatio, double precision);
    bool test_dRfdRho(double money, double bumpRatio, double precision);
    bool test_dRfdVovol(double money, double bumpRatio, double precision);

    // Tests volatility
    bool test_dSdAlpha(double money, double bumpRatio, double precision);
    bool test_dSdBeta(double money, double bumpRatio, double precision);
    bool test_dSdRho(double money, double bumpRatio, double precision);
    bool test_dSdVovol(double money, double bumpRatio, double precision);
    bool test_dSdForward(double money, double bumpRatio, double precision);
};

#endif // TEST_SMILE_CURVE_H
