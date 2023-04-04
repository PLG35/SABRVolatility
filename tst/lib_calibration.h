#ifndef LIB_CALIBRATION_H
#define LIB_CALIBRATION_H

#include "lib_smile_curve.h"
#include <vector>
#include <iostream>

using namespace std;

typedef double gradient[numberOfParameters];
typedef double hessian[numberOfParameters][numberOfParameters];

class Calibrator
{
public:
    explicit Calibrator();
    explicit Calibrator(Sabr* mySabr, vector<double>* strikes, vector<double>* targetVolatilities, double targetMatchingPrecision, int limitOfIterations);
    virtual ~Calibrator();

    void calibrate();
    void fillStructure(vector<double> *strikes, int calibrationDate, int optionMaturity);
    void fillMarketData(double forward, vector<double> *targetVolatilities);

protected:
    // Market data
    Sabr* mySabr;
    vector<double>* strikes;
    vector<double>* targetVolatilities;

    // Calibration parameters
    double targetMatchingPrecision;
    int limitOfIterations;
    int numberOfPoints;

    // Used by calibration
    void increment(gradient* gradIncrement);

    // Measures
    double computeError();
    void gradientOfError(gradient* gradientError);
    double normOfGradientOfError();
    void hessianOfError(hessian* hessianError);
    void inverseHessianOfError(hessian* inverseHessian);
};

#endif // LIB_CALIBRATION_H
