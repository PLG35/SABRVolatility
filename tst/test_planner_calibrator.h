#ifndef TEST_PLANNER_CALIBRATOR_H
#define TEST_PLANNER_CALIBRATOR_H

#include "lib_calibration.h"
#include <iostream>
#include <string>

class Calibrator;

class TestGradient : public Calibrator
{
public:
    explicit TestGradient();
    explicit TestGradient(Sabr* mySabr, string name, paramGetter pGetter, paramSetter pSetter, int positionInGradient, double relativeBump, double absoluteThreshold, double relativeThreshold);
    virtual ~TestGradient();

    bool launch();

protected:
    string name;
    paramGetter pGetter;
    paramSetter pSetter;
    int positionInGradient;
    double relativeBump;
    double absoluteThreshold;
    double relativeThreshold;
};

class TestHessian : public Calibrator
{
public:
    explicit TestHessian();
    explicit TestHessian(Sabr* mySabr, string name, paramGetter pGetter, paramSetter pSetter, int positionInGradient, int column, double relativeBump, double absoluteThreshold, double relativeThreshold);
    virtual ~TestHessian();

    bool launch();

protected:
    string name;
    paramGetter pGetter;
    paramSetter pSetter;
    int positionInGradient;
    int column;
    double relativeBump;
    double absoluteThreshold;
    double relativeThreshold;
};

class TestInverseHessian : public Calibrator
{
public:
    explicit TestInverseHessian();
    explicit TestInverseHessian(Sabr* mySabr, double absoluteThreshold);
    virtual ~TestInverseHessian();

    void launch();

protected:
    bool testElement(hessian *hess, int line, int column, double absoluteThreshold);
    double absoluteThreshold;
    double relativeThreshold;
};

class TestCalibrator
{
public:
    explicit TestCalibrator();
    explicit TestCalibrator(double targetMatchingPrecision, int limitOfIterations);
    virtual ~TestCalibrator();

    void fillStructure(vector<double> *strikes, int calibrationDate, int optionMaturity);
    void fillMarketData(double forward, vector<double> *targetVolatilities);
    void run();

protected:
    // Market data
    Sabr* mySabr;
    vector<double>* strikes;
    vector<double>* targetVolatilities;

    // Calibration parameters
    double targetMatchingPrecision;
    int limitOfIterations;
    int numberOfPoints;

    // Tests
    std::vector<TestGradient *> myTestsGradient;
    void fillTestsGradient();
    std::vector<TestHessian *> myTestsHessian;
    void fillTestsHessian();
    double absoluteThresholdForInverseHessian = 0.000001;
    TestInverseHessian* myTestInverseHessian;
};

#endif // TEST_PLANNER_CALIBRATOR_H
