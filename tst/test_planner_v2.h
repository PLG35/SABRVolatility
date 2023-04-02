#ifndef TEST_PLANNER_H
#define TEST_PLANNER_H

#include <vector>
#include <iostream>
#include <string>
#include "lib_smile_curve.h"

class Sabr;

class TestDerivative : public Sabr
{
public:
    explicit TestDerivative();
    explicit TestDerivative(string name, paramGetter pGetter, paramSetter pSetter, computer sComputer, computer vComputer, double money, double relativeBump, double absoluteThreshold, double relativeThreshold);
    TestDerivative(const TestDerivative &Source);
    ~TestDerivative();

    bool launch();
    string getName();

protected:
    string name;
    paramGetter pGetter;
    paramSetter pSetter;
    computer sComputer;
    computer vComputer;
    double money;
    double relativeBump;
    double absoluteThreshold;
    double relativeThreshold;
};

class TestPlanDerivatives
{
public:
    explicit TestPlanDerivatives();
    TestPlanDerivatives(const TestPlanDerivatives &Source);
    TestPlanDerivatives &operator=(const TestPlanDerivatives &source);
    virtual ~TestPlanDerivatives();

    void addTest(TestDerivative *);
    void run();

protected:
    std::vector<TestDerivative *> myTests;
    void clean();
    void addCopyOfTest(TestDerivative *);
};

class TestIntermediates : public Sabr
{
public:
    explicit TestIntermediates();
    explicit TestIntermediates(string name, computer pComputer, computerIntermediate sComputer, computerIntermediate vComputer, double money, double relativeBump, double absoluteThreshold, double relativeThreshold);
    TestIntermediates(const TestIntermediates &Source);
    ~TestIntermediates();

    bool launch();
    string getName();

protected:
    string name;
    computer pComputer;
    computerIntermediate sComputer;
    computerIntermediate vComputer;
    double money;
    double relativeBump;
    double absoluteThreshold;
    double relativeThreshold;
};

class TestPlanIntermediates
{
public:
    explicit TestPlanIntermediates();
    TestPlanIntermediates(const TestPlanIntermediates &Source);
    TestPlanIntermediates &operator=(const TestPlanIntermediates &source);
    virtual ~TestPlanIntermediates();

    void addTest(TestIntermediates *);
    void run();

protected:
    std::vector<TestIntermediates *> myTests;
    void clean();
    void addCopyOfTest(TestIntermediates *);
};

class MainTestPlan
{
public:
    explicit MainTestPlan();
    virtual ~MainTestPlan();

    void run();

protected:
    TestPlanDerivatives myTestsDerivatives;
    void fillTestsDerivatives();

    TestPlanIntermediates myTestsIntermediates;
    void fillTestsIntermediates();

    TestPlanDerivatives myTestsFallbackAlpha;
    void fillTestsFallbackAlpha();

    TestPlanDerivatives myTestsFallbackBeta;
    void fillTestsFallbackBeta();

    TestPlanDerivatives myTestsFallbackRho;
    void fillTestsFallbackRho();

    TestPlanDerivatives myTestsFallbackVovol;
    void fillTestsFallbackVovol();

    TestPlanDerivatives myTestsFallbackForward;
    void fillTestsFallbackForward();
};

#endif // TEST_PLANNER_H
