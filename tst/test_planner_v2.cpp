#include "test_planner_v2.h"

////////////////////
// TestDerivative //
////////////////////

TestDerivative::TestDerivative(){};

TestDerivative::TestDerivative(string n, paramGetter pG, paramSetter pS, computer sC, computer vC, double m, double rB, double aT, double rT){
    this->name = n;
    this->pGetter = pG;
    this->pSetter = pS;
    this->sComputer = sC;
    this->vComputer = vC;
    this->money = m;
    this->relativeBump = rB;
    this->absoluteThreshold = aT;
    this->relativeThreshold = rT;
}

TestDerivative::TestDerivative(const TestDerivative &Source){
    this->name = Source.name;
    this->pGetter = Source.pGetter;
    this->pSetter = Source.pSetter;
    this->sComputer = Source.sComputer;
    this->vComputer = Source.vComputer;
    this->money = Source.money;
    this->relativeBump = Source.relativeBump;
    this->absoluteThreshold = Source.absoluteThreshold;
    this->relativeThreshold = Source.relativeThreshold;
}

TestDerivative::~TestDerivative(){}

bool TestDerivative::launch(){
    double strike = this->getForward() + this->money;
    double sensitivity = (this->*(this->sComputer))(strike);

    double z1 = (this->*vComputer)(strike);

    double parameter = (this->*(this->pGetter))();
    double dparameter = parameter * this->relativeBump;
    (this->*(this->pSetter))(parameter + dparameter);
    double z2 = (this->*vComputer)(strike);

    double variation = (z2 - z1) / dparameter;

    // Restore parameter
    (this->*(this->pSetter))(parameter);
    double relativeError = std::abs(sensitivity / variation - 1);
    double absoluteError = std::abs(sensitivity - variation);

    bool testResult = relativeError < this->relativeThreshold && absoluteError < this->absoluteThreshold;
    if (testResult || (abs(variation) <= precisionMachine && abs(sensitivity) <= precisionMachine)){
        cout << "Test " + this->name + " passed." << endl;
    }
    else{
        cout << "Test " + this->name + " failed." << endl;
        cout << "   Sensitivity : " << sensitivity << "." << endl;
        cout << "   Variation reached : " << variation << "." << endl;
        cout << "   Relative bump by " << (this->relativeBump * 100) << "%." << endl;
        cout << "   Relative error expected : " << (this->relativeThreshold * 100) << "%." << endl;
        cout << "   Relative error reached : " << (relativeError * 100) << "%." << endl;
        cout << "   Absolute error expected : " << this->absoluteThreshold << "." << endl;
        cout << "   Absolute error reached : " << absoluteError << "." << endl;
    }
    return testResult;
}

string TestDerivative::getName(){
    return this->name;
}

/////////////////////////
// TestPlanDerivatives //
/////////////////////////

TestPlanDerivatives::TestPlanDerivatives(){}

void TestPlanDerivatives::addCopyOfTest(TestDerivative * test){
    TestDerivative *myTest = new TestDerivative(*test);
    this->myTests.push_back(myTest);
}

TestPlanDerivatives::TestPlanDerivatives(const TestPlanDerivatives &Source){
    for(auto it = Source.myTests.begin(); it != Source.myTests.end(); it++){
        this->addCopyOfTest(*it);
    }
}

TestPlanDerivatives &TestPlanDerivatives::operator=(const TestPlanDerivatives &Source)
{
    if (&Source != this)
    {
        for(auto it = Source.myTests.begin(); it != Source.myTests.end(); it++){
            this->addCopyOfTest(*it);
        }
    }
    return *this;
}

TestPlanDerivatives::~TestPlanDerivatives(){
    this->clean();
}

void TestPlanDerivatives::addTest(TestDerivative * test){
    this->myTests.push_back(test);
}

void TestPlanDerivatives::run(){
    cout << "START TEST PLAN DERIVATIVES" << endl;
    for(std::vector<TestDerivative*>::iterator it = this->myTests.begin(); it != this->myTests.end(); ++it){
        (*it)->launch();
    }
    cout << "END TEST PLAN DERIVATIVES" << endl;
}

void TestPlanDerivatives::clean(){
    for(std::vector<TestDerivative*>::iterator it = this->myTests.begin(); it != this->myTests.end(); ++it){
        delete *it;
    }
}

///////////////////////
// TestIntermediates //
///////////////////////

TestIntermediates::TestIntermediates(){};

TestIntermediates::TestIntermediates(string n, computer pC, computerIntermediate sCI, computerIntermediate vCI, double m, double rB, double aT, double rT){
    this->name = n;
    this->pComputer = pC;
    this->sComputer = sCI;
    this->vComputer = vCI;
    this->money = m;
    this->relativeBump = rB;
    this->absoluteThreshold = aT;
    this->relativeThreshold = rT;
}

TestIntermediates::TestIntermediates(const TestIntermediates &Source){
    this->name = Source.name;
    this->pComputer = Source.pComputer;
    this->sComputer = Source.sComputer;
    this->vComputer = Source.vComputer;
    this->money = Source.money;
    this->relativeBump = Source.relativeBump;
    this->absoluteThreshold = Source.absoluteThreshold;
    this->relativeThreshold = Source.relativeThreshold;
}

TestIntermediates::~TestIntermediates(){}

bool TestIntermediates::launch(){
    double strike = this->getForward() + this->money;

    double parameter = (this->*pComputer)(strike);
    double sensitivity = (this->*(this->sComputer))(strike, parameter);
    double value1 = (this->*vComputer)(strike, parameter);

    double dparameter = parameter * this->relativeBump;
    double value2 = (this->*vComputer)(strike, parameter + dparameter);

    double variation = (value2 - value1) / dparameter;

    double relativeError = std::abs(sensitivity / variation - 1);
    double absoluteError = std::abs(sensitivity - variation);

    bool testResult = relativeError < this->relativeThreshold && absoluteError < this->absoluteThreshold;
    if (testResult || (abs(variation) <= precisionMachine && abs(sensitivity) <= precisionMachine) ){
        cout << "Test " + this->name + " passed." << endl;
    }
    else{
        cout << "Test " + this->name + " failed." << endl;
        cout << "   Sensitivity : " << sensitivity << "." << endl;
        cout << "   Variation reached : " << variation << "." << endl;
        cout << "   Relative bump by " << (this->relativeBump * 100) << "%." << endl;
        cout << "   Parameter bumped from " << parameter << " to " << (parameter + dparameter) << "." << endl;
        cout << "   Relative error expected : " << (this->relativeThreshold * 100) << "%." << endl;
        cout << "   Relative error reached : " << (relativeError * 100) << "%." << endl;
        cout << "   Absolute error expected : " << this->absoluteThreshold << "." << endl;
        cout << "   Absolute error reached : " << absoluteError << "." << endl;
    }
    return testResult;
}

string TestIntermediates::getName(){
    return this->name;
}

///////////////////////////
// TestPlanIntermediates //
///////////////////////////

TestPlanIntermediates::TestPlanIntermediates(){}

void TestPlanIntermediates::addCopyOfTest(TestIntermediates * test){
    TestIntermediates *myTest = new TestIntermediates(*test);
    this->myTests.push_back(myTest);
}

TestPlanIntermediates::TestPlanIntermediates(const TestPlanIntermediates &Source){
    for(auto it = Source.myTests.begin(); it != Source.myTests.end(); it++){
        this->addCopyOfTest(*it);
    }
}

TestPlanIntermediates &TestPlanIntermediates::operator=(const TestPlanIntermediates &Source)
{
    if (&Source != this)
    {
        for(auto it = Source.myTests.begin(); it != Source.myTests.end(); it++){
            this->addCopyOfTest(*it);
        }
    }
    return *this;
}

TestPlanIntermediates::~TestPlanIntermediates(){
    this->clean();
}

void TestPlanIntermediates::addTest(TestIntermediates * test){
    this->myTests.push_back(test);
}

void TestPlanIntermediates::run(){
    cout << "START TEST PLAN INTERMEDIATES" << endl;
    for(std::vector<TestIntermediates*>::iterator it = this->myTests.begin(); it != this->myTests.end(); ++it){
        (*it)->launch();
    }
    cout << "START TEST PLAN INTERMEDIATES" << endl;
}

void TestPlanIntermediates::clean(){
    for(std::vector<TestIntermediates*>::iterator it = this->myTests.begin(); it != this->myTests.end(); ++it){
        delete *it;
    }
}

//////////////////
// MainTestPlan //
//////////////////

MainTestPlan::MainTestPlan(){
    this->fillTestsDerivatives();
    this->fillTestsIntermediates();

    this->fillTestsFallbackAlpha();
    this->fillTestsFallbackBeta();
    this->fillTestsFallbackRho();
    this->fillTestsFallbackVovol();
    this->fillTestsFallbackForward();
}

MainTestPlan::~MainTestPlan(){}

void MainTestPlan::fillTestsDerivatives(){
    computer volatility = &Sabr::volatility;

    computer dSdAlpha = &Sabr::dSdAlpha;
    paramGetter getAlpha = &Sabr::getAlpha;
    paramSetter setAlpha = &Sabr::setAlpha;

    computer dSdBeta = &Sabr::dSdBeta;
    paramGetter getBeta = &Sabr::getBeta;
    paramSetter setBeta = &Sabr::setBeta;

    computer dSdRho = &Sabr::dSdRho;
    paramGetter getRho = &Sabr::getRho;
    paramSetter setRho = &Sabr::setRho;

    computer dSdVovol = &Sabr::dSdVovol;
    paramGetter getVovol = &Sabr::getVovol;
    paramSetter setVovol = &Sabr::setVovol;

    computer dSdForward = &Sabr::dSdForward;
    paramGetter getForward = &Sabr::getForward;
    paramSetter setForward = &Sabr::setForward;

    TestDerivative * test1 = new TestDerivative("dSdAlpha", getAlpha, setAlpha, dSdAlpha, volatility, 0.001, 0.00001, 0.00001, 0.0001);
    this->myTestsDerivatives.addTest(test1);
    TestDerivative * test2 = new TestDerivative("dSdBeta", getBeta, setBeta, dSdBeta, volatility, 0.001, 0.00001, 0.00001, 0.0001);
    this->myTestsDerivatives.addTest(test2);
    TestDerivative * test3 = new TestDerivative("dSdRho", getRho, setRho, dSdRho, volatility, 0.001, 0.00001, 0.00001, 0.0001);
    this->myTestsDerivatives.addTest(test3);
    TestDerivative * test4 = new TestDerivative("dSdVovol", getVovol, setVovol, dSdVovol, volatility, 0.001, 0.00001, 0.00001, 0.0001);
    this->myTestsDerivatives.addTest(test4);
    TestDerivative * test5 = new TestDerivative("dSdForward", getForward, setForward, dSdForward, volatility, 0.001, 0.00001, 0.0001, 0.0001);
    this->myTestsDerivatives.addTest(test5);

    // Second order Alpha
    computer d2SdAlphadAlpha = &Sabr::d2SdAlphadAlpha;
    computer d2SdBetadAlpha = &Sabr::d2SdBetadAlpha;
    computer d2SdRhodAlpha = &Sabr::d2SdRhodAlpha;
    computer d2SdVovoldAlpha = &Sabr::d2SdVovoldAlpha;
    computer d2SdForwarddAlpha = &Sabr::d2SdForwarddAlpha;

    TestDerivative * test6 = new TestDerivative("d2SdAlphadAlpha", getAlpha, setAlpha, d2SdAlphadAlpha, dSdAlpha, 0.001, 0.000001, 0.0001, 0.000001);
    this->myTestsDerivatives.addTest(test6);
    TestDerivative * test7 = new TestDerivative("d2SdBetadAlpha", getAlpha, setAlpha, d2SdBetadAlpha, dSdBeta, 0.001, 0.000001, 0.0001, 0.000001);
    this->myTestsDerivatives.addTest(test7);
    TestDerivative * test8 = new TestDerivative("d2SdRhodAlpha", getAlpha, setAlpha, d2SdRhodAlpha, dSdRho, 0.001, 0.000001, 0.0001, 0.000001);
    this->myTestsDerivatives.addTest(test8);
    TestDerivative * test9 = new TestDerivative("d2SdVovoldAlpha", getAlpha, setAlpha, d2SdVovoldAlpha, dSdVovol, 0.001, 0.000001, 0.0001, 0.000001);
    this->myTestsDerivatives.addTest(test9);
    TestDerivative * test10 = new TestDerivative("d2SdForwarddAlpha", getAlpha, setAlpha, d2SdForwarddAlpha, dSdForward, 0.001, 0.000001, 0.0001, 0.000001);
    this->myTestsDerivatives.addTest(test10);

    // Second order Beta
    computer d2SdAlphadBeta = &Sabr::d2SdAlphadBeta;
    computer d2SdBetadBeta = &Sabr::d2SdBetadBeta;
    computer d2SdRhodBeta = &Sabr::d2SdRhodBeta;
    computer d2SdVovoldBeta = &Sabr::d2SdVovoldBeta;
    computer d2SdForwarddBeta = &Sabr::d2SdForwarddBeta;
    TestDerivative * test11 = new TestDerivative("d2SdAlphadBeta", getBeta, setBeta, d2SdAlphadBeta, dSdAlpha, 0.001, 0.000001, 0.0001, 0.000001);
    this->myTestsDerivatives.addTest(test11);
    TestDerivative * test12 = new TestDerivative("d2SdBetadBeta", getBeta, setBeta, d2SdBetadBeta, dSdBeta, 0.001, 0.000001, 0.0001, 0.000001);
    this->myTestsDerivatives.addTest(test12);
    TestDerivative * test13 = new TestDerivative("d2SdRhodBeta", getBeta, setBeta, d2SdRhodBeta, dSdRho, 0.001, 0.000001, 0.0001, 0.000001);
    this->myTestsDerivatives.addTest(test13);
    TestDerivative * test14 = new TestDerivative("d2SdVovoldBeta", getBeta, setBeta, d2SdVovoldBeta, dSdVovol, 0.001, 0.000001, 0.0001, 0.000001);
    this->myTestsDerivatives.addTest(test14);
    TestDerivative * test15 = new TestDerivative("d2SdForwarddBeta", getBeta, setBeta, d2SdForwarddBeta, dSdForward, 0.001, 0.000001, 0.0001, 0.000001);
    this->myTestsDerivatives.addTest(test15);

    // Second order Rho
    computer d2SdAlphadRho = &Sabr::d2SdAlphadRho;
    computer d2SdBetadRho = &Sabr::d2SdBetadRho;
    computer d2SdRhodRho = &Sabr::d2SdRhodRho;
    computer d2SdVovoldRho = &Sabr::d2SdVovoldRho;
    computer d2SdForwarddRho = &Sabr::d2SdForwarddRho;
    TestDerivative * test16 = new TestDerivative("d2SdAlphadRho", getRho, setRho, d2SdAlphadRho, dSdAlpha, 0.001, 0.000001, 0.0001, 0.000001);
    this->myTestsDerivatives.addTest(test16);
    TestDerivative * test17 = new TestDerivative("d2SdBetadRho", getRho, setRho, d2SdBetadRho, dSdBeta, 0.001, 0.000001, 0.0001, 0.000001);
    this->myTestsDerivatives.addTest(test17);
    TestDerivative * test18 = new TestDerivative("d2SdRhodRho", getRho, setRho, d2SdRhodRho, dSdRho, 0.001, 0.000001, 0.0001, 0.000001);
    this->myTestsDerivatives.addTest(test18);
    TestDerivative * test19 = new TestDerivative("d2SdVovoldRho", getRho, setRho, d2SdVovoldRho, dSdVovol, 0.001, 0.000001, 0.0001, 0.000001);
    this->myTestsDerivatives.addTest(test19);
    TestDerivative * test20 = new TestDerivative("d2SdForwarddRho", getRho, setRho, d2SdForwarddRho, dSdForward, 0.001, 0.000001, 0.0001, 0.000001);
    this->myTestsDerivatives.addTest(test20);

    // Second order Vovol
    computer d2SdAlphadVovol = &Sabr::d2SdAlphadVovol;
    computer d2SdBetadVovol = &Sabr::d2SdBetadVovol;
    computer d2SdRhodVovol = &Sabr::d2SdRhodVovol;
    computer d2SdVovoldVovol = &Sabr::d2SdVovoldVovol;
    computer d2SdForwarddVovol = &Sabr::d2SdForwarddVovol;
    TestDerivative * test21 = new TestDerivative("d2SdAlphadVovol", getVovol, setVovol, d2SdAlphadVovol, dSdAlpha, 0.001, 0.000001, 0.0001, 0.000001);
    this->myTestsDerivatives.addTest(test21);
    TestDerivative * test22 = new TestDerivative("d2SdBetadVovol", getVovol, setVovol, d2SdBetadVovol, dSdBeta, 0.001, 0.000001, 0.0001, 0.000001);
    this->myTestsDerivatives.addTest(test22);
    TestDerivative * test23 = new TestDerivative("d2SdRhodVovol", getVovol, setVovol, d2SdRhodVovol, dSdRho, 0.001, 0.000001, 0.0001, 0.000001);
    this->myTestsDerivatives.addTest(test23);
    TestDerivative * test24 = new TestDerivative("d2SdVovoldVovol", getVovol, setVovol, d2SdVovoldVovol, dSdVovol, 0.001, 0.000001, 0.0001, 0.000001);
    this->myTestsDerivatives.addTest(test24);
    TestDerivative * test25 = new TestDerivative("d2SdForwarddVovol", getVovol, setVovol, d2SdForwarddVovol, dSdForward, 0.001, 0.000001, 0.0001, 0.000001);
    this->myTestsDerivatives.addTest(test25);

    // Second order Forward
    computer d2SdAlphadForward = &Sabr::d2SdAlphadForward;
    computer d2SdBetadForward = &Sabr::d2SdBetadForward;
    computer d2SdRhodForward = &Sabr::d2SdRhodForward;
    computer d2SdVovoldForward = &Sabr::d2SdVovoldForward;
    computer d2SdForwarddForward = &Sabr::d2SdForwarddForward;
    TestDerivative * test26 = new TestDerivative("d2SdAlphadForward", getForward, setForward, d2SdAlphadForward, dSdAlpha, 0.001, 0.000001, 0.0001, 0.000001);
    this->myTestsDerivatives.addTest(test26);
    TestDerivative * test27 = new TestDerivative("d2SdBetadForward", getForward, setForward, d2SdBetadForward, dSdBeta, 0.001, 0.000001, 0.0001, 0.000001);
    this->myTestsDerivatives.addTest(test27);
    TestDerivative * test28 = new TestDerivative("d2SdRhodForward", getForward, setForward, d2SdRhodForward, dSdRho, 0.001, 0.000001, 0.0001, 0.000001);
    this->myTestsDerivatives.addTest(test28);
    TestDerivative * test29 = new TestDerivative("d2SdVovoldForward", getForward, setForward, d2SdVovoldForward, dSdVovol, 0.001, 0.000001, 0.0001, 0.000001);
    this->myTestsDerivatives.addTest(test29);
    TestDerivative * test30 = new TestDerivative("d2SdForwarddForward", getForward, setForward, d2SdForwarddForward, dSdForward, 0.001, 0.000001, 0.0001, 0.000001);
    this->myTestsDerivatives.addTest(test30);
}

void MainTestPlan::fillTestsIntermediates(){
    computer z = &Sabr::z;

    computerIntermediate d2xzdzdz = &Sabr::d2xzdzdz;
    computerIntermediate dxzdz = &Sabr::dxzdz;
    computerIntermediate xz = &Sabr::xz;

    TestIntermediates * test1 = new TestIntermediates("d2xzdzdz", z, d2xzdzdz, dxzdz, 0.001, 0.00001, 0.00001, 0.0001);
    this->myTestsIntermediates.addTest(test1);
    TestIntermediates * test2 = new TestIntermediates("dxzdz", z, dxzdz, xz, 0.001, 0.00001, 0.00001, 0.0001);
    this->myTestsIntermediates.addTest(test2);
}

void MainTestPlan::fillTestsFallbackAlpha(){
    paramGetter getAlpha = &Sabr::getAlpha;
    paramSetter setAlpha = &Sabr::setAlpha;

    computer d2zxzdAlphadAlpha = &Sabr::d2zxzdAlphadAlpha;
    computer d2zxzdBetadAlpha = &Sabr::d2zxzdBetadAlpha;
    computer d2zxzdRhodAlpha = &Sabr::d2zxzdRhodAlpha;
    computer d2zxzdVovoldAlpha = &Sabr::d2zxzdVovoldAlpha;
    computer d2zxzdForwarddAlpha = &Sabr::d2zxzdForwarddAlpha;
    computer dzxzdAlpha = &Sabr::dzxzdAlpha;
    computer dzxzdBeta = &Sabr::dzxzdBeta;
    computer dzxzdRho = &Sabr::dzxzdRho;
    computer dzxzdVovol = &Sabr::dzxzdVovol;
    computer dzxzdForward = &Sabr::dzxzdForward;
    computer zxz = &Sabr::zxz;
    TestDerivative * test7 = new TestDerivative("dzxzdAlpha", getAlpha, setAlpha, dzxzdAlpha, zxz, 0.001, 0.000001, 0.001, 0.00001);
    this->myTestsFallbackAlpha.addTest(test7);
    TestDerivative * test8 = new TestDerivative("d2zxzdAlphadAlpha", getAlpha, setAlpha, d2zxzdAlphadAlpha, dzxzdAlpha, 0.001, 0.000001, 1, 0.00001);
    this->myTestsFallbackAlpha.addTest(test8);
    TestDerivative * test13 = new TestDerivative("d2zxzdBetadAlpha", getAlpha, setAlpha, d2zxzdBetadAlpha, dzxzdBeta, 0.001, 0.000001, 1, 0.00001);
    this->myTestsFallbackAlpha.addTest(test13);
    TestDerivative * test17 = new TestDerivative("d2zxzdRhodAlpha", getAlpha, setAlpha, d2zxzdRhodAlpha, dzxzdRho, 0.001, 0.000001, 1, 0.00001);
    this->myTestsFallbackAlpha.addTest(test17);
    TestDerivative * test21 = new TestDerivative("d2zxzdVovoldAlpha", getAlpha, setAlpha, d2zxzdVovoldAlpha, dzxzdVovol, 0.001, 0.000001, 1, 0.00001);
    this->myTestsFallbackAlpha.addTest(test21);
    TestDerivative * test25 = new TestDerivative("d2zxzdForwarddAlpha", getAlpha, setAlpha, d2zxzdForwarddAlpha, dzxzdForward, 0.001, 0.000001, 1, 0.00001);
    this->myTestsFallbackAlpha.addTest(test25);

    computer d2xzdAlphadAlpha = &Sabr::d2xzdAlphadAlpha;
    computer dxzdAlpha = &Sabr::dxzdAlpha;
    computer d2xzdBetadAlpha = &Sabr::d2xzdBetadAlpha;
    computer dxzdBeta = &Sabr::dxzdBeta;
    computer d2xzdRhodAlpha = &Sabr::d2xzdRhodAlpha;
    computer dxzdRho = &Sabr::dxzdRho;
    computer d2xzdVovoldAlpha = &Sabr::d2xzdVovoldAlpha;
    computer dxzdVovol = &Sabr::dxzdVovol;
    computer d2xzdForwarddAlpha = &Sabr::d2xzdForwarddAlpha;
    computer dxzdForward = &Sabr::dxzdForward;
    /*
    computer xz = &Sabr::xz;
    TestDerivative * test29 = new TestDerivative("dxzdAlpha", getAlpha, setAlpha, dxzdAlpha, xz, 0.001, 0.000001, 10, 0.00001);
    this->myTestsFallbackAlpha.addTest(test29);
    */
    TestDerivative * test9 = new TestDerivative("d2xzdAlphadAlpha", getAlpha, setAlpha, d2xzdAlphadAlpha, dxzdAlpha, 0.001, 0.000001, 10, 0.00001);
    this->myTestsFallbackAlpha.addTest(test9);
    TestDerivative * test14 = new TestDerivative("d2xzdBetadAlpha", getAlpha, setAlpha, d2xzdBetadAlpha, dxzdBeta, 0.001, 0.000001, 10, 0.00001);
    this->myTestsFallbackAlpha.addTest(test14);
    TestDerivative * test18 = new TestDerivative("d2xzdRhodAlpha", getAlpha, setAlpha, d2xzdRhodAlpha, dxzdRho, 0.001, 0.000001, 10, 0.00001);
    this->myTestsFallbackAlpha.addTest(test18);
    TestDerivative * test22 = new TestDerivative("d2xzdVovoldAlpha", getAlpha, setAlpha, d2xzdVovoldAlpha, dxzdVovol, 0.001, 0.000001, 10, 0.00001);
    this->myTestsFallbackAlpha.addTest(test22);
    TestDerivative * test26 = new TestDerivative("d2xzdForwarddAlpha", getAlpha, setAlpha, d2xzdForwarddAlpha, dxzdForward, 0.001, 0.000001, 10, 0.00001);
    this->myTestsFallbackAlpha.addTest(test26);

    computer d2zdAlphadAlpha = &Sabr::d2zdAlphadAlpha;
    computer dzdAlpha = &Sabr::dzdAlpha;
    computer d2zdBetadAlpha = &Sabr::d2zdBetadAlpha;
    computer dzdBeta = &Sabr::dzdBeta;
    computer d2zdRhodAlpha = &Sabr::d2zdRhodAlpha;
    computer dzdRho = &Sabr::dzdRho;
    computer d2zdVovoldAlpha = &Sabr::d2zdVovoldAlpha;
    computer dzdVovol = &Sabr::dzdVovol;
    computer d2zdForwarddAlpha = &Sabr::d2zdForwarddAlpha;
    computer dzdForward = &Sabr::dzdForward;
    computer z = &Sabr::z;
    TestDerivative * test11 = new TestDerivative("dzdAlpha", getAlpha, setAlpha, dzdAlpha, z, 0.001, 0.0000001, 0.0001, 0.0000001);
    this->myTestsFallbackAlpha.addTest(test11);
    TestDerivative * test10 = new TestDerivative("d2zdAlphadAlpha", getAlpha, setAlpha, d2zdAlphadAlpha, dzdAlpha, 0.001, 0.0000001, 1, 0.000001);
    this->myTestsFallbackAlpha.addTest(test10);
    TestDerivative * test15 = new TestDerivative("d2zdBetadAlpha", getAlpha, setAlpha, d2zdBetadAlpha, dzdBeta, 0.001, 0.0000001, 1, 0.000001);
    this->myTestsFallbackAlpha.addTest(test15);
    TestDerivative * test19 = new TestDerivative("d2zdRhodAlpha", getAlpha, setAlpha, d2zdRhodAlpha, dzdRho, 0.001, 0.0000001, 1, 0.000001);
    this->myTestsFallbackAlpha.addTest(test19);
    TestDerivative * test23 = new TestDerivative("d2zdVovoldAlpha", getAlpha, setAlpha, d2zdVovoldAlpha, dzdVovol, 0.001, 0.0000001, 1, 0.000001);
    this->myTestsFallbackAlpha.addTest(test23);
    TestDerivative * test27 = new TestDerivative("d2zdForwarddAlpha", getAlpha, setAlpha, d2zdForwarddAlpha, dzdForward, 0.001, 0.0000001, 1, 0.000001);
    this->myTestsFallbackAlpha.addTest(test27);

    computer d2RfdAlphadAlpha = &Sabr::d2RfdAlphadAlpha;
    computer dRfdAlpha = &Sabr::dRfdAlpha;
    computer d2RfdBetadAlpha = &Sabr::d2RfdBetadAlpha;
    computer dRfdBeta = &Sabr::dRfdBeta;
    computer d2RfdRhodAlpha = &Sabr::d2RfdRhodAlpha;
    computer dRfdRho = &Sabr::dRfdRho;
    computer d2RfdVovoldAlpha = &Sabr::d2RfdVovoldAlpha;
    computer dRfdVovol = &Sabr::dRfdVovol;
    computer d2RfdForwarddAlpha = &Sabr::d2RfdForwarddAlpha;
    computer dRfdForward = &Sabr::dRfdForward;
    computer Rf = &Sabr::Rf;
    TestDerivative * test30 = new TestDerivative("dRfdAlpha", getAlpha, setAlpha, dRfdAlpha, Rf, 0.001, 0.00001, 0.0001, 0.0001);
    this->myTestsFallbackAlpha.addTest(test30);
    TestDerivative * test12 = new TestDerivative("d2RfdAlphadAlpha", getAlpha, setAlpha, d2RfdAlphadAlpha, dRfdAlpha, 0.001, 0.00001, 0.0001, 0.0001);
    this->myTestsFallbackAlpha.addTest(test12);
    TestDerivative * test16 = new TestDerivative("d2RfdBetadAlpha", getAlpha, setAlpha, d2RfdBetadAlpha, dRfdBeta, 0.001, 0.00001, 0.0001, 0.0001);
    this->myTestsFallbackAlpha.addTest(test16);
    TestDerivative * test20 = new TestDerivative("d2RfdRhodAlpha", getAlpha, setAlpha, d2RfdRhodAlpha, dRfdRho, 0.001, 0.00001, 0.0001, 0.0001);
    this->myTestsFallbackAlpha.addTest(test20);
    TestDerivative * test24 = new TestDerivative("d2RfdVovoldAlpha", getAlpha, setAlpha, d2RfdVovoldAlpha, dRfdVovol, 0.001, 0.00001, 0.0001, 0.0001);
    this->myTestsFallbackAlpha.addTest(test24);
    TestDerivative * test28 = new TestDerivative("d2RfdForwarddAlpha", getAlpha, setAlpha, d2RfdForwarddAlpha, dRfdForward, 0.001, 0.00001, 0.0001, 0.0001);
    this->myTestsFallbackAlpha.addTest(test28);
}

void MainTestPlan::fillTestsFallbackBeta(){
    paramGetter getBeta = &Sabr::getBeta;
    paramSetter setBeta = &Sabr::setBeta;

    computer d2zxzdAlphadBeta = &Sabr::d2zxzdAlphadBeta;
    computer d2zxzdBetadBeta = &Sabr::d2zxzdBetadBeta;
    computer d2zxzdRhodBeta = &Sabr::d2zxzdRhodBeta;
    computer d2zxzdVovoldBeta = &Sabr::d2zxzdVovoldBeta;
    computer d2zxzdForwarddBeta = &Sabr::d2zxzdForwarddBeta;

    computer dzxzdAlpha = &Sabr::dzxzdAlpha;
    computer dzxzdBeta = &Sabr::dzxzdBeta;
    computer dzxzdRho = &Sabr::dzxzdRho;
    computer dzxzdVovol = &Sabr::dzxzdVovol;
    computer dzxzdForward = &Sabr::dzxzdForward;
    computer zxz = &Sabr::zxz;

    TestDerivative * test1 = new TestDerivative("dzxzdBeta", getBeta, setBeta, dzxzdBeta, zxz, 0.001, 0.000001, 1, 0.00001);
    this->myTestsFallbackBeta.addTest(test1);
    TestDerivative * test2 = new TestDerivative("d2zxzdAlphadBeta", getBeta, setBeta, d2zxzdAlphadBeta, dzxzdAlpha, 0.001, 0.000001, 1, 0.00001);
    this->myTestsFallbackBeta.addTest(test2);
    TestDerivative * test3 = new TestDerivative("d2zxzdBetadBeta", getBeta, setBeta, d2zxzdBetadBeta, dzxzdBeta, 0.001, 0.000001, 0.001, 0.00001);
    this->myTestsFallbackBeta.addTest(test3);
    TestDerivative * test4 = new TestDerivative("d2zxzdRhodBeta", getBeta, setBeta, d2zxzdRhodBeta, dzxzdRho, 0.001, 0.000001, 1, 0.00001);
    this->myTestsFallbackBeta.addTest(test4);
    TestDerivative * test5 = new TestDerivative("d2zxzdVovoldBeta", getBeta, setBeta, d2zxzdVovoldBeta, dzxzdVovol, 0.001, 0.000001, 1, 0.00001);
    this->myTestsFallbackBeta.addTest(test5);
    TestDerivative * test6 = new TestDerivative("d2zxzdForwarddBeta", getBeta, setBeta, d2zxzdForwarddBeta, dzxzdForward, 0.001, 0.000001, 1, 0.00001);
    this->myTestsFallbackBeta.addTest(test6);

    computer d2xzdAlphadBeta = &Sabr::d2xzdAlphadBeta;
    computer dxzdAlpha = &Sabr::dxzdAlpha;
    computer d2xzdBetadBeta = &Sabr::d2xzdBetadBeta;
    computer dxzdBeta = &Sabr::dxzdBeta;
    computer d2xzdRhodBeta = &Sabr::d2xzdRhodBeta;
    computer dxzdRho = &Sabr::dxzdRho;
    computer d2xzdVovoldBeta = &Sabr::d2xzdVovoldBeta;
    computer dxzdVovol = &Sabr::dxzdVovol;
    computer d2xzdForwarddBeta = &Sabr::d2xzdForwarddBeta;
    computer dxzdForward = &Sabr::dxzdForward;
    /*
    computer xz = &Sabr::xz;
    TestDerivative * test7 = new TestDerivative("dxzdBeta", getBeta, setBeta, dxzdBeta, xz, 0.001, 0.000001, 10, 0.00001);
    this->myTestsFallbackBeta.addTest(test7);
    */
    TestDerivative * test8 = new TestDerivative("d2xzdAlphadBeta", getBeta, setBeta, d2xzdAlphadBeta, dxzdAlpha, 0.001, 0.000001, 10, 0.00001);
    this->myTestsFallbackBeta.addTest(test8);
    TestDerivative * test9 = new TestDerivative("d2xzdBetadBeta", getBeta, setBeta, d2xzdBetadBeta, dxzdBeta, 0.001, 0.000001, 10, 0.00001);
    this->myTestsFallbackBeta.addTest(test9);
    TestDerivative * test10 = new TestDerivative("d2xzdRhodBeta", getBeta, setBeta, d2xzdRhodBeta, dxzdRho, 0.001, 0.000001, 10, 0.00001);
    this->myTestsFallbackBeta.addTest(test10);
    TestDerivative * test11 = new TestDerivative("d2xzdVovoldBeta", getBeta, setBeta, d2xzdVovoldBeta, dxzdVovol, 0.001, 0.000001, 10, 0.00001);
    this->myTestsFallbackBeta.addTest(test11);
    TestDerivative * test12 = new TestDerivative("d2xzdForwarddBeta", getBeta, setBeta, d2xzdForwarddBeta, dxzdForward, 0.001, 0.000001, 10, 0.00001);
    this->myTestsFallbackBeta.addTest(test12);

    computer d2zdAlphadBeta = &Sabr::d2zdAlphadBeta;
    computer dzdAlpha = &Sabr::dzdAlpha;
    computer d2zdBetadBeta = &Sabr::d2zdBetadBeta;
    computer dzdBeta = &Sabr::dzdBeta;
    computer d2zdRhodBeta = &Sabr::d2zdRhodBeta;
    computer dzdRho = &Sabr::dzdRho;
    computer d2zdVovoldBeta = &Sabr::d2zdVovoldBeta;
    computer dzdVovol = &Sabr::dzdVovol;
    computer d2zdForwarddBeta = &Sabr::d2zdForwarddBeta;
    computer dzdForward = &Sabr::dzdForward;
    computer z = &Sabr::z;
    TestDerivative * test13 = new TestDerivative("dzdBeta", getBeta, setBeta, dzdBeta, z, 0.001, 0.0000001, 0.0001, 0.0000001);
    this->myTestsFallbackBeta.addTest(test13);
    TestDerivative * test14 = new TestDerivative("d2zdAlphadBeta", getBeta, setBeta, d2zdAlphadBeta, dzdAlpha, 0.001, 0.0000001, 1, 0.000001);
    this->myTestsFallbackBeta.addTest(test14);
    TestDerivative * test15 = new TestDerivative("d2zdBetadBeta", getBeta, setBeta, d2zdBetadBeta, dzdBeta, 0.001, 0.0000001, 1, 0.000001);
    this->myTestsFallbackBeta.addTest(test15);
    TestDerivative * test16 = new TestDerivative("d2zdRhodBeta", getBeta, setBeta, d2zdRhodBeta, dzdRho, 0.001, 0.0000001, 1, 0.000001);
    this->myTestsFallbackBeta.addTest(test16);
    TestDerivative * test17 = new TestDerivative("d2zdVovoldBeta", getBeta, setBeta, d2zdVovoldBeta, dzdVovol, 0.001, 0.0000001, 1, 0.000001);
    this->myTestsFallbackBeta.addTest(test17);
    TestDerivative * test18 = new TestDerivative("d2zdForwarddBeta", getBeta, setBeta, d2zdForwarddBeta, dzdForward, 0.001, 0.0000001, 1, 0.000001);
    this->myTestsFallbackBeta.addTest(test18);

    computer d2RfdAlphadBeta = &Sabr::d2RfdAlphadBeta;
    computer dRfdAlpha = &Sabr::dRfdAlpha;
    computer d2RfdBetadBeta = &Sabr::d2RfdBetadBeta;
    computer dRfdBeta = &Sabr::dRfdBeta;
    computer d2RfdRhodBeta = &Sabr::d2RfdRhodBeta;
    computer dRfdRho = &Sabr::dRfdRho;
    computer d2RfdVovoldBeta = &Sabr::d2RfdVovoldBeta;
    computer dRfdVovol = &Sabr::dRfdVovol;
    computer d2RfdForwarddBeta = &Sabr::d2RfdForwarddBeta;
    computer dRfdForward = &Sabr::dRfdForward;
    computer Rf = &Sabr::Rf;
    TestDerivative * test19 = new TestDerivative("dRfdBeta", getBeta, setBeta, dRfdBeta, Rf, 0.001, 0.00001, 0.0001, 0.0001);
    this->myTestsFallbackBeta.addTest(test19);
    TestDerivative * test20 = new TestDerivative("d2RfdAlphadBeta", getBeta, setBeta, d2RfdAlphadBeta, dRfdAlpha, 0.001, 0.00001, 0.0001, 0.0001);
    this->myTestsFallbackBeta.addTest(test20);
    TestDerivative * test21 = new TestDerivative("d2RfdBetadBeta", getBeta, setBeta, d2RfdBetadBeta, dRfdBeta, 0.001, 0.00001, 0.0001, 0.0001);
    this->myTestsFallbackBeta.addTest(test21);
    TestDerivative * test22 = new TestDerivative("d2RfdRhodBeta", getBeta, setBeta, d2RfdRhodBeta, dRfdRho, 0.001, 0.00001, 0.0001, 0.0001);
    this->myTestsFallbackBeta.addTest(test22);
    TestDerivative * test23 = new TestDerivative("d2RfdVovoldBeta", getBeta, setBeta, d2RfdVovoldBeta, dRfdVovol, 0.001, 0.00001, 0.0001, 0.0001);
    this->myTestsFallbackBeta.addTest(test23);
    TestDerivative * test24 = new TestDerivative("d2RfdForwarddBeta", getBeta, setBeta, d2RfdForwarddBeta, dRfdForward, 0.001, 0.00001, 0.0001, 0.0001);
    this->myTestsFallbackBeta.addTest(test24);

    computer W = &Sabr::W;
    computer dWdBeta = &Sabr::dWdBeta;
    computer d2WdBetadBeta = &Sabr::d2WdBetadBeta;
    computer dWdForward = &Sabr::dWdForward;
    computer d2WdForwarddBeta = &Sabr::d2WdForwarddBeta;
    TestDerivative * test25 = new TestDerivative("dWdBeta", getBeta, setBeta, dWdBeta, W, 0.001, 0.00001, 0.0001, 0.0001);
    this->myTestsFallbackBeta.addTest(test25);
    TestDerivative * test26 = new TestDerivative("d2WdBetadBeta", getBeta, setBeta, d2WdBetadBeta, dWdBeta, 0.001, 0.00001, 0.0001, 0.0001);
    this->myTestsFallbackBeta.addTest(test26);
    TestDerivative * test27 = new TestDerivative("d2WdForwarddBeta", getBeta, setBeta, d2WdForwarddBeta, dWdForward, 0.001, 0.00001, 0.0001, 0.00001);
    this->myTestsFallbackBeta.addTest(test27);

    computer Lf = &Sabr::Lf;
    computer dLfdBeta = &Sabr::dLfdBeta;
    computer d2LfdBetadBeta = &Sabr::d2LfdBetadBeta;
    computer dLfdForward = &Sabr::dLfdForward;
    computer d2LfdForwarddBeta = &Sabr::d2LfdForwarddBeta;
    TestDerivative * test28 = new TestDerivative("dLfdBeta", getBeta, setBeta, dLfdBeta, Lf, 0.001, 0.00001, 0.0001, 0.0001);
    this->myTestsFallbackBeta.addTest(test28);
    TestDerivative * test29 = new TestDerivative("d2LfdBetadBeta", getBeta, setBeta, d2LfdBetadBeta, dLfdBeta, 0.001, 0.00001, 0.0001, 0.0001);
    this->myTestsFallbackBeta.addTest(test29);
    TestDerivative * test30 = new TestDerivative("d2LfdForwarddBeta", getBeta, setBeta, d2LfdForwarddBeta, dLfdForward, 0.001, 0.00001, 0.0001, 0.0001);
    this->myTestsFallbackBeta.addTest(test30);
}

void MainTestPlan::fillTestsFallbackRho(){
    paramGetter getRho = &Sabr::getRho;
    paramSetter setRho = &Sabr::setRho;

    computer d2zxzdAlphadRho = &Sabr::d2zxzdAlphadRho;
    computer d2zxzdBetadRho = &Sabr::d2zxzdBetadRho;
    computer d2zxzdRhodRho = &Sabr::d2zxzdRhodRho;
    computer d2zxzdVovoldRho = &Sabr::d2zxzdVovoldRho;
    computer d2zxzdForwarddRho = &Sabr::d2zxzdForwarddRho;

    computer dzxzdAlpha = &Sabr::dzxzdAlpha;
    computer dzxzdBeta = &Sabr::dzxzdBeta;
    computer dzxzdRho = &Sabr::dzxzdRho;
    computer dzxzdVovol = &Sabr::dzxzdVovol;
    computer dzxzdForward = &Sabr::dzxzdForward;
    computer zxz = &Sabr::zxz;

    TestDerivative * test1 = new TestDerivative("dzxzdRho", getRho, setRho, dzxzdRho, zxz, 0.001, 0.000001, 1, 0.00001);
    this->myTestsFallbackRho.addTest(test1);
    TestDerivative * test2 = new TestDerivative("d2zxzdAlphadRho", getRho, setRho, d2zxzdAlphadRho, dzxzdAlpha, 0.001, 0.000001, 1, 0.00001);
    this->myTestsFallbackRho.addTest(test2);
    TestDerivative * test3 = new TestDerivative("d2zxzdBetadRho", getRho, setRho, d2zxzdBetadRho, dzxzdBeta, 0.001, 0.000001, 0.001, 0.00001);
    this->myTestsFallbackRho.addTest(test3);
    TestDerivative * test4 = new TestDerivative("d2zxzdRhodRho", getRho, setRho, d2zxzdRhodRho, dzxzdRho, 0.001, 0.000001, 1, 0.00001);
    this->myTestsFallbackRho.addTest(test4);
    TestDerivative * test5 = new TestDerivative("d2zxzdVovoldRho", getRho, setRho, d2zxzdVovoldRho, dzxzdVovol, 0.001, 0.000001, 1, 0.00001);
    this->myTestsFallbackRho.addTest(test5);
    TestDerivative * test6 = new TestDerivative("d2zxzdForwarddRho", getRho, setRho, d2zxzdForwarddRho, dzxzdForward, 0.001, 0.000001, 1, 0.00001);
    this->myTestsFallbackRho.addTest(test6);

    computer d2xzdAlphadRho = &Sabr::d2xzdAlphadRho;
    computer dxzdAlpha = &Sabr::dxzdAlpha;
    computer d2xzdBetadRho = &Sabr::d2xzdBetadRho;
    computer dxzdBeta = &Sabr::dxzdBeta;
    computer d2xzdRhodRho = &Sabr::d2xzdRhodRho;
    computer dxzdRho = &Sabr::dxzdRho;
    computer d2xzdVovoldRho = &Sabr::d2xzdVovoldRho;
    computer dxzdVovol = &Sabr::dxzdVovol;
    computer d2xzdForwarddRho = &Sabr::d2xzdForwarddRho;
    computer dxzdForward = &Sabr::dxzdForward;
    /*
    computer xz = &Sabr::xz;
    TestDerivative * test7 = new TestDerivative("dxzdBeta", getBeta, setBeta, dxzdBeta, xz, 0.001, 0.000001, 10, 0.00001);
    this->myTestsFallbackRho.addTest(test7);
    */
    TestDerivative * test8 = new TestDerivative("d2xzdAlphadRho", getRho, setRho, d2xzdAlphadRho, dxzdAlpha, 0.001, 0.000001, 10, 0.00001);
    this->myTestsFallbackRho.addTest(test8);
    TestDerivative * test9 = new TestDerivative("d2xzdBetadRho", getRho, setRho, d2xzdBetadRho, dxzdBeta, 0.001, 0.000001, 10, 0.00001);
    this->myTestsFallbackRho.addTest(test9);
    TestDerivative * test10 = new TestDerivative("d2xzdRhodRho", getRho, setRho, d2xzdRhodRho, dxzdRho, 0.001, 0.000001, 10, 0.00001);
    this->myTestsFallbackRho.addTest(test10);
    TestDerivative * test11 = new TestDerivative("d2xzdVovoldRho", getRho, setRho, d2xzdVovoldRho, dxzdVovol, 0.001, 0.000001, 10, 0.00001);
    this->myTestsFallbackRho.addTest(test11);
    TestDerivative * test12 = new TestDerivative("d2xzdForwarddRho", getRho, setRho, d2xzdForwarddRho, dxzdForward, 0.001, 0.000001, 10, 0.00001);
    this->myTestsFallbackRho.addTest(test12);

    computer d2zdAlphadRho = &Sabr::d2zdAlphadRho;
    computer dzdAlpha = &Sabr::dzdAlpha;
    computer d2zdBetadRho = &Sabr::d2zdBetadRho;
    computer dzdBeta = &Sabr::dzdBeta;
    computer d2zdRhodRho = &Sabr::d2zdRhodRho;
    computer dzdRho = &Sabr::dzdRho;
    computer d2zdVovoldRho = &Sabr::d2zdVovoldRho;
    computer dzdVovol = &Sabr::dzdVovol;
    computer d2zdForwarddRho = &Sabr::d2zdForwarddRho;
    computer dzdForward = &Sabr::dzdForward;
    computer z = &Sabr::z;
    TestDerivative * test13 = new TestDerivative("dzdRho", getRho, setRho, dzdRho, z, 0.001, 0.0000001, 0.0001, 0.0000001);
    this->myTestsFallbackRho.addTest(test13);
    TestDerivative * test14 = new TestDerivative("d2zdAlphadRho", getRho, setRho, d2zdAlphadRho, dzdAlpha, 0.001, 0.0000001, 1, 0.000001);
    this->myTestsFallbackRho.addTest(test14);
    TestDerivative * test15 = new TestDerivative("d2zdBetadRho", getRho, setRho, d2zdBetadRho, dzdBeta, 0.001, 0.0000001, 1, 0.000001);
    this->myTestsFallbackRho.addTest(test15);
    TestDerivative * test16 = new TestDerivative("d2zdRhodRho", getRho, setRho, d2zdRhodRho, dzdRho, 0.001, 0.0000001, 1, 0.000001);
    this->myTestsFallbackRho.addTest(test16);
    TestDerivative * test17 = new TestDerivative("d2zdVovoldRho", getRho, setRho, d2zdVovoldRho, dzdVovol, 0.001, 0.0000001, 1, 0.000001);
    this->myTestsFallbackRho.addTest(test17);
    TestDerivative * test18 = new TestDerivative("d2zdForwarddRho", getRho, setRho, d2zdForwarddRho, dzdForward, 0.001, 0.0000001, 1, 0.000001);
    this->myTestsFallbackRho.addTest(test18);

    computer d2RfdAlphadRho = &Sabr::d2RfdAlphadRho;
    computer dRfdAlpha = &Sabr::dRfdAlpha;
    computer d2RfdBetadRho = &Sabr::d2RfdBetadRho;
    computer dRfdBeta = &Sabr::dRfdBeta;
    computer d2RfdRhodRho = &Sabr::d2RfdRhodRho;
    computer dRfdRho = &Sabr::dRfdRho;
    computer d2RfdVovoldRho = &Sabr::d2RfdVovoldRho;
    computer dRfdVovol = &Sabr::dRfdVovol;
    computer d2RfdForwarddRho = &Sabr::d2RfdForwarddRho;
    computer dRfdForward = &Sabr::dRfdForward;
    computer Rf = &Sabr::Rf;
    TestDerivative * test19 = new TestDerivative("dRfdRho", getRho, setRho, dRfdRho, Rf, 0.001, 0.00001, 0.0001, 0.0001);
    this->myTestsFallbackRho.addTest(test19);
    TestDerivative * test20 = new TestDerivative("d2RfdAlphadRho", getRho, setRho, d2RfdAlphadRho, dRfdAlpha, 0.001, 0.00001, 0.0001, 0.0001);
    this->myTestsFallbackRho.addTest(test20);
    TestDerivative * test21 = new TestDerivative("d2RfdBetadRho", getRho, setRho, d2RfdBetadRho, dRfdBeta, 0.001, 0.00001, 0.0001, 0.0001);
    this->myTestsFallbackRho.addTest(test21);
    TestDerivative * test22 = new TestDerivative("d2RfdRhodRho", getRho, setRho, d2RfdRhodRho, dRfdRho, 0.001, 0.00001, 0.0001, 0.0001);
    this->myTestsFallbackRho.addTest(test22);
    TestDerivative * test23 = new TestDerivative("d2RfdVovoldRho", getRho, setRho, d2RfdVovoldRho, dRfdVovol, 0.001, 0.00001, 0.0001, 0.0001);
    this->myTestsFallbackRho.addTest(test23);
    TestDerivative * test24 = new TestDerivative("d2RfdForwarddRho", getRho, setRho, d2RfdForwarddRho, dRfdForward, 0.001, 0.00001, 0.0001, 0.0001);
    this->myTestsFallbackRho.addTest(test24);
}

void MainTestPlan::fillTestsFallbackVovol(){}

void MainTestPlan::fillTestsFallbackForward(){}

void MainTestPlan::run(){
    // Run the tests on first order derivatives
    this->myTestsDerivatives.run();

    // Run the tests on intermediates
    this->myTestsIntermediates.run();

    // this->myTestsFallbackAlpha.run();
    // this->myTestsFallbackBeta.run();
    // this->myTestsFallbackRho.run();
    // this->myTestsFallbackVovol.run();
    // this->myTestsFallbackForward.run();
}
