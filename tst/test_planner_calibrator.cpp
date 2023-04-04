#include "test_planner_calibrator.h";

///////////////////
// TEST GRADIENT //
///////////////////

TestGradient::TestGradient(){}

TestGradient::TestGradient(Sabr* mS, string n, paramGetter pG, paramSetter pS, int pIG, double rB, double aT, double rT){
    this->mySabr = mS;
    this->name = n;
    this->pGetter = pG;
    this->pSetter = pS;
    this->positionInGradient = pIG;
    this->relativeBump = rB;
    this->absoluteThreshold = aT;
    this->relativeThreshold = rT;
}

TestGradient::~TestGradient(){}

bool TestGradient::launch(){
    // Initial value
    double initialError = this->computeError();

    // Bump parameter
    double parameter = (this->mySabr->*(this->pGetter))();
    double dparameter = parameter * this->relativeBump;
    (this->mySabr->*(this->pSetter))(parameter + dparameter);

    // Compute variations
    double bumpError = this->computeError();
    double variation = bumpError - initialError;
    gradient gradOfError = {0};
    this->gradientOfError(&gradOfError);
    double sensitivity = gradOfError[this->positionInGradient] * dparameter;
    double relativeError = std::abs(sensitivity / variation - 1);
    double absoluteError = std::abs(sensitivity - variation);

    // Revert
    (this->mySabr->*(this->pSetter))(parameter + dparameter);

    // Test
    bool testResult = relativeError < this->relativeThreshold && absoluteError < this->absoluteThreshold;
    if (testResult || (abs(variation) <= precisionMachine && abs(sensitivity) <= precisionMachine)){
        std::cout << "Test " + this->name + " passed." << std::endl;
    }
    else{
        std::cout << "Test " + this->name + " failed." << std::endl;
        std::cout << "   Gradient value : " << sensitivity << "." << std::endl;
        std::cout << "   Variation reached : " << variation << "." << std::endl;
        std::cout << "   Relative bump by " << (this->relativeBump * 100) << "%." << std::endl;
        std::cout << "   Relative error expected : " << (this->relativeThreshold * 100) << "%." << std::endl;
        std::cout << "   Relative error reached : " << (relativeError * 100) << "%." << std::endl;
        std::cout << "   Absolute error expected : " << this->absoluteThreshold << "." << std::endl;
        std::cout << "   Absolute error reached : " << absoluteError << "." << std::endl;
    }
    return testResult;
}

//////////////////
// TEST HESSIAN //
//////////////////

TestHessian::TestHessian(){}

TestHessian::TestHessian(Sabr* mS, string n, paramGetter pG, paramSetter pS, int pIG, int c, double rB, double aT, double rT){
    this->mySabr = mS;
    this->name = n;
    this->pGetter = pG;
    this->pSetter = pS;
    this->positionInGradient = pIG;
    this->column = c;
    this->relativeBump = rB;
    this->absoluteThreshold = aT;
    this->relativeThreshold = rT;
}

TestHessian::~TestHessian(){}

bool TestHessian::launch(){
    // Initial value
    gradient gradInitial;
    this->gradientOfError(&gradInitial);

    // Bump parameter
    double parameter = (this->mySabr->*(this->pGetter))();
    double dparameter = parameter * this->relativeBump;
    (this->mySabr->*(this->pSetter))(parameter + dparameter);

    // Compute variations
    gradient gradBumped;
    this->gradientOfError(&gradBumped);

    // Compute sensitivity
    hessian hess;
    this->hessianOfError(&hess);

    double sensitivity = hess[this->positionInGradient][this->column] * dparameter;
    double variation = gradBumped[this->column] - gradInitial[this->column];
    double relativeError = std::abs(sensitivity / variation - 1);
    double absoluteError = std::abs(sensitivity - variation);

    // Revert
    (this->mySabr->*(this->pSetter))(parameter + dparameter);

    // Test
    bool testResult = relativeError < this->relativeThreshold && absoluteError < this->absoluteThreshold;
    if (testResult || (abs(variation) <= precisionMachine && abs(sensitivity) <= precisionMachine)){
        std::cout << "Test " + this->name + " passed." << std::endl;
    }
    else{
        std::cout << "Test " + this->name + " failed." << std::endl;
        std::cout << "   Hessian value : " << sensitivity << "." << std::endl;
        std::cout << "   Variation reached : " << variation << "." << std::endl;
        std::cout << "   Relative bump by " << (this->relativeBump * 100) << "%." << std::endl;
        std::cout << "   Relative error expected : " << (this->relativeThreshold * 100) << "%." << std::endl;
        std::cout << "   Relative error reached : " << (relativeError * 100) << "%." << std::endl;
        std::cout << "   Absolute error expected : " << this->absoluteThreshold << "." << std::endl;
        std::cout << "   Absolute error reached : " << absoluteError << "." << std::endl;
    }
    return testResult;
}

//////////////////////////
// TEST INVERSE HESSIAN //
//////////////////////////

TestInverseHessian::TestInverseHessian(){}

TestInverseHessian::TestInverseHessian(Sabr* mS, double aT){
    this->mySabr = mS;
    this->absoluteThreshold = aT;
}

TestInverseHessian::~TestInverseHessian(){}

void TestInverseHessian::launch(){
    hessian hess;
    this->hessianOfError(&hess);

    hessian invHess;
    this->inverseHessianOfError(&invHess);

    hessian product;
    for(int i=0; i<numberOfParameters; i++){
        for(int j=0; j<numberOfParameters; j++){
            product[i][j] = 0;
            for(int k=0; k<numberOfParameters; k++){
                product[i][j] += hess[i][k] * invHess[k][j];
            }
        }
    }

    this->testElement(&product, 0, 0, this->absoluteThreshold);
    this->testElement(&product, 0, 1, this->absoluteThreshold);
    this->testElement(&product, 0, 2, this->absoluteThreshold);
    this->testElement(&product, 0, 3, this->absoluteThreshold);

    this->testElement(&product, 1, 0, this->absoluteThreshold);
    this->testElement(&product, 1, 1, this->absoluteThreshold);
    this->testElement(&product, 1, 2, this->absoluteThreshold);
    this->testElement(&product, 1, 3, this->absoluteThreshold);

    this->testElement(&product, 2, 0, this->absoluteThreshold);
    this->testElement(&product, 2, 1, this->absoluteThreshold);
    this->testElement(&product, 2, 2, this->absoluteThreshold);
    this->testElement(&product, 2, 3, this->absoluteThreshold);

    this->testElement(&product, 3, 0, this->absoluteThreshold);
    this->testElement(&product, 3, 1, this->absoluteThreshold);
    this->testElement(&product, 3, 2, this->absoluteThreshold);
    this->testElement(&product, 3, 3, this->absoluteThreshold);
}

bool TestInverseHessian::testElement(hessian* hess, int line, int column, double aT){
    // Values
    double reachedValue = (*hess)[line][column];
    double targetValue = line == column ? 1 : 0;
    double absoluteError = std::abs(targetValue - reachedValue);

    // Test
    bool testResult = absoluteError < aT;
    if (testResult){
        std::cout << "Test InverseHessian " << line << "x" << column << " passed." << std::endl;
    }
    else{
        std::cout << "Test InverseHessian failed." << std::endl;
        std::cout << "   Product value : " << reachedValue << "." << std::endl;
        std::cout << "   Absolute error expected : " << aT << "." << std::endl;
    }
    return testResult;
}

/////////////////////
// TEST CALIBRATOR //
/////////////////////

TestCalibrator::TestCalibrator(){
    this->mySabr = new Sabr();
    this->targetMatchingPrecision = 1;
    this->limitOfIterations = 10;

    this->fillTestsGradient();
    this->fillTestsHessian();
    this->myTestInverseHessian = new TestInverseHessian(this->mySabr, this->absoluteThresholdForInverseHessian);
}

TestCalibrator::TestCalibrator(double tMP, int lI){
    this->mySabr = new Sabr();
    this->targetMatchingPrecision = tMP;
    this->limitOfIterations = lI;

    this->fillTestsGradient();
    this->fillTestsHessian();
    this->myTestInverseHessian = new TestInverseHessian(this->mySabr, this->absoluteThresholdForInverseHessian);
}

void TestCalibrator::fillTestsGradient(){
    paramGetter getAlpha = &Sabr::getAlpha;
    paramSetter setAlpha = &Sabr::setAlpha;
    paramGetter getBeta = &Sabr::getBeta;
    paramSetter setBeta = &Sabr::setBeta;
    paramGetter getRho = &Sabr::getRho;
    paramSetter setRho = &Sabr::setRho;
    paramGetter getVovol = &Sabr::getVovol;
    paramSetter setVovol = &Sabr::setVovol;

    TestGradient * test1 = new TestGradient(this->mySabr, "gradient on Alpha", getAlpha, setAlpha, 0, 0.000001, 0.001, 0.000001);
    this->myTestsGradient.push_back(test1);
    TestGradient * test2 = new TestGradient(this->mySabr, "gradient on Beta", getBeta, setBeta, 1, 0.000001, 0.001, 0.000001);
    this->myTestsGradient.push_back(test2);
    TestGradient * test3 = new TestGradient(this->mySabr, "gradient on Rho", getRho, setRho, 2, 0.000001, 0.001, 0.000001);
    this->myTestsGradient.push_back(test3);
    TestGradient * test4 = new TestGradient(this->mySabr, "gradient on Vovol", getVovol, setVovol, 3, 0.000001, 0.001, 0.000001);
    this->myTestsGradient.push_back(test4);
}

void TestCalibrator::fillTestsHessian(){
    paramGetter getAlpha = &Sabr::getAlpha;
    paramSetter setAlpha = &Sabr::setAlpha;
    paramGetter getBeta = &Sabr::getBeta;
    paramSetter setBeta = &Sabr::setBeta;
    paramGetter getRho = &Sabr::getRho;
    paramSetter setRho = &Sabr::setRho;
    paramGetter getVovol = &Sabr::getVovol;
    paramSetter setVovol = &Sabr::setVovol;

    TestHessian * test1 = new TestHessian(this->mySabr, "hessian on gradient variation on alpha induced by a variation of alpha", getAlpha, setAlpha, 0, 0, 0.000001, 0.01, 0.00001);
    this->myTestsHessian.push_back(test1);
    TestHessian * test2 = new TestHessian(this->mySabr, "hessian on gradient variation on beta induced by a variation of alpha", getBeta, setBeta, 1, 0, 0.000001, 0.01, 0.00001);
    this->myTestsHessian.push_back(test2);
    TestHessian * test3 = new TestHessian(this->mySabr, "hessian on gradient variation on rho induced by a variation of alpha", getRho, setRho, 2, 0, 0.000001, 0.01, 0.00001);
    this->myTestsHessian.push_back(test3);
    TestHessian * test4 = new TestHessian(this->mySabr, "hessian on gradient variation on vovol induced by a variation of alpha", getVovol, setVovol, 3, 0, 0.000001, 0.01, 0.00001);
    this->myTestsHessian.push_back(test4);
}

TestCalibrator::~TestCalibrator(){
    for(std::vector<TestGradient*>::iterator it = this->myTestsGradient.begin(); it != this->myTestsGradient.end(); ++it){
        delete *it;
    }

    for(std::vector<TestHessian*>::iterator it = this->myTestsHessian.begin(); it != this->myTestsHessian.end(); ++it){
        delete *it;
    }
    delete this->myTestInverseHessian;
}

void TestCalibrator::fillStructure(vector<double> *s, int calibrationDate, int optionMaturity){
    this->strikes = s;
    this->mySabr->setCalDate(calibrationDate);
    this->mySabr->setMaturity(optionMaturity);
    this->numberOfPoints = s->size();

    // Update underlying tests
    for(std::vector<TestGradient*>::iterator it = this->myTestsGradient.begin(); it != this->myTestsGradient.end(); ++it){
        (*it)->fillStructure(s, calibrationDate, optionMaturity);
    }
    for(std::vector<TestHessian*>::iterator it = this->myTestsHessian.begin(); it != this->myTestsHessian.end(); ++it){
        (*it)->fillStructure(s, calibrationDate, optionMaturity);
    }
    this->myTestInverseHessian->fillStructure(s, calibrationDate, optionMaturity);
}

void TestCalibrator::fillMarketData(double forward, vector<double> *tV){
    this->targetVolatilities = tV;
    this->mySabr->setForward(forward);

    // Update underlying tests
    for(std::vector<TestGradient*>::iterator it = this->myTestsGradient.begin(); it != this->myTestsGradient.end(); ++it){
        (*it)->fillMarketData(forward, tV);
    }
    for(std::vector<TestHessian*>::iterator it = this->myTestsHessian.begin(); it != this->myTestsHessian.end(); ++it){
        (*it)->fillMarketData(forward, tV);
    }
    this->myTestInverseHessian->fillMarketData(forward, tV);
}

void TestCalibrator::run(){
    bool testOutcome;

    // Test Gradient
    cout << "START TEST PLAN GRADIENT" << endl;
    for(std::vector<TestGradient*>::iterator it = this->myTestsGradient.begin(); it != this->myTestsGradient.end(); ++it){
        testOutcome = (*it)->launch();
    }
    cout << "END TEST PLAN GRADIENT" << endl;

    // Test Hessian
    cout << "START TEST PLAN HESSIAN" << endl;
    for(std::vector<TestHessian*>::iterator it = this->myTestsHessian.begin(); it != this->myTestsHessian.end(); ++it){
        testOutcome = (*it)->launch();
    }
    cout << "END TEST PLAN HESSIAN" << endl;

    // Test Inverse
    cout << "START TEST PLAN INVERSE HESSIAN" << endl;
    this->myTestInverseHessian->launch();
    cout << "END TEST PLAN INVERSE HESSIAN" << endl;

    // Test reduction of error
}
