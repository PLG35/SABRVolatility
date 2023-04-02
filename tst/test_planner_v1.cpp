#include "test_smile_curve.h"
#include <iostream>

void launch_tests(){

    Test_Sabr myTest_Sabr = Test_Sabr();
    double money = 0.001;

    // TESTING Z
    cout << "Test z" << endl;

    if(myTest_Sabr.test_dzdAlpha(money, 0.00001, 0.0001)){
        cout << "Test Sabr::dzdAlpha() : PASSED" << endl;
    }
    else{
        cout << "Test Sabr::dzdAlpha() : FAILED" << endl;
    }

    // z has convexity wrt beta
    if(myTest_Sabr.test_dzdBeta(money, 0.00001, 0.001)){
        cout << "Test Sabr::dzdBeta() : PASSED" << endl;
    }
    else{
        cout << "Test Sabr::dzdBeta() : FAILED" << endl;
    }

    if(myTest_Sabr.test_dzdVovol(money, 0.00001, 0.0001)){
        cout << "Test Sabr::dzdVovol() : PASSED" << endl;
    }
    else{
        cout << "Test Sabr::dzdVovol() : FAILED" << endl;
    }

    // TESTING XZ
    cout << "Test xz" << endl;

    if(myTest_Sabr.test_dxzdz(money, 0.00001, 0.0001)){
        cout << "Test Sabr::dxzdz() : PASSED" << endl;
    }
    else{
        cout << "Test Sabr::dxzdz() : FAILED" << endl;
    }

    if(myTest_Sabr.test_dxzdAlpha(money, 0.00001, 0.0001)){
        cout << "Test Sabr::dxzdAlpha() : PASSED" << endl;
    }
    else{
        cout << "Test Sabr::dxzdAlpha() : FAILED" << endl;
    }

    if(myTest_Sabr.test_dxzdBeta(money, 0.00001, 0.0001)){
        cout << "Test Sabr::dxzdBeta() : PASSED" << endl;
    }
    else{
        cout << "Test Sabr::dxzdBeta() : FAILED" << endl;
    }

    if(myTest_Sabr.test_dxzdRho(money, 0.00001, 0.0001)){
        cout << "Test Sabr::dxzdRho() : PASSED" << endl;
    }
    else{
        cout << "Test Sabr::dxzdRho() : FAILED" << endl;
    }

    if(myTest_Sabr.test_dxzdVovol(money, 0.00001, 0.0001)){
        cout << "Test Sabr::dxzdVovol() : PASSED" << endl;
    }
    else{
        cout << "Test Sabr::dxzdVovol() : FAILED" << endl;
    }

    // TESTING ZXZ
    cout << "Test zxz" << endl;

    if(myTest_Sabr.test_dzxzdAlpha(money, 0.00001, 0.0001)){
        cout << "Test Sabr::dzxzdAlpha() general case : PASSED" << endl;
    }
    else{
        cout << "Test Sabr::dzxzdAlpha() general case : FAILED" << endl;
    }

    if(myTest_Sabr.test_dzxzdBeta(money, 0.00001, 0.0001)){
        cout << "Test Sabr::dzxzdBeta() general case : PASSED" << endl;
    }
    else{
        cout << "Test Sabr::dzxzdBeta() general case : FAILED" << endl;
    }

    if(myTest_Sabr.test_dzxzdRho(money, 0.00001, 0.0001)){
        cout << "Test Sabr::dzxzdRho() general case : PASSED" << endl;
    }
    else{
        cout << "Test Sabr::dzxzdRho() general case : FAILED" << endl;
    }

    if(myTest_Sabr.test_dzxzdVovol(money, 0.00001, 0.0001)){
        cout << "Test Sabr::dzxzdVovol() general case : PASSED" << endl;
    }
    else{
        cout << "Test Sabr::dzxzdVovol() general case : FAILED" << endl;
    }

    if(myTest_Sabr.test_dzxzdForward(money, 0.00001, 0.0001)){
        cout << "Test Sabr::dzxzdForward() general case : PASSED" << endl;
    }
    else{
        cout << "Test Sabr::dzxzdForward() general case : FAILED" << endl;
    }

    // TESTS W
    cout << "Test W" << endl;

    if(myTest_Sabr.test_dWdBeta(money, 0.00001, 0.0001)){
        cout << "Test Sabr::dWdBeta() general case : PASSED" << endl;
    }
    else{
        cout << "Test Sabr::dWdBeta() general case : FAILED" << endl;
    }

    if(myTest_Sabr.test_dWdForward(money, 0.00001, 0.0001)){
        cout << "Test Sabr::dWdForward() general case : PASSED" << endl;
    }
    else{
        cout << "Test Sabr::dWdForward() general case : FAILED" << endl;
    }

    if(myTest_Sabr.test_dLfdBeta(money, 0.00001, 0.001)){
        cout << "Test Sabr::dLfdBeta() general case : PASSED" << endl;
    }
    else{
        cout << "Test Sabr::dLfdBeta() general case : FAILED" << endl;
    }

    if(myTest_Sabr.test_dLfdForward(money, 0.00001, 0.001)){
        cout << "Test Sabr::dLfdForward() general case : PASSED" << endl;
    }
    else{
        cout << "Test Sabr::dLfdForward() general case : FAILED" << endl;
    }

    // TESTS Rf
    cout << "Test Rf" << endl;

    if(myTest_Sabr.test_dRfdAlpha(money, 0.00001, 0.0001)){
        cout << "Test Sabr::dRfdAlpha() general case : PASSED" << endl;
    }
    else{
        cout << "Test Sabr::dRfdAlpha() general case : FAILED" << endl;
    }

    if(myTest_Sabr.test_dRfdBeta(money, 0.00001, 0.0001)){
        cout << "Test Sabr::dRfdBeta() general case : PASSED" << endl;
    }
    else{
        cout << "Test Sabr::dRfdBeta() general case : FAILED" << endl;
    }

    if(myTest_Sabr.test_dRfdRho(money, 0.00001, 0.0001)){
        cout << "Test Sabr::dRfdRho() general case : PASSED" << endl;
    }
    else{
        cout << "Test Sabr::dRfdRho() general case : FAILED" << endl;
    }

    if(myTest_Sabr.test_dRfdVovol(money, 0.00001, 0.0001)){
        cout << "Test Sabr::dRfdVovol() general case : PASSED" << endl;
    }
    else{
        cout << "Test Sabr::dRfdVovol() general case : FAILED" << endl;
    }

    if(myTest_Sabr.test_dRfdForward(money, 0.00001, 0.0001)){
        cout << "Test Sabr::dRfdForward() general case : PASSED" << endl;
    }
    else{
        cout << "Test Sabr::dRfdForward() general case : FAILED" << endl;
    }

    // Tests volatility
    cout << "Test Volatility" << endl;

    if(myTest_Sabr.test_dSdAlpha(money, 0.00001, 0.0001)){
        cout << "Test Sabr::dSdAlpha() general case : PASSED" << endl;
    }
    else{
        cout << "Test Sabr::dSdAlpha() general case : FAILED" << endl;
    }

    if(myTest_Sabr.test_dSdBeta(money, 0.00001, 0.0001)){
        cout << "Test Sabr::dSdBeta() general case : PASSED" << endl;
    }
    else{
        cout << "Test Sabr::dSdBeta() general case : FAILED" << endl;
    }

    if(myTest_Sabr.test_dSdRho(money, 0.00001, 0.0001)){
        cout << "Test Sabr::dSdRho() general case : PASSED" << endl;
    }
    else{
        cout << "Test Sabr::dSdRho() general case : FAILED" << endl;
    }

    if(myTest_Sabr.test_dSdVovol(money, 0.00001, 0.0001)){
        cout << "Test Sabr::dSdVovol() general case : PASSED" << endl;
    }
    else{
        cout << "Test Sabr::dSdVovol() general case : FAILED" << endl;
    }

    if(myTest_Sabr.test_dSdForward(money, 0.00001, 0.0001)){
        cout << "Test Sabr::dSdForward() general case : PASSED" << endl;
    }
    else{
        cout << "Test Sabr::dSdForward() general case : FAILED" << endl;
    }
}
