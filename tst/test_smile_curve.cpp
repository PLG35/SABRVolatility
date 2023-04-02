#include "test_smile_curve.h"
#include <iostream>
#include <math.h>
#include <cmath>

using namespace std;

// Constructor & destructor
Test_Sabr::Test_Sabr() : Sabr(){};
Test_Sabr::~Test_Sabr(){};

// Tests Z
bool Test_Sabr::test_dzdAlpha(double money, double bumpRatio, double precision){
    double strike = this->forward + money;
    double sensitivity = this->dzdAlpha(strike);

    double z1 = this->z(strike);

    double dalpha = this->alpha * bumpRatio;
    this->setAlpha(this->alpha + dalpha);
    double z2 = this->z(strike);

    double variation = (z2 - z1) / dalpha;
    
    // Restore alpha
    this->setAlpha(this->alpha - dalpha);
    double error = std::abs(sensitivity / variation - 1);

    return  error < precision;
}

bool Test_Sabr::test_dzdBeta(double money, double bumpRatio, double precision){
    double strike = this->forward + money;
    double sensitivity = this->dzdBeta(strike);

    double z1 = this->z(strike);

    double dbeta = this->beta * bumpRatio;
    this->setBeta(this->beta + dbeta);
    double z2 = this->z(strike);

    double variation = (z2 - z1) / dbeta;
    
    // Restore beta
    this->setBeta(this->beta - dbeta);

    double error = std::abs(sensitivity / variation - 1);
    return error < precision;
}

bool Test_Sabr::test_dzdVovol(double money, double bumpRatio, double precision){
    double strike = this->forward + money;
    double sensitivity = this->dzdVovol(strike);

    double z1 = this->z(strike);

    double dvovol = this->vovol * bumpRatio;
    this->setVovol(this->vovol + dvovol);
    double z2 = this->z(strike);

    double variation = (z2 - z1) / dvovol;

    // Restore beta
    this->setVovol(this->vovol - dvovol);

    double error = std::abs(sensitivity / variation - 1);
    return error < precision;
}

// Tests XZ
bool Test_Sabr::test_dxzdz(double money, double bumpRatio, double precision){
    double strike = this->forward + money;

    double z1 = this->z(strike);
    double sensitivity = this->dxzdz(strike, z1);
    double xz1 = this->xz(strike, z1);

    double dstrike = strike * bumpRatio;
    double z2 = this->z(strike + dstrike);
    double xz2 = this->xz(strike + dstrike, z2);

    double variation = (xz2 - xz1) / (z2 - z1);

    double error = std::abs(sensitivity / variation - 1);
    return error < precision;
}

bool Test_Sabr::test_dxzdAlpha(double money, double bumpRatio, double precision){
    double strike = this->forward + money;
    double sensitivity = this->dxzdAlpha(strike);

    double z1 = this->z(strike);
    double xz1 = this->xz(strike, z1);

    double dalpha = alpha * bumpRatio;
    this->setAlpha(this->alpha + dalpha);
    double z2 = this->z(strike);
    double xz2 = this->xz(strike, z2);

    double variation = (xz2 - xz1) / dalpha;

    // Restore beta
    this->setAlpha(this->alpha - dalpha);

    double error = std::abs(sensitivity / variation - 1);
    return error < precision;
}

bool Test_Sabr::test_dxzdBeta(double money, double bumpRatio, double precision){
    double strike = this->forward + money;
    double sensitivity = this->dxzdBeta(strike);

    double z1 = this->z(strike);
    double xz1 = this->xz(strike, z1);

    double dbeta = beta * bumpRatio;
    this->setBeta(this->beta + dbeta);
    double z2 = this->z(strike);
    double xz2 = this->xz(strike, z2);

    double variation = (xz2 - xz1) / dbeta;

    // Restore beta
    this->setBeta(this->beta - dbeta);

    double error = std::abs(sensitivity / variation - 1);
    return error < precision;
}

bool Test_Sabr::test_dxzdRho(double money, double bumpRatio, double precision){
    double strike = this->forward + money;
    double sensitivity = this->dxzdRho(strike);

    double z1 = this->z(strike);
    double xz1 = this->xz(strike, z1);

    double drho = rho * bumpRatio;
    this->setRho(this->rho + drho);
    double z2 = this->z(strike);
    double xz2 = this->xz(strike, z2);

    double variation = (xz2 - xz1) / drho;

    // Restore beta
    this->setRho(this->rho - drho);

    double error = std::abs(sensitivity / variation - 1);
    return error < precision;
}

bool Test_Sabr::test_dxzdVovol(double money, double bumpRatio, double precision){
    double strike = this->forward + money;
    double sensitivity = this->dxzdVovol(strike);

    double z1 = this->z(strike);
    double xz1 = this->xz(strike, z1);

    double dvovol = vovol * bumpRatio;
    this->setVovol(this->vovol + dvovol);
    double z2 = this->z(strike);
    double xz2 = this->xz(strike, z2);

    double variation = (xz2 - xz1) / dvovol;

    // Restore beta
    this->setVovol(this->vovol - dvovol);

    double error = std::abs(sensitivity / variation - 1);
    return error < precision;
}

// Tests ZXZ
bool Test_Sabr::test_dzxzdAlpha(double money, double bumpRatio, double precision){
    double strike = this->forward + money;
    double sensitivity = this->dzxzdAlpha(strike);

    double zxz1 = this->zxz(strike);

    double dalpha = alpha * bumpRatio;
    this->setAlpha(this->alpha + dalpha);
    double zxz2 = this->zxz(strike);

    double variation = (zxz2 - zxz1) / dalpha;

    // Restore beta
    this->setAlpha(this->alpha - dalpha);

    double error = std::abs(sensitivity / variation - 1);
    return error < precision;
}

bool Test_Sabr::test_dzxzdBeta(double money, double bumpRatio, double precision){
    double strike = this->forward + money;
    double sensitivity = this->dzxzdBeta(strike);

    double zxz1 = this->zxz(strike);

    double dbeta = beta * bumpRatio;
    this->setBeta(this->beta + dbeta);
    double zxz2 = this->zxz(strike);

    double variation = (zxz2 - zxz1) / dbeta;

    // Restore beta
    this->setBeta(this->beta - dbeta);

    double error = std::abs(sensitivity / variation - 1);
    return error < precision;
}

bool Test_Sabr::test_dzxzdRho(double money, double bumpRatio, double precision){
    double strike = this->forward + money;
    double sensitivity = this->dzxzdRho(strike);

    double zxz1 = this->zxz(strike);

    double drho = rho * bumpRatio;
    this->setRho(this->rho + drho);
    double zxz2 = this->zxz(strike);

    double variation = (zxz2 - zxz1) / drho;

    // Restore beta
    this->setRho(this->rho - drho);

    double error = std::abs(sensitivity / variation - 1);
    return error < precision;
}

bool Test_Sabr::test_dzxzdVovol(double money, double bumpRatio, double precision){
    double strike = this->forward + money;
    double sensitivity = this->dzxzdVovol(strike);

    double zxz1 = this->zxz(strike);

    double dvovol = vovol * bumpRatio;
    this->setVovol(this->vovol + dvovol);
    double zxz2 = this->zxz(strike);

    double variation = (zxz2 - zxz1) / dvovol;

    // Restore beta
    this->setVovol(this->vovol - dvovol);

    double error = std::abs(sensitivity / variation - 1);
    return error < precision;
}

bool Test_Sabr::test_dzxzdForward(double money, double bumpRatio, double precision){
    double strike = this->forward + money;
    double sensitivity = this->dzxzdForward(strike);

    double zxz1 = this->zxz(strike);

    double dForward = this->forward * bumpRatio;
    this->setForward(this->forward + dForward);
    double zxz2 = this->zxz(strike);

    double variation = (zxz2 - zxz1) / dForward;

    // Restore beta
    this->setForward(this->forward - dForward);

    double error = std::abs(sensitivity / variation - 1);
    return error < precision;
}

// Tests W
bool Test_Sabr::test_dWdBeta(double money, double bumpRatio, double precision){
    double strike = this->forward + money;
    double sensitivity = this->dWdBeta(strike);

    double W1 = this->W(strike);

    double dbeta = beta * bumpRatio;
    this->setBeta(this->beta + dbeta);
    double W2 = this->W(strike);

    double variation = (W2 - W1) / dbeta;

    // Restore beta
    this->setBeta(this->beta - dbeta);

    double error = std::abs(sensitivity / variation - 1);
    return error < precision;
}

bool Test_Sabr::test_dWdForward(double money, double bumpRatio, double precision){
    double strike = this->forward + money;
    double sensitivity = this->dWdForward(strike);

    double W1 = this->W(strike);

    double dForward = this->forward * bumpRatio;
    this->setForward(this->forward + dForward);
    double W2 = this->W(strike);

    double variation = (W2 - W1) / dForward;

    // Restore beta
    this->setForward(this->forward - dForward);

    double error = std::abs(sensitivity / variation - 1);
    return error < precision;
}

// Tests Lf
bool Test_Sabr::test_dLfdBeta(double money, double bumpRatio, double precision){
    double strike = this->forward + money;
    double sensitivity = this->dLfdBeta(strike);

    double Lf1 = this->Lf(strike);

    double dbeta = beta * bumpRatio;
    this->setBeta(this->beta + dbeta);
    double Lf2 = this->Lf(strike);

    double variation = (Lf2 - Lf1) / dbeta;

    // Restore beta
    this->setBeta(this->beta - dbeta);

    double error = std::abs(sensitivity / variation - 1);
    return error < precision;
}

bool Test_Sabr::test_dLfdForward(double money, double bumpRatio, double precision){
    double strike = this->forward + money;
    double sensitivity = this->dLfdForward(strike);

    double Lf1 = this->Lf(strike);

    double dForward = forward * bumpRatio;
    this->setForward(this->forward + dForward);
    double Lf2 = this->Lf(strike);

    double variation = (Lf2 - Lf1) / dForward;

    // Restore beta
    this->setForward(this->forward - dForward);

    double error = std::abs(sensitivity / variation - 1);
    return error < precision;
}

// Tests Rf
bool Test_Sabr::test_dRfdAlpha(double money, double bumpRatio, double precision){
    double strike = this->forward + money;
    double sensitivity = this->dRfdAlpha(strike);

    double Rf1 = this->Rf(strike);

    double dalpha = alpha * bumpRatio;
    this->setAlpha(this->alpha + dalpha);
    double Rf2 = this->Rf(strike);

    double variation = (Rf2 - Rf1) / dalpha;

    // Restore beta
    this->setAlpha(this->alpha - dalpha);

    double error = std::abs(sensitivity / variation - 1);
    return error < precision;
}

bool Test_Sabr::test_dRfdBeta(double money, double bumpRatio, double precision){
    double strike = this->forward + money;
    double sensitivity = this->dRfdBeta(strike);

    double Rf1 = this->Rf(strike);

    double dbeta = beta * bumpRatio;
    this->setBeta(this->beta + dbeta);
    double Rf2 = this->Rf(strike);

    double variation = (Rf2 - Rf1) / dbeta;

    // Restore beta
    this->setBeta(this->beta - dbeta);

    double error = std::abs(sensitivity / variation - 1);
    return error < precision;
}

bool Test_Sabr::test_dRfdRho(double money, double bumpRatio, double precision){
    double strike = this->forward + money;
    double sensitivity = this->dRfdRho(strike);

    double Rf1 = this->Rf(strike);

    double drho = rho * bumpRatio;
    this->setRho(this->rho + drho);
    double Rf2 = this->Rf(strike);

    double variation = (Rf2 - Rf1) / drho;

    // Restore beta
    this->setRho(this->rho - drho);

    double error = std::abs(sensitivity / variation - 1);
    return error < precision;
}

bool Test_Sabr::test_dRfdVovol(double money, double bumpRatio, double precision){
    double strike = this->forward + money;
    double sensitivity = this->dRfdVovol(strike);

    double Rf1 = this->Rf(strike);

    double dvovol = vovol * bumpRatio;
    this->setVovol(this->vovol + dvovol);
    double Rf2 = this->Rf(strike);

    double variation = (Rf2 - Rf1) / dvovol;

    // Restore beta
    this->setVovol(this->vovol - dvovol);

    double error = std::abs(sensitivity / variation - 1);
    return error < precision;
}

bool Test_Sabr::test_dRfdForward(double money, double bumpRatio, double precision){
    double strike = this->forward + money;
    double sensitivity = this->dRfdForward(strike);

    double Rf1 = this->Rf(strike);

    double dForward = forward * bumpRatio;
    this->setForward(this->forward + dForward);
    double Rf2 = this->Rf(strike);

    double variation = (Rf2 - Rf1) / dForward;

    // Restore beta
    this->setForward(this->forward - dForward);

    double error = std::abs(sensitivity / variation - 1);
    return error < precision;
}

// Tests volatility
bool Test_Sabr::test_dSdAlpha(double money, double bumpRatio, double precision){
    double strike = this->forward + money;
    double sensitivity = this->dSdAlpha(strike);

    double S1 = this->volatility(strike);

    double dalpha = alpha * bumpRatio;
    this->setAlpha(this->alpha + dalpha);
    double S2 = this->volatility(strike);

    double variation = (S2 - S1) / dalpha;

    // Restore beta
    this->setAlpha(this->alpha - dalpha);

    double error = std::abs(sensitivity / variation - 1);
    return error < precision;
}

bool Test_Sabr::test_dSdBeta(double money, double bumpRatio, double precision){
    double strike = this->forward + money;
    double sensitivity = this->dSdBeta(strike);

    double S1 = this->volatility(strike);

    double dbeta = beta * bumpRatio;
    this->setBeta(this->beta + dbeta);
    double S2 = this->volatility(strike);

    double variation = (S2 - S1) / dbeta;

    // Restore beta
    this->setBeta(this->beta - dbeta);

    double error = std::abs(sensitivity / variation - 1);
    return error < precision;
}

bool Test_Sabr::test_dSdRho(double money, double bumpRatio, double precision){
    double strike = this->forward + money;
    double sensitivity = this->dSdRho(strike);

    double S1 = this->volatility(strike);

    double drho = rho * bumpRatio;
    this->setRho(this->rho + drho);
    double S2 = this->volatility(strike);

    double variation = (S2 - S1) / drho;

    // Restore beta
    this->setRho(this->rho - drho);

    double error = std::abs(sensitivity / variation - 1);
    return error < precision;
}

bool Test_Sabr::test_dSdVovol(double money, double bumpRatio, double precision){
    double strike = this->forward + money;
    double sensitivity = this->dSdVovol(strike);

    double S1 = this->volatility(strike);

    double dvovol = vovol * bumpRatio;
    this->setVovol(this->vovol + dvovol);
    double S2 = this->volatility(strike);

    double variation = (S2 - S1) / dvovol;

    // Restore beta
    this->setVovol(this->vovol - dvovol);

    double error = std::abs(sensitivity / variation - 1);
    return error < precision;
}

bool Test_Sabr::test_dSdForward(double money, double bumpRatio, double precision){
    double strike = this->forward + money;
    double sensitivity = this->dSdForward(strike);

    double S1 = this->volatility(strike);

    double dforward = this->forward * bumpRatio;
    this->setForward(this->forward + dforward);
    double S2 = this->volatility(strike);

    double variation = (S2 - S1) / dforward;

    // Restore beta
    this->setForward(this->forward - dforward);

    double error = std::abs(sensitivity / variation - 1);
    return error < precision;
}
