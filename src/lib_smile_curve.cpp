#include "lib_smile_curve.h"
#include <math.h>

using namespace std;

// Approximation of the Gaussian function and its first order derivative
double Nprime(double x){
    double value = 1 / sqrt(2 * M_PI) * exp(-pow(x, 2) / 2);
    return value;
}

double N(double x){
    double y = 1 / (1 + 0.2316419 * fabs(x));
    double z = 1.330274429 * pow(y, 5) - 1.821255978 * pow(y, 4) + 1.781477937 * pow(y, 3) - 0.356563782 * pow(y, 2) + 0.31938153 * y;
    double n;

    if(x >= 0){
        n = 1 - z * exp(-0.5 * x * x) / sqrt(2 * M_PI);
    }
    else{
        n = z * exp(-0.5 * x * x) / sqrt(2 * M_PI);
    }

    return n;
}

// Constructor & destructor
Sabr::Sabr(){}

Sabr::Sabr(double a, double b, double r, double v, double f, int t0, int t){
    alpha = a;
    beta = b;
    rho = r;
    vovol = v;
    forward = f;
    calDate = t0;
    maturity = t;
}

Sabr::~Sabr(){}

// Getters
int Sabr::getCalDate() const{
    return this->calDate;
}

double Sabr::getAlpha() const{
    return this->alpha;
}

double Sabr::getBeta() const{
    return this->beta;
}

double Sabr::getRho() const{
    return this->rho;
}

double Sabr::getVovol() const{
    return this->vovol;
}

double Sabr::getForward() const{
    return this->forward;
}

double Sabr::getMaturity() const{
    return this->maturity;
}

// Setters
void Sabr::setCalDate(int date){
    this->calDate = date;
}

void Sabr::setAlpha(double a){
    this->alpha = a;
}

void Sabr::setBeta(double b){
    this->beta = b;
}

void Sabr::setRho(double r){
    this->rho = r;
}

void Sabr::setVovol(double v){
    this->vovol = v;
}

void Sabr::setForward(double f){
    this->forward = f;
}

void Sabr::setMaturity(int t){
    this->maturity = t;
}

// Main functions
double Sabr::volatility(double strike){
    double dcf = static_cast<double>(maturity - calDate) / 365.0;
    double W = this->W(strike);
    double zxz = this->zxz(strike);
    double Rf = this->Rf(strike);
    
    double sigma = alpha / W * zxz * (1 + dcf * Rf);
    return sigma;
}

double Sabr::dSdAlpha(double strike){
    double dcf = static_cast<double>(maturity - calDate) / 365.0;
    double W = this->W(strike);
    double zxz = this->zxz(strike);
    double dzxzdAlpha = this->dzxzdAlpha(strike);
    double dRfdAlpha = this->dRfdAlpha(strike);
    double Rf = this->Rf(strike);

    double dSdAlpha = alpha / W * (dzxzdAlpha * (1 + dcf * Rf) + zxz * dcf * dRfdAlpha) + 1 / W * zxz * (1 + dcf * Rf);
    return dSdAlpha;
}

double Sabr::dSdBeta(double strike){
    double dcf = static_cast<double>(maturity - calDate) / 365.0;
    double W = this->W(strike);
    double dWdBeta = this->dWdBeta(strike);
    double zxz = this->zxz(strike);
    double dzxzdBeta = this->dzxzdBeta(strike);
    double Rf = this->Rf(strike);
    double dRfdBeta = this->dRfdBeta(strike);

    double dSdBeta = - alpha * dWdBeta / (W * W) * zxz * (1 + dcf * Rf) + alpha / W * dzxzdBeta * (1 + dcf * Rf) + alpha / W * zxz * dcf * dRfdBeta;
    return dSdBeta;
}

double Sabr::dSdRho(double strike){
    double dcf = static_cast<double>(maturity - calDate) / 365.0;
    double W = this->W(strike);
    double zxz = this->zxz(strike);
    double dzxzdRho = this->dzxzdRho(strike);
    double Rf = this->Rf(strike);
    double dRfdRho = this->dRfdRho(strike);

    double dSdRho = alpha / W * (dzxzdRho * (1 + dcf * Rf) + zxz * dcf * dRfdRho);
    return dSdRho;
}

double Sabr::dSdVovol(double strike){
    double dcf = static_cast<double>(maturity - calDate) / 365.0;
    double W = this->W(strike);
    double zxz = this->zxz(strike);
    double dzxzdVovol = this->dzxzdVovol(strike);
    double Rf = this->Rf(strike);
    double dRfdVovol = this->dRfdVovol(strike);

    double dSdVovol = alpha / W * (dzxzdVovol * (1 + dcf * Rf) + zxz * dcf * dRfdVovol);
    return dSdVovol;
}

double Sabr::dSdForward(double strike){
    double dcf = static_cast<double>(maturity - calDate) / 365.0;
    double W = this->W(strike);
    double dWdForward = this->dWdForward(strike);
    double zxz = this->zxz(strike);
    double dzxzdForward = this->dzxzdForward(strike);
    double Rf = this->Rf(strike);
    double dRfdForward = this->dRfdForward(strike);

    double dSdForward = -alpha * dWdForward / (W * W) * zxz * (1 + dcf * Rf);
    dSdForward += alpha / W * dzxzdForward * (1 + dcf * Rf);
    dSdForward += alpha / W * zxz * dcf * dRfdForward;

    return dSdForward;
}

// Alpha
double Sabr::d2SdAlphadAlpha(double strike){
    double dcf = static_cast<double>(maturity - calDate) / 365.0;
    double W = this->W(strike);

    double zxz = this->zxz(strike);
    double dzxzdAlpha = this->dzxzdAlpha(strike);
    double d2zxzdAlphadAlpha = this->d2zxzdAlphadAlpha(strike);

    double Rf = this->Rf(strike);
    double dRfdAlpha = this->dRfdAlpha(strike);
    double d2RfdAlphadAlpha = this->d2RfdAlphadAlpha(strike);

    double d2SdAlphadAlpha = 1 / W * (dzxzdAlpha * (1 + dcf * Rf) + zxz * dcf * dRfdAlpha);
    d2SdAlphadAlpha += alpha / W * (d2zxzdAlphadAlpha * (1 + dcf * Rf) + dzxzdAlpha * (1 + dcf * dRfdAlpha) + dzxzdAlpha * dcf * dRfdAlpha + zxz * dcf * d2RfdAlphadAlpha);
    d2SdAlphadAlpha += 1 / W * dzxzdAlpha * (1 + dcf * Rf) + 1 / W * zxz * dcf * dRfdAlpha;
    return d2SdAlphadAlpha;
}

double Sabr::d2SdBetadAlpha(double strike){
    double dcf = static_cast<double>(maturity - calDate) / 365.0;
    double W = this->W(strike); // does not depend on alpha
    double dWdBeta = this->dWdBeta(strike); // does not depend on alpha

    double zxz = this->zxz(strike); // depend on alpha
    double dzxzdAlpha = this->dzxzdAlpha(strike);
    double dzxzdBeta = this->dzxzdBeta(strike); // depend on alpha
    double d2zxzdBetadAlpha = this->d2zxzdBetadAlpha(strike);

    double Rf = this->Rf(strike); // depend on alpha
    double dRfdAlpha = this->dRfdAlpha(strike);
    double dRfdBeta = this->dRfdBeta(strike);
    double d2RfdBetadAlpha = this->d2RfdBetadAlpha(strike);

    double d2SdBetadAlpha = - 1 * dWdBeta / (W * W) * zxz * (1 + dcf * Rf) - alpha * dWdBeta / (W * W) * dzxzdAlpha * (1 + dcf * Rf) - alpha * dWdBeta / (W * W) * zxz * dcf * dRfdAlpha;
    d2SdBetadAlpha += 1 / W * dzxzdBeta * (1 + dcf * Rf) + alpha / W * d2zxzdBetadAlpha * (1 + dcf * Rf) + alpha / W * dzxzdBeta * dcf * dRfdAlpha;
    d2SdBetadAlpha += 1 / W * zxz * dcf * dRfdBeta + alpha / W * dzxzdAlpha * dcf * dRfdBeta + alpha / W * zxz * dcf * d2RfdBetadAlpha;
    return d2SdBetadAlpha;
}

double Sabr::d2SdRhodAlpha(double strike){
    double dcf = static_cast<double>(maturity - calDate) / 365.0;
    double W = this->W(strike);

    double zxz = this->zxz(strike);
    double dzxzdRho = this->dzxzdRho(strike);
    double dzxzdAlpha = this->dzxzdAlpha(strike);
    double d2zxzdRhodAlpha = this->d2zxzdRhodAlpha(strike);
    double Rf = this->Rf(strike);
    double dRfdRho = this->dRfdRho(strike);
    double dRfdAlpha = this->dRfdAlpha(strike);
    double d2RfdRhodAlpha = this->d2RfdRhodAlpha(strike);

    double d2SdRhodAlpha = 1 / W * (dzxzdRho * (1 + dcf * Rf) + zxz * dcf * dRfdRho);
    d2SdRhodAlpha += alpha / W * (d2zxzdRhodAlpha * (1 + dcf * Rf) + dzxzdRho * dcf * dRfdAlpha + dzxzdAlpha * dcf * dRfdRho + zxz * dcf * d2RfdRhodAlpha);
    return d2SdRhodAlpha;
}

double Sabr::d2SdVovoldAlpha(double strike){
    double dcf = static_cast<double>(maturity - calDate) / 365.0;
    double W = this->W(strike);

    double zxz = this->zxz(strike);
    double dzxzdVovol = this->dzxzdVovol(strike);
    double dzxzdAlpha = this->dzxzdAlpha(strike);
    double d2zxzdVovoldAlpha = this->d2zxzdVovoldAlpha(strike);
    double Rf = this->Rf(strike);
    double dRfdVovol = this->dRfdVovol(strike);
    double dRfdAlpha = this->dRfdAlpha(strike);
    double d2RfdVovoldAlpha = this->d2RfdVovoldAlpha(strike);

    double d2SdVovoldAlpha = 1 / W * (dzxzdVovol * (1 + dcf * Rf) + zxz * dcf * dRfdVovol);
    d2SdVovoldAlpha += alpha / W * (d2zxzdVovoldAlpha * (1 + dcf * Rf) + dzxzdVovol * (1 + dcf * dRfdAlpha) + dzxzdAlpha * dcf * dRfdVovol + zxz * dcf * d2RfdVovoldAlpha);
    return d2SdVovoldAlpha;
}

double Sabr::d2SdForwarddAlpha(double strike){
    double dcf = static_cast<double>(maturity - calDate) / 365.0;
    double W = this->W(strike);
    double dWdForward = this->dWdForward(strike);

    double zxz = this->zxz(strike);
    double dzxzdAlpha = this->dzxzdAlpha(strike);
    double dzxzdForward = this->dzxzdForward(strike);
    double d2zxzdForwarddAlpha = this->d2zxzdForwarddAlpha(strike);
    double Rf = this->Rf(strike);
    double dRfdAlpha = this->dRfdAlpha(strike);
    double dRfdForward = this->dRfdForward(strike);
    double d2RfdForwarddAlpha = this->d2RfdForwarddAlpha(strike);

    double d2SdForwarddAlpha = -dWdForward / W / W * zxz * (1 + dcf * Rf);
    d2SdForwarddAlpha += -alpha * dWdForward / W / W * dzxzdAlpha * (1 + dcf * Rf);
    d2SdForwarddAlpha += -alpha * dWdForward / W / W * zxz * (dcf * dRfdAlpha);

    d2SdForwarddAlpha += 1 / W * dzxzdForward * (1 + dcf * Rf);
    d2SdForwarddAlpha += alpha / W * d2zxzdForwarddAlpha * (1 + dcf * Rf);
    d2SdForwarddAlpha += alpha / W * dzxzdForward * (dcf * dRfdAlpha);

    d2SdForwarddAlpha += 1 / W * zxz * dcf * dRfdForward;
    d2SdForwarddAlpha += alpha / W * dzxzdAlpha * dcf * dRfdForward;
    d2SdForwarddAlpha += alpha / W * zxz * dcf * d2RfdForwarddAlpha;

    return d2SdForwarddAlpha;
}

// Beta
double Sabr::d2SdAlphadBeta(double strike){
    return this->d2SdBetadAlpha(strike);
}

double Sabr::d2SdBetadBeta(double strike){
    double dcf = static_cast<double>(maturity - calDate) / 365.0;
    double W = this->W(strike);
    double dWdBeta = this->dWdBeta(strike);
    double d2WdBetadBeta = this->d2WdBetadBeta(strike);
    double zxz = this->zxz(strike);
    double dzxzdBeta = this->dzxzdBeta(strike);
    double d2zxzdBetadBeta = this->d2zxzdBetadBeta(strike);
    double Rf = this->Rf(strike);
    double dRfdBeta = this->dRfdBeta(strike);
    double d2RfdBetadBeta = this->d2RfdBetadBeta(strike);

    double d2SdBetadBeta = - alpha * (d2WdBetadBeta * W * W - dWdBeta * 2 * dWdBeta * W) / (W * W * W * W) * zxz * (1 + dcf * Rf) - alpha * dWdBeta / (W * W) * dzxzdBeta * (1 + dcf * Rf) - alpha * dWdBeta / (W * W) * zxz * dcf * dRfdBeta;
    d2SdBetadBeta += -alpha * dWdBeta / W / W * dzxzdBeta * (1 + dcf * Rf) + alpha / W * d2zxzdBetadBeta * (1 + dcf * Rf) + alpha / W * dzxzdBeta * dcf * dRfdBeta;
    d2SdBetadBeta += -alpha * dWdBeta / W / W * zxz * dcf * dRfdBeta + alpha / W * dzxzdBeta * dcf * dRfdBeta + alpha / W * zxz * dcf * d2RfdBetadBeta;
    return d2SdBetadBeta;
}

double Sabr::d2SdRhodBeta(double strike){
    double dcf = static_cast<double>(maturity - calDate) / 365.0;
    double W = this->W(strike);
    double dWdBeta = this->dWdBeta(strike);
    double zxz = this->zxz(strike);
    double dzxzdBeta = this->dzxzdBeta(strike);
    double dzxzdRho = this->dzxzdRho(strike);
    double d2zxzdRhodBeta = this->d2zxzdRhodBeta(strike);
    double Rf = this->Rf(strike);
    double dRfdBeta = this->dRfdBeta(strike);
    double dRfdRho = this->dRfdRho(strike);
    double d2RfdRhodBeta = this->d2RfdRhodBeta(strike);

    double d2SdRhodBeta = -alpha * dWdBeta / (W * W) * (dzxzdRho * (1 + dcf * Rf) + zxz * dcf * dRfdRho);
    d2SdRhodBeta += alpha / W * (d2zxzdRhodBeta * (1 + dcf * Rf) + dzxzdRho * (dcf * dRfdBeta) + dzxzdBeta * dcf * dRfdRho + zxz * dcf * d2RfdRhodBeta);
    return d2SdRhodBeta;
}

double Sabr::d2SdVovoldBeta(double strike){
    double dcf = static_cast<double>(maturity - calDate) / 365.0;
    double W = this->W(strike);
    double dWdBeta = this->dWdBeta(strike);
    double zxz = this->zxz(strike);
    double dzxzdBeta = this->dzxzdBeta(strike);
    double dzxzdVovol = this->dzxzdVovol(strike);
    double d2zxzdVovoldBeta = this->d2zxzdVovoldBeta(strike);
    double Rf = this->Rf(strike);
    double dRfdBeta = this->dRfdBeta(strike);
    double dRfdVovol = this->dRfdVovol(strike);
    double d2RfdVovoldBeta = this->d2RfdVovoldBeta(strike);

    double d2SdVovoldBeta = -alpha * dWdBeta / (W * W) * (dzxzdVovol * (1 + dcf * Rf) + zxz * dcf * dRfdVovol);
    d2SdVovoldBeta += alpha / W * (d2zxzdVovoldBeta * (1 + dcf * Rf) + dzxzdVovol * (dcf * dRfdBeta) + dzxzdBeta * dcf * dRfdVovol + zxz * dcf * d2RfdVovoldBeta);
    return d2SdVovoldBeta;
}

double Sabr::d2SdForwarddBeta(double strike){
    double dcf = static_cast<double>(maturity - calDate) / 365.0;
    double W = this->W(strike);
    double dWdBeta = this->dWdBeta(strike);
    double dWdForward = this->dWdForward(strike);
    double d2WdForwarddBeta = this->d2WdForwarddBeta(strike);
    double zxz = this->zxz(strike);
    double dzxzdBeta = this->dzxzdBeta(strike);
    double dzxzdForward = this->dzxzdForward(strike);
    double d2zxzdForwarddBeta = this->d2zxzdForwarddBeta(strike);
    double Rf = this->Rf(strike);
    double dRfdBeta = this->dRfdBeta(strike);
    double dRfdForward = this->dRfdForward(strike);
    double d2RfdForwarddBeta = this->d2RfdForwarddBeta(strike);

    double dSdForward = -alpha * (d2WdForwarddBeta * W * W - dWdForward * 2 * dWdBeta * W) / (W * W * W * W) * zxz * (1 + dcf * Rf);
    dSdForward -= alpha * dWdForward / (W * W) * dzxzdBeta * (1 + dcf * Rf);
    dSdForward -= alpha * dWdForward / (W * W) * zxz * (dcf * dRfdBeta);
    dSdForward -= alpha * dWdBeta / (W * W) * dzxzdForward * (1 + dcf * Rf);
    dSdForward += alpha / W * d2zxzdForwarddBeta * (1 + dcf * Rf);
    dSdForward += alpha / W * dzxzdForward * (dcf * dRfdBeta);
    dSdForward -= alpha * dWdBeta / (W * W) * zxz * dcf * dRfdForward;
    dSdForward += alpha / W * dzxzdBeta * dcf * dRfdForward;
    dSdForward += alpha / W * zxz * dcf * d2RfdForwarddBeta;

    return dSdForward;
}

// Rho
double Sabr::d2SdAlphadRho(double strike){
    return this->d2SdRhodAlpha(strike);
}

double Sabr::d2SdBetadRho(double strike){
    return this->d2SdRhodBeta(strike);
}

double Sabr::d2SdRhodRho(double strike){
    double dcf = static_cast<double>(maturity - calDate) / 365.0;
    double W = this->W(strike);
    double zxz = this->zxz(strike);
    double dzxzdRho = this->dzxzdRho(strike);
    double d2zxzdRhodRho = this->d2zxzdRhodRho(strike);
    double Rf = this->Rf(strike);
    double dRfdRho = this->dRfdRho(strike);
    double d2RfdRhodRho = this->d2RfdRhodRho(strike);

    double d2SdRhodRho = alpha / W * (d2zxzdRhodRho * (1 + dcf * Rf) + dzxzdRho * (dcf * dRfdRho) + dzxzdRho * dcf * dRfdRho + zxz * dcf * d2RfdRhodRho);
    return d2SdRhodRho;
}

double Sabr::d2SdVovoldRho(double strike){
    double dcf = static_cast<double>(maturity - calDate) / 365.0;
    double W = this->W(strike);
    double zxz = this->zxz(strike);
    double dzxzdRho = this->dzxzdRho(strike);
    double dzxzdVovol = this->dzxzdVovol(strike);
    double d2zxzdVovoldRho = this->dzxzdVovol(strike);
    double Rf = this->Rf(strike);
    double dRfdRho = this->dRfdRho(strike);
    double dRfdVovol = this->dRfdVovol(strike);
    double d2RfdVovoldRho = this->d2RfdVovoldRho(strike);

    double d2SdVovoldRho = alpha / W * (d2zxzdVovoldRho * (1 + dcf * Rf) + dzxzdVovol * dcf * dRfdRho + dzxzdRho * dcf * dRfdVovol + zxz * dcf * d2RfdVovoldRho);
    return d2SdVovoldRho;
}

double Sabr::d2SdForwarddRho(double strike){
    double dcf = static_cast<double>(maturity - calDate) / 365.0;
    double W = this->W(strike);
    double dWdForward = this->dWdForward(strike);
    double zxz = this->zxz(strike);
    double dzxzdRho = this->dzxzdRho(strike);
    double dzxzdForward = this->dzxzdForward(strike);
    double d2zxzdForwarddRho = this->d2zxzdForwarddRho(strike);
    double Rf = this->Rf(strike);
    double dRfdRho = this->dRfdRho(strike);
    double dRfdForward = this->dRfdForward(strike);
    double d2RfdForwarddRho = this->d2RfdForwarddRho(strike);

    double d2SdForwarddRho = -alpha * dWdForward / (W * W) * dzxzdRho * (1 + dcf * Rf);
    d2SdForwarddRho += -alpha * dWdForward / (W * W) * zxz * dcf * dRfdRho;
    d2SdForwarddRho += alpha / W * d2zxzdForwarddRho * (1 + dcf * Rf);
    d2SdForwarddRho += alpha / W * dzxzdForward * dcf * dRfdRho;
    d2SdForwarddRho += alpha / W * dzxzdRho * dcf * dRfdForward;
    d2SdForwarddRho += alpha / W * zxz * dcf * d2RfdForwarddRho;

    return d2SdForwarddRho;
}

// Vovol
double Sabr::d2SdAlphadVovol(double strike){
    return this->d2SdVovoldAlpha(strike);
}

double Sabr::d2SdBetadVovol(double strike){
    return this->d2SdVovoldBeta(strike);
}

double Sabr::d2SdRhodVovol(double strike){
    return this->d2SdVovoldRho(strike);
}

double Sabr::d2SdVovoldVovol(double strike){
    double dcf = static_cast<double>(maturity - calDate) / 365.0;
    double W = this->W(strike);
    double zxz = this->zxz(strike);
    double dzxzdVovol = this->dzxzdVovol(strike);
    double d2zxzdVovoldVovol = this->d2zxzdVovoldVovol(strike);
    double Rf = this->Rf(strike);
    double dRfdVovol = this->dRfdVovol(strike);
    double d2RfdVovoldVovol = this->d2RfdVovoldVovol(strike);

    double dSdVovol = alpha / W * (d2zxzdVovoldVovol * (1 + dcf * Rf) + dzxzdVovol * (dcf * dRfdVovol) + dzxzdVovol * dcf * dRfdVovol + zxz * dcf * d2RfdVovoldVovol);
    return dSdVovol;
}

double Sabr::d2SdForwarddVovol(double strike){
    double dcf = static_cast<double>(maturity - calDate) / 365.0;
    double W = this->W(strike);
    double dWdForward = this->dWdForward(strike);
    double zxz = this->zxz(strike);
    double dzxzdVovol = this->dzxzdVovol(strike);
    double dzxzdForward = this->dzxzdForward(strike);
    double d2zxzdForwarddVovol = this->d2zxzdForwarddVovol(strike);
    double Rf = this->Rf(strike);
    double dRfdVovol = this->dRfdVovol(strike);
    double dRfdForward = this->dRfdForward(strike);
    double d2RfdForwarddVovol = this->d2RfdForwarddVovol(strike);

    double d2SdForwarddVovol = -alpha * dWdForward / (W * W) * dzxzdVovol * (1 + dcf * Rf);
    d2SdForwarddVovol += -alpha * dWdForward / (W * W) * zxz * (dcf * dRfdVovol);
    d2SdForwarddVovol += alpha / W * d2zxzdForwarddVovol * (1 + dcf * Rf);
    d2SdForwarddVovol += alpha / W * dzxzdForward * (dcf * dRfdVovol);
    d2SdForwarddVovol += alpha / W * dzxzdVovol * dcf * dRfdForward;
    d2SdForwarddVovol += alpha / W * zxz * dcf * d2RfdForwarddVovol;

    return d2SdForwarddVovol;
}

// Forward
double Sabr::d2SdAlphadForward(double strike){
    return this->d2SdForwarddAlpha(strike);
}

double Sabr::d2SdBetadForward(double strike){
    return this->d2SdForwarddBeta(strike);
}

double Sabr::d2SdRhodForward(double strike){
    return this->d2SdForwarddRho(strike);
}

double Sabr::d2SdVovoldForward(double strike){
    return this->d2SdForwarddVovol(strike);
}

double Sabr::d2SdForwarddForward(double strike){
    double dcf = static_cast<double>(maturity - calDate) / 365.0;
    double W = this->W(strike);
    double dWdForward = this->dWdForward(strike);
    double d2WdForwarddForward = this->d2WdForwarddForward(strike);
    double zxz = this->zxz(strike);
    double dzxzdForward = this->dzxzdForward(strike);
    double d2zxzdForwarddForward = this->d2zxzdForwarddForward(strike);
    double Rf = this->Rf(strike);
    double dRfdForward = this->dRfdForward(strike);
    double d2RfdForwarddForward = this->d2RfdForwarddForward(strike);

    double dSdForward = -alpha * (d2WdForwarddForward * W * W - dWdForward * 2 * dWdForward * W) / (W * W * W * W) * zxz * (1 + dcf * Rf);
    dSdForward += -alpha * dWdForward / (W * W) * dzxzdForward * (1 + dcf * Rf);
    dSdForward += -alpha * dWdForward / (W * W) * zxz * (dcf * dRfdForward);
    dSdForward -= alpha / (W * W) * dzxzdForward * (1 + dcf * Rf);
    dSdForward += alpha / W * d2zxzdForwarddForward * (1 + dcf * Rf);
    dSdForward += alpha / W * dzxzdForward * (dcf * dRfdForward);
    dSdForward -= alpha * dWdForward / (W * W) * zxz * dcf * dRfdForward;
    dSdForward += alpha / W * dzxzdForward * dcf * dRfdForward;
    dSdForward += alpha / W * zxz * dcf * d2RfdForwarddForward;

    return dSdForward;
}

// Private intermediates
double Sabr::W(double strike){
    double W = pow(forward * strike, (1 - beta) / 2) * (1 + this->Lf(strike));
    return W;
}

double Sabr::dWdBeta(double strike){
    double dWdBeta = pow(forward * strike, (1 - beta) / 2) * this->dLfdBeta(strike) + -0.5 * log(forward * strike) * pow(forward * strike, (1 - beta) / 2) * (1 + this->Lf(strike));
    return dWdBeta;
}

double Sabr::d2WdBetadBeta(double strike){
    double dLfdBeta = this->dLfdBeta(strike);
    double d2WdBetadBeta = -0.5 * log(forward * strike) * pow(forward * strike, (1 - beta) / 2) * dLfdBeta;
    d2WdBetadBeta += pow(forward * strike, (1 - beta) / 2) * this->d2LfdBetadBeta(strike);
    d2WdBetadBeta += 0.25 * log(forward * strike) * log(forward * strike) * pow(forward * strike, (1 - beta) / 2) * (1 + this->Lf(strike));
    d2WdBetadBeta += -0.5 * log(forward * strike) * pow(forward * strike, (1 - beta) / 2) * dLfdBeta;
    return d2WdBetadBeta;
}

double Sabr::dWdForward(double strike){
    double dWdForward = (1 - beta) / 2 * pow(forward * strike, (-1 - beta) / 2) * strike * (1 + this->Lf(strike));
    dWdForward += pow(forward * strike, (1 - beta) / 2) * this->dLfdForward(strike);
    return dWdForward;
}

double Sabr::d2WdForwarddBeta(double strike){
    double Lf = this->Lf(strike);
    double d2WdForwarddBeta = -1 / 2 * pow(forward * strike, (-1 - beta) / 2) * strike * (1 + Lf);
    d2WdForwarddBeta += (1 - beta) / 2 * -0.5 * log(forward * strike) * pow(forward * strike, (-1 - beta) / 2) * strike * (1 + Lf);
    d2WdForwarddBeta += (1 - beta) / 2 * pow(forward * strike, (-1 - beta) / 2) * strike * this->dLfdBeta(strike);
    d2WdForwarddBeta += -0.5 * log(forward * strike) * pow(forward * strike, (1 - beta) / 2) * this->dLfdForward(strike);
    d2WdForwarddBeta += pow(forward * strike, (1 - beta) / 2) * this->d2LfdForwarddBeta(strike);
    return d2WdForwarddBeta;
}

double Sabr::d2WdForwarddForward(double strike){
    double dLfdForward = this->dLfdForward(strike);
    double dWdForward = (1 - beta) / 2 * strike * (-1 - beta) / 2 * pow(forward * strike, (-3 - beta) / 2) * strike * (1 + this->Lf(strike));
    dWdForward += (1 - beta) / 2 * pow(forward * strike, (-1 - beta) / 2) * strike * dLfdForward;
    dWdForward += strike * (1 - beta) / 2 * pow(forward * strike, (-1 - beta) / 2) * dLfdForward;
    dWdForward += pow(forward * strike, (1 - beta) / 2) * this->d2LfdForwarddForward(strike);
    return dWdForward;
}

double Sabr::Lf(double strike){
    double ratio = forward / strike;
    double Lf = (1 - beta)*(1 - beta) / 24 * log(ratio)*log(ratio) + pow(1 - beta, 4) / 1920 * pow(log(ratio), 4);
    return Lf;
}

double Sabr::dLfdBeta(double strike){
    double ratio = forward / strike;
    double dLfdBeta = (beta - 1) / 12 * log(ratio)*log(ratio) - pow(1 - beta, 3) / 480 * pow(log(ratio), 4);
    return dLfdBeta;
}

double Sabr::d2LfdBetadBeta(double strike){
    double ratio = forward / strike;
    double dLfdBeta = 1 / 12 * log(ratio)*log(ratio) + 3 * pow(1 - beta, 2) / 480 * pow(log(ratio), 4);
    return dLfdBeta;
}

double Sabr::dLfdForward(double strike){
    double ratio = forward / strike;
    double dratio = 1 / strike;
    double dLfdForward = (1 - beta)*(1 - beta) / 24 * ( 2 * log(ratio) * dratio / ratio);
    dLfdForward += pow(1 - beta, 4) / 1920 *  4 * pow(log(ratio), 3) * dratio / ratio;
    return dLfdForward;
}

double Sabr::d2LfdForwarddBeta(double strike){
    double ratio = forward / strike;
    double dratio = 1 / strike;
    double dLfdForward = -2 * (1 - beta) / 24 * ( 2 * log(ratio) * dratio / ratio);
    dLfdForward -= 4 * pow(1 - beta, 3) / 1920 *  4 * pow(log(ratio), 3) * dratio / ratio;
    return dLfdForward;
}

double Sabr::d2LfdForwarddForward(double strike){
    double ratio = forward / strike;
    double dratio = 1 / strike;
    double d2LfdForwarddForward = (1 - beta)*(1 - beta) / 24 * ( 2 * dratio / ratio * dratio / ratio - 2 * log(ratio) * dratio * dratio / ratio / ratio);
    d2LfdForwarddForward += pow(1 - beta, 4) / 1920 *  4 * dratio / ratio * 3 * pow(log(ratio), 2) * dratio / ratio;
    d2LfdForwarddForward -= pow(1 - beta, 4) / 1920 *  4 * pow(log(ratio), 3) * dratio * dratio / ratio / ratio;
    return d2LfdForwarddForward;
}

double Sabr::z(double strike){
    double z = vovol / alpha * pow(forward * strike, (1 - beta) / 2) * log(forward / strike);
    return z;
}

double Sabr::dzdAlpha(double strike){
    double dzdAlpha =  - 1 / alpha * this->z(strike);
    return dzdAlpha;
}

double Sabr::d2zdAlphadAlpha(double strike){
    double d2zdAlphadAlpha =  - 1 / (alpha * alpha) * (this->dzdAlpha(strike) * alpha - this->z(strike));
    return d2zdAlphadAlpha;
}

double Sabr::d2zdAlphadBeta(double strike){
    double d2zdAlphadBeta = this->d2zdBetadAlpha(strike);
    return d2zdAlphadBeta;
}

double Sabr::d2zdAlphadRho(double strike){
    double d2zdAlphadRho =  - 1 / alpha * this->dzdRho(strike);
    return d2zdAlphadRho;
}

double Sabr::dzdBeta(double strike){
    double dzdBeta = this->z(strike) * log(forward * strike) * (-0.5);
    return dzdBeta;
}

double Sabr::d2zdBetadAlpha(double strike){
    double d2zdBetadAlpha = this->dzdAlpha(strike) * log(forward * strike) * (-0.5);
    return d2zdBetadAlpha;
}

double Sabr::d2zdBetadBeta(double strike){
    double dzdBeta = this->dzdBeta(strike) * log(forward * strike) * (-0.5);
    return dzdBeta;
}

double Sabr::d2zdBetadRho(double strike){
    return this->d2zdRhodBeta(strike);
}

double Sabr::dzdRho(double strike){
    return 0;
}

double Sabr::d2zdRhodAlpha(double strike){
    return 0;
}

double Sabr::d2zdRhodBeta(double strike){
    return 0;
}

double Sabr::d2zdRhodRho(double strike){
    return 0;
}

double Sabr::dzdVovol(double strike){
    double dzdVovol = 1 / alpha * pow(forward * strike, (1 - beta) / 2) * log(forward / strike);
    return dzdVovol;
}

double Sabr::d2zdVovoldAlpha(double strike){
    double d2zdVovoldAlpha = -1 / alpha / alpha * pow(forward * strike, (1 - beta) / 2) * log(forward / strike);
    return d2zdVovoldAlpha;
}

double Sabr::d2zdVovoldBeta(double strike){
    double dzdVovol = -1 / alpha * 0.5 * log(forward * strike) * pow(forward * strike, (1 - beta) / 2) * log(forward / strike);
    return dzdVovol;
}

double Sabr::d2zdVovoldRho(double strike){
    return 0;
}

double Sabr::d2zdVovoldVovol(double strike){
    return 0;
}

double Sabr::dzdForward(double strike){
    double dzdForward = vovol / alpha * pow(forward * strike, (1 - beta) / 2) / forward;
    dzdForward += vovol / alpha * (1 - beta) / 2 * pow(forward * strike, (- 1 - beta) / 2) * strike * log(forward / strike);
    return dzdForward;
}

double Sabr::d2zdForwarddAlpha(double strike){
    double d2zdForwarddAlpha = -vovol / (alpha * alpha) * pow(forward * strike, (1 - beta) / 2) / forward;
    d2zdForwarddAlpha -= vovol / (alpha * alpha) * (1 - beta) / 2 * pow(forward * strike, (- 1 - beta) / 2) * strike * log(forward / strike);
    return d2zdForwarddAlpha;
}

double Sabr::d2zdForwarddBeta(double strike){
    double d2zdForwarddBeta = vovol / alpha * (-0.5) * log(forward * strike) * pow(forward * strike, (1 - beta) / 2) / forward;
    d2zdForwarddBeta += vovol / alpha * (-1) / 2 * pow(forward * strike, (- 1 - beta) / 2) * strike * log(forward / strike);
    d2zdForwarddBeta += vovol / alpha * (1 - beta) / 2 * (-0.5) * log(forward * strike) * pow(forward * strike, (- 1 - beta) / 2) * strike * log(forward / strike);
    return d2zdForwarddBeta;
}

double Sabr::d2zdForwarddRho(double strike){
    return 0;
}

double Sabr::d2zdForwarddVovol(double strike){
    double d2zdForwarddVovol = 1 / alpha * pow(forward * strike, (1 - beta) / 2) / forward;
    d2zdForwarddVovol += 1 / alpha * (1 - beta) / 2 * pow(forward * strike, (- 1 - beta) / 2) * strike * log(forward / strike);
    return d2zdForwarddVovol;
}

double Sabr::d2zdForwarddForward(double strike){
    double d2zdForwarddForward = vovol / alpha * (strike * (1 - beta) / 2 * pow(forward * strike, (-1 - beta) / 2) * forward - pow(forward * strike, (1 - beta) / 2)) / (forward * forward);
    d2zdForwarddForward += vovol / alpha * (1 - beta) / 2 * strike * (strike * (- 1 - beta) / 2 * pow(forward * strike, (- 3 - beta) / 2) * log(forward / strike) + pow(forward * strike, (- 1 - beta) / 2) * (1 / forward));
    return d2zdForwarddForward;
}

double Sabr::xz(double strike, double z){
    double xz = log((sqrt(1 - 2 * rho * z + z * z) + z - rho) / (1 - rho));
    return xz;
}

double Sabr::dxzdz(double strike, double z){
    double dxzdz = (1 + 0.5 * (2 * z - 2 * rho) / sqrt(1 - 2 * rho * z + z * z) ) / (sqrt(1 - 2 * rho * z + z * z) + z - rho);
    return dxzdz;
}

double Sabr::d2xzdzdz(double strike, double z){
    double d2xzdzdz = ( (sqrt(1 - 2 * rho * z + z * z) - 0.5 * (- 2 * rho +  2 * z) / sqrt(1 - 2 * rho * z + z * z) * (z - rho)) / (1 - 2 * rho * z + z * z) ) * (sqrt(1 - 2 * rho * z + z * z) + z - rho);
    d2xzdzdz -= (0.5 * (- 2 * rho +  2 * z) / sqrt(1 - 2 * rho * z + z * z) + 1) * (1 + (z - rho) / sqrt(1 - 2 * rho * z + z * z) );
    d2xzdzdz /= (sqrt(1 - 2 * rho * z + z * z) + z - rho) * (sqrt(1 - 2 * rho * z + z * z) + z - rho);
    return d2xzdzdz;
}

double Sabr::d2xzdzdRho(double strike, double z){
    double d2xzdzdRho = ((-sqrt(1 - 2 * rho * z + z * z) - (z - rho) * 0.5 * (-2 * z) / sqrt(1 - 2 * rho * z + z * z)) / (1 - 2 * rho * z + z * z) ) * (sqrt(1 - 2 * rho * z + z * z) + z - rho);
    d2xzdzdRho -= (1 + 0.5 * (2 * z - 2 * rho) / sqrt(1 - 2 * rho * z + z * z) ) * (0.5 * (-2 * z) / sqrt(1 - 2 * rho * z + z * z) - 1);
    d2xzdzdRho /= ((sqrt(1 - 2 * rho * z + z * z) + z - rho) * (sqrt(1 - 2 * rho * z + z * z) + z - rho));
    return d2xzdzdRho;
}

double Sabr::dxzdAlpha(double strike){
    double z = this->z(strike);
    double dxzdAlpha = this->dxzdz(strike, z) * this->dzdAlpha(strike);
    return dxzdAlpha;
}

double Sabr::d2xzdAlphadAlpha(double strike){
    double z = this->z(strike);
    double t1 = this->d2xzdzdz(strike, z) * this->dzdAlpha(strike) * this->dzdAlpha(strike);
    double t2 = this->dxzdz(strike, z) * this->d2zdAlphadAlpha(strike);
    // double t3 = this->d2xzdzdz(strike, z) * this->d2zdAlphadAlpha(strike);
    // double d2xzdAlphadAlpha = 0.5 * this->d2xzdzdz(strike, z) * (this->dzdAlpha(strike) * this->dzdAlpha(strike) + this->d2zdAlphadAlpha(strike));
    // d2xzdAlphadAlpha += this->dxzdz(strike, z) * this->d2zdAlphadAlpha(strike);
    return t1 + t2;
}

double Sabr::d2xzdAlphadRho(double strike){
    double z = this->z(strike);
    double d2xzdAlphadRho = this->dxzdz(strike, z) * this->d2zdAlphadRho(strike);
    d2xzdAlphadRho += this->d2xzdRhodz(strike, z) * this->dzdAlpha(strike);
    return d2xzdAlphadRho;
}

double Sabr::dxzdBeta(double strike){
    double z = this->z(strike);
    double dxzBbeta = this->dxzdz(strike, z) * this->dzdBeta(strike);
    return dxzBbeta;
}

double Sabr::d2xzdBetadAlpha(double strike){
    double z = this->z(strike);
    double d2xzdBetadAlpha = this->d2xzdzdz(strike, z) * this->dzdBeta(strike) * this->dzdAlpha(strike);
    d2xzdBetadAlpha += this->dxzdz(strike, z) * this->d2zdBetadAlpha(strike);
    return d2xzdBetadAlpha;
}

double Sabr::d2xzdBetadRho(double strike){
    double z = this->z(strike);
    double d2xzdBetadRho = this->d2xzdRhodz(strike, z) * this->dzdBeta(strike);
    d2xzdBetadRho += this->dxzdz(strike, z) * this->d2zdBetadRho(strike);
    return d2xzdBetadRho;
}

double Sabr::d2xzdAlphadBeta(double strike){
    return this->d2xzdBetadAlpha(strike);
}

double Sabr::d2xzdBetadBeta(double strike){
    double z = this->z(strike);
    double d2xzdBetadBeta = this->d2xzdzdz(strike, z) * this->dzdBeta(strike) * this->dzdBeta(strike);
    d2xzdBetadBeta += this->dxzdz(strike, z) * this->d2zdBetadBeta(strike);
    return d2xzdBetadBeta;
}

double Sabr::dxzdRho(double strike){
    double z = this->z(strike);
    double dxzdRho;

    // dxzdRho = (((-1 + -z / sqrt(1 - 2 * rho * z + z * z)) * (1 - rho) + (sqrt(1 - 2 * rho * z + z * z) + z - rho)) / ((1 - rho) * (1 - rho))) * (1 - rho) / (sqrt(1 - 2 * rho * z + z * z) + z - rho);

    // Rewriting
    double fRho = 1 - 2 * rho * z + z * z;
    double dfRhodRho = -2 * z;
    double racineFRho = sqrt(fRho);
    double dracineFRhodRho = dfRhodRho / fRho * 0.5 / racineFRho;
    dxzdRho = (-1 + -z / racineFRho) / (racineFRho + z - rho);
    dxzdRho += (racineFRho + z - rho) / ((1 - rho) * (racineFRho + z - rho));

    return dxzdRho;
}

double Sabr::d2xzdRhodz(double strike, double z){
    double d2xzdRhodz = d2xzdzdRho(strike, z);
    return d2xzdRhodz;
}

double Sabr::d2xzdRhodAlpha(double strike){
    double d2xzdRhodAlpha = this->d2xzdAlphadRho(strike);
    return d2xzdRhodAlpha;
}

double Sabr::d2xzdRhodBeta(double strike){
    return this->d2xzdBetadRho(strike);
}

double Sabr::d2xzdRhodRho(double strike){
    double z = this->z(strike);
    /*
    double part1 = ((-z * (-1/2) * (-2 * z) * pow(1 - 2 * rho * z + z * z, -3/2)) * (1 - rho) + (-1 + -z / sqrt(1 - 2 * rho * z + z * z)) * (-1) + 0.5 * (-2 * z) / (sqrt(1 - 2 * rho * z + z * z) + z - rho)) * ((1 - rho) * (1 - rho));
    part1 -= ((-1 + -z / sqrt(1 - 2 * rho * z + z * z)) * (1 - rho) + (sqrt(1 - 2 * rho * z + z * z) + z - rho)) * 2 * (-1) * (1 - rho);
    part1 *= (1 - rho) / (sqrt(1 - 2 * rho * z + z * z) + z - rho) / ((1 - rho) * (1 - rho) * (1 - rho) * (1 - rho));
    double part2 = (-1 * (sqrt(1 - 2 * rho * z + z * z) + z - rho) - (1 - rho) * (0.5 * (2 * z) / sqrt(1 - 2 * rho * z + z * z) - 1)) / ((sqrt(1 - 2 * rho * z + z * z) + z - rho) * (sqrt(1 - 2 * rho * z + z * z) + z - rho));
    part2 *= (((-1 + -z / sqrt(1 - 2 * rho * z + z * z)) * (1 - rho) + (sqrt(1 - 2 * rho * z + z * z) + z - rho)) / ((1 - rho) * (1 - rho)));
    */

    // Rewriting
    double fRho = 1 - 2 * rho * z + z * z;
    double dfRhodRho = -2 * z;
    double racineFRho = sqrt(fRho);
    double dracineFRhodRho = 0.5 * dfRhodRho / racineFRho;
    //dxzdRho = (-1 + -z / racineFRho) / (racineFRho + z - rho);
    //dxzdRho += (racineFRho + z - rho) / ((1 - rho) * (racineFRho + z - rho));

    double d2xzdRhodRho = ((z * dracineFRhodRho / racineFRho / racineFRho) * (racineFRho + z - rho) - (-1 + -z / racineFRho) * (dracineFRhodRho - 1)) / ((racineFRho + z - rho) * (racineFRho + z - rho));
    d2xzdRhodRho += ((dracineFRhodRho - 1) * ((1 - rho) * (racineFRho + z - rho)) - (racineFRho + z - rho) * ((-1 * (racineFRho + z - rho) + (1 - rho) * (dracineFRhodRho - 1)))) / ((1 - rho) * (racineFRho + z - rho) * (1 - rho) * (racineFRho + z - rho));

    // double d2xzdRhodRho = part1 + part2;
    return d2xzdRhodRho;
}

double Sabr::dxzdVovol(double strike){
    double z = this->z(strike);
    double dxzdVovol = this->dxzdz(strike, z) * this->dzdVovol(strike);
    return dxzdVovol;
}

double Sabr::d2xzdVovoldAlpha(double strike){
    double z = this->z(strike);
    double d2xzdVovoldAlpha = this->d2xzdzdz(strike, z) * this->dzdVovol(strike) * this->dzdAlpha(strike);
    d2xzdVovoldAlpha += this->dxzdz(strike, z) * this->d2zdVovoldAlpha(strike);
    return d2xzdVovoldAlpha;
}

double Sabr::d2xzdVovoldBeta(double strike){
    double z = this->z(strike);
    double d2xzdVovoldBeta = this->d2xzdzdz(strike, z) * this->dzdVovol(strike) * this->dzdBeta(strike);
    d2xzdVovoldBeta += this->dxzdz(strike, z) * this->d2zdVovoldBeta(strike);
    return d2xzdVovoldBeta;
}

double Sabr::d2xzdVovoldRho(double strike){
    double z = this->z(strike);
    double dxzdVovol = this->d2xzdzdRho(strike, z) * this->dzdVovol(strike);
    dxzdVovol += this->dxzdz(strike, z) * this->d2zdVovoldRho(strike); // =0
    return dxzdVovol;
}

double Sabr::d2xzdVovoldVovol(double strike){
    double z = this->z(strike);
    double d2xzdVovoldVovol = this->d2xzdzdz(strike, z) * this->dzdVovol(strike) * this->dzdVovol(strike);
    d2xzdVovoldVovol += this->dxzdz(strike, z) * this->d2zdVovoldVovol(strike);
    return d2xzdVovoldVovol;
}

double Sabr::dxzdForward(double strike){
    double z = this->z(strike);
    double dxzdForward = this->dxzdz(strike, z) * this->dzdForward(strike);
    return dxzdForward;
}

double Sabr::d2xzdForwarddAlpha(double strike){
    double z = this->z(strike);
    double d2xzdForwarddAlpha = this->d2xzdzdz(strike, z) * this->dzdForward(strike) * this->dzdAlpha(strike);
    d2xzdForwarddAlpha += this->dxzdz(strike, z) * this->d2zdForwarddAlpha(strike);
    return d2xzdForwarddAlpha;
}

double Sabr::d2xzdForwarddBeta(double strike){
    double z = this->z(strike);
    double dxzdForward = this->d2xzdzdz(strike, z) * this->dzdForward(strike) * this->dzdBeta(strike);
    dxzdForward += this->dxzdz(strike, z) * this->d2zdForwarddBeta(strike);
    return dxzdForward;
}

double Sabr::d2xzdForwarddRho(double strike){
    double z = this->z(strike);
    double dxzdForward = this->d2xzdzdRho(strike, z) * this->dzdForward(strike);
    dxzdForward += this->dxzdz(strike, z) * this->d2zdForwarddRho(strike); // =0
    return dxzdForward;
}

double Sabr::d2xzdForwarddVovol(double strike){
    double z = this->z(strike);
    double d2xzdForwarddVovol = this->d2xzdzdz(strike, z) * this->dzdForward(strike) * this->dzdVovol(strike);
    d2xzdForwarddVovol += this->dxzdz(strike, z) * this->d2zdForwarddVovol(strike);
    return d2xzdForwarddVovol;
}

double Sabr::d2xzdForwarddForward(double strike){
    double z = this->z(strike);
    double d2xzdForwarddForward = this->d2xzdzdz(strike, z) * this->dzdForward(strike) * this->dzdForward(strike);
    d2xzdForwarddForward += this->dxzdz(strike, z) * this->d2zdForwarddForward(strike);
    return d2xzdForwarddForward;
}

double Sabr::zxz(double strike){
    double z = this->z(strike);

    double zxz;
    if(fabs(forward - strike) > 0.00000000000001){
        zxz = z / this->xz(strike, z);
    }
    else{
        zxz = 1 - rho / 2 * z + (2 - 3 * rho * rho) / 12 * z * z;
    }

    return zxz;
}

double Sabr::dzxzdAlpha(double strike){
    double dzxzdAlpha;
    double z = this->z(strike);
    double dzdAlpha = this->dzdAlpha(strike);
    double xz = this->xz(strike, z);
    double dxzdAlpha = this->dxzdAlpha(strike);

    if(fabs(forward - strike) > 0.00000000000001){
        dzxzdAlpha = (dzdAlpha * xz - dxzdAlpha * z) / (xz * xz);
    }
    else{
        dzxzdAlpha = - rho / 2 * dzdAlpha + (2 - 3 * rho * rho) / 6 * z * dzdAlpha;
    }

    return dzxzdAlpha;
}

double Sabr::d2zxzdAlphadAlpha(double strike){
    double d2zxzdAlphadAlpha;
    double z = this->z(strike);
    double dzdAlpha = this->dzdAlpha(strike);
    double d2zdAlphadAlpha = this->d2zdAlphadAlpha(strike);
    double xz = this->xz(strike, z);
    double dxzdAlpha = this->dxzdAlpha(strike);
    double d2xzdAlphadAlpha = this->d2xzdAlphadAlpha(strike);

    if(fabs(forward - strike) > 0.00000000000001){
        d2zxzdAlphadAlpha = d2zdAlphadAlpha / xz;
        d2zxzdAlphadAlpha -= d2xzdAlphadAlpha * z / (xz * xz);
        d2zxzdAlphadAlpha -= 2 * dzdAlpha * dxzdAlpha / (xz * xz);
        d2zxzdAlphadAlpha += 2 * dxzdAlpha * dxzdAlpha * z / (xz * xz * xz);
    }
    else{
        d2zxzdAlphadAlpha = - rho / 2 * d2zdAlphadAlpha + (2 - 3 * rho * rho) / 6 * (z * d2zdAlphadAlpha + dzdAlpha * dzdAlpha);
    }

    return d2zxzdAlphadAlpha;
}

double Sabr::dzxzdBeta(double strike){
    double dzxzdBeta;

    double z = this->z(strike);
    double dzdBeta = this->dzdBeta(strike);

    double xz = this->xz(strike, z);
    double dxzdBeta = this->dxzdBeta(strike);

    if(fabs(forward - strike) > 0.00000000000001){
        dzxzdBeta = (dzdBeta * xz - dxzdBeta * z) / (xz * xz);
    }
    else{
        dzxzdBeta =  - rho / 2 * dzdBeta + (2 - 3 * rho * rho) / 6 * z * dzdBeta;
    }

    return dzxzdBeta;
}

double Sabr::d2zxzdBetadAlpha(double strike){
    double d2zxzdBetadAlpha;

    double z = this->z(strike);
    double dzdBeta = this->dzdBeta(strike);
    double dzdAlpha = this->dzdAlpha(strike);
    double d2zdBetadAlpha = this->d2zdBetadAlpha(strike);

    double xz = this->xz(strike, z);
    double dxzdBeta = this->dxzdBeta(strike);
    double dxzdAlpha = this->dxzdAlpha(strike);
    double d2xzdBetadAlpha = this->d2xzdBetadAlpha(strike);

    if(fabs(forward - strike) > 0.00000000000001){
        d2zxzdBetadAlpha = ((d2zdBetadAlpha * xz + dzdBeta * dxzdAlpha - d2xzdBetadAlpha * z - dxzdBeta * dzdAlpha) * (xz * xz) - (dzdBeta * xz - dxzdBeta * z) * ( 2 * dxzdAlpha * xz)) / (xz * xz * xz * xz);
    }
    else{
        d2zxzdBetadAlpha =  - rho / 2 * d2zdBetadAlpha + (2 - 3 * rho * rho) / 6 * z * d2zdBetadAlpha;
    }

    return d2zxzdBetadAlpha;
}

double Sabr::d2zxzdAlphadBeta(double strike){
    return this->d2zxzdBetadAlpha(strike);
}

double Sabr::d2zxzdAlphadRho(double strike){
    return this->d2zxzdRhodAlpha(strike);
}

double Sabr::d2zxzdBetadBeta(double strike){
    double d2zxzdBetadBeta;

    double z = this->z(strike);
    double dzdBeta = this->dzdBeta(strike);
    double d2zdBetadBeta = this->d2zdBetadBeta(strike);

    double xz = this->xz(strike, z);
    double dxzdBeta = this->dxzdBeta(strike);
    double d2xzdBetadBeta = this->d2xzdBetadBeta(strike);

    if(fabs(forward - strike) > moneynessDiscontinuityThreshold){
        d2zxzdBetadBeta = ((d2zdBetadBeta * xz + dzdBeta * dxzdBeta - d2xzdBetadBeta * z - dxzdBeta * dzdBeta) * (xz * xz) - (dzdBeta * xz - dxzdBeta * z) * ( 2 * dxzdBeta * xz))/ (xz * xz * xz * xz);
    }
    else{
        d2zxzdBetadBeta =  - rho / 2 * d2zdBetadBeta + (2 - 3 * rho * rho) / 6 * z * d2zdBetadBeta;
    }

    return d2zxzdBetadBeta;
}

double Sabr::d2zxzdBetadRho(double strike){
    return this->d2zxzdRhodBeta(strike);
}

double Sabr::dzxzdRho(double strike){
    double dzxzdRho;

    double z = this->z(strike);

    double xz = this->xz(strike, z);
    double dxzdRho = this->dxzdRho(strike);

    if(fabs(forward - strike) > moneynessDiscontinuityThreshold){
        dzxzdRho = - dxzdRho * z / (xz * xz);
    }
    else{
        dzxzdRho =  - z / 2 + (2 - 6 * rho) / 12 * z * z;
    }

    return dzxzdRho;
}

double Sabr::d2zxzdRhodAlpha(double strike){
    double d2zxzdRhodAlpha;

    double z = this->z(strike);
    double dzdAlpha = this->dzdAlpha(strike);

    double xz = this->xz(strike, z);
    double dxzdRho = this->dxzdRho(strike);
    double dxzdAlpha = this->dxzdAlpha(strike);
    double d2xzdRhodAlpha = this->d2xzdRhodAlpha(strike);

    if(fabs(forward - strike) > moneynessDiscontinuityThreshold){
        d2zxzdRhodAlpha = - ((d2xzdRhodAlpha * z + dxzdRho * dzdAlpha) * xz * xz - dxzdRho * z * 2 * xz * dxzdAlpha) / (xz * xz * xz * xz);
    }
    else{
        d2zxzdRhodAlpha = - dzdAlpha / 2 + (2 - 6 * rho) / 12 * 2 * z * dzdAlpha;
    }

    return d2zxzdRhodAlpha;
}

double Sabr::d2zxzdRhodBeta(double strike){
    double d2zxzdRhodBeta;

    double z = this->z(strike);
    double dzdBeta = this->dzdBeta(strike);

    double xz = this->xz(strike, z);
    double dxzdBeta = this->dxzdBeta(strike);
    double dxzdRho = this->dxzdRho(strike);
    double d2xzdRhodBeta = this->d2xzdRhodBeta(strike);

    if(fabs(forward - strike) > moneynessDiscontinuityThreshold){
        d2zxzdRhodBeta = - ((d2xzdRhodBeta * z + dxzdRho * dzdBeta) * xz * xz - dxzdRho * z * 2 * dxzdBeta * xz) / (xz * xz * xz * xz);
    }
    else{
        d2zxzdRhodBeta =  - dzdBeta / 2 + (2 - 6 * rho) / 12 * 2 * dzdBeta * z;
    }

    return d2zxzdRhodBeta;
}

double Sabr::d2zxzdRhodRho(double strike){
    double d2zxzdRhodRho;

    double z = this->z(strike);
    double dzdRho = this->dzdRho(strike);

    double xz = this->xz(strike, z);
    double dxzdRho = this->dxzdRho(strike);
    double d2xzdRhodRho = this->d2xzdRhodRho(strike);

    if(fabs(forward - strike) > moneynessDiscontinuityThreshold){
        d2zxzdRhodRho = - ((d2xzdRhodRho * z + dxzdRho * dzdRho) * xz * xz - dxzdRho * z * 2 * dxzdRho * xz) / (xz * xz * xz * xz);
    }
    else{
        d2zxzdRhodRho =  - dzdRho / 2 + (-6 * z * z + (2 - 6 * rho) * 2 * dzdRho * z) / 12;
    }

    return d2zxzdRhodRho;
}

double Sabr::dzxzdVovol(double strike){
    double dzxzdVovol;

    double z = this->z(strike);
    double dzdVovol = this->dzdVovol(strike);

    double xz = this->xz(strike, z);
    double dxzdVovol = this->dxzdVovol(strike);

    if(fabs(forward - strike) > moneynessDiscontinuityThreshold){
        dzxzdVovol = (dzdVovol * xz - dxzdVovol * z) / (xz * xz);
    }
    else{
        dzxzdVovol =  - rho / 2 * dzdVovol + (2 - 3 * rho * rho) / 6 * z * dzdVovol;
    }

    return dzxzdVovol;
}

double Sabr::d2zxzdVovoldAlpha(double strike){
    double d2zxzdVovoldAlpha;

    double z = this->z(strike);
    double dzdAlpha = this->dzdAlpha(strike);
    double dzdVovol = this->dzdVovol(strike);
    double d2zdVovoldAlpha = this->d2zdVovoldAlpha(strike);

    double xz = this->xz(strike, z);
    double dxzdAlpha = this->dxzdAlpha(strike);
    double dxzdVovol = this->dxzdVovol(strike);
    double d2xzdVovoldAlpha = this->d2xzdVovoldAlpha(strike);

    if(fabs(forward - strike) > moneynessDiscontinuityThreshold){
        d2zxzdVovoldAlpha = (d2zdVovoldAlpha * xz + dzdVovol * dxzdAlpha - d2xzdVovoldAlpha * z - dxzdVovol * dzdAlpha) * (xz * xz);
        d2zxzdVovoldAlpha -= (dzdVovol * xz - dxzdVovol * z) * (2 * xz * dxzdAlpha);
        d2zxzdVovoldAlpha /= (xz * xz * xz * xz);
    }
    else{
        d2zxzdVovoldAlpha = - rho / 2 * d2xzdVovoldAlpha;
        d2zxzdVovoldAlpha += (2 - 3 * rho * rho) / 6 * z * d2zdVovoldAlpha;
        d2zxzdVovoldAlpha += (2 - 3 * rho * rho) / 6 * dzdAlpha * dzdVovol;
    }

    return d2zxzdVovoldAlpha;
}

double Sabr::d2zxzdVovoldBeta(double strike){
    double d2zxzdVovoldBeta;

    double z = this->z(strike);
    double dzdBeta = this->dzdBeta(strike);
    double dzdVovol = this->dzdVovol(strike);
    double d2zdVovoldBeta = this->d2zdVovoldBeta(strike);

    double xz = this->xz(strike, z);
    double dxzdBeta = this->dxzdBeta(strike);
    double dxzdVovol = this->dxzdVovol(strike);
    double d2xzdVovoldBeta = this->d2xzdVovoldBeta(strike);

    if(fabs(forward - strike) > moneynessDiscontinuityThreshold){
        d2zxzdVovoldBeta = ((d2zdVovoldBeta * xz  + dzdVovol * dxzdBeta - d2xzdVovoldBeta * z - dxzdVovol * dzdBeta) * xz * xz - (dzdVovol * xz - dxzdVovol * z) * 2 * dxzdBeta * xz) / (xz * xz * xz * xz);
    }
    else{
        d2zxzdVovoldBeta =  - rho / 2 * d2zdVovoldBeta + (2 - 3 * rho * rho) / 6 * z * d2zdVovoldBeta;
    }

    return d2zxzdVovoldBeta;
}

double Sabr::d2zxzdVovoldRho(double strike){
    double dzxzdVovol;

    double z = this->z(strike);
    double dzdRho = this->dzdRho(strike);
    double dzdVovol = this->dzdVovol(strike);
    double d2zdVovoldRho = this->d2zdVovoldRho(strike);

    double xz = this->xz(strike, z);
    double dxzdRho = this->dxzdRho(strike);
    double dxzdVovol = this->dxzdVovol(strike);
    double d2xzdVovoldRho = this->d2xzdVovoldRho(strike);

    if(fabs(forward - strike) > moneynessDiscontinuityThreshold){
        dzxzdVovol = ((d2zdVovoldRho * xz + dzdVovol * dxzdRho - d2xzdVovoldRho * z - dxzdVovol * dzdRho) * xz * xz - (dzdVovol * xz - dxzdVovol * z) * 2 * dxzdRho * xz) / (xz * xz * xz * xz);
    }
    else{
        dzxzdVovol =  - rho / 2 * d2xzdVovoldRho + (2 - 3 * rho * rho) / 6 * (dzdRho * dzdVovol + z * d2zdVovoldRho);
    }

    return dzxzdVovol;
}

double Sabr::d2zxzdVovoldVovol(double strike){
    double d2zxzdVovoldVovol;

    double z = this->z(strike);
    double dzdVovol = this->dzdVovol(strike);
    double d2zdVovoldVovol = this->d2zdVovoldVovol(strike);

    double xz = this->xz(strike, z);
    double dxzdVovol = this->dxzdVovol(strike);
    double d2xzdVovoldVovol = this->d2xzdVovoldVovol(strike);

    if(fabs(forward - strike) > moneynessDiscontinuityThreshold){
        d2zxzdVovoldVovol = ((d2zdVovoldVovol * xz + dzdVovol * dxzdVovol - d2xzdVovoldVovol * z - dxzdVovol * dzdVovol) * xz * xz - (dzdVovol * xz - dxzdVovol * z) * 2 * dxzdVovol * xz) / (xz * xz * xz * xz);
    }
    else{
        d2zxzdVovoldVovol =  - rho / 2 * d2zdVovoldVovol + (2 - 3 * rho * rho) / 6 * z * d2zdVovoldVovol;
    }

    return d2zxzdVovoldVovol;
}

double Sabr::dzxzdForward(double strike){
    double z = this->z(strike);
    double dzdForward = this->dzdForward(strike);

    double xz = this->xz(strike, z);
    double dxzdForward = this->dxzdForward(strike);

    double dzxzdForward;
    if(fabs(forward - strike) > moneynessDiscontinuityThreshold){
        dzxzdForward = (dzdForward * xz - z * dxzdForward) / (xz * xz);
    }
    else{
        dzxzdForward = - rho / 2 * dzdForward + (2 - 3 * rho * rho) / 12 * 2 * dzdForward * z;
    }

    return dzxzdForward;
}

double Sabr::d2zxzdForwarddAlpha(double strike){
    double z = this->z(strike);
    double dzdForward = this->dzdForward(strike);
    double dzdAlpha = this->dzdAlpha(strike);
    double d2zdForwarddAlpha = this->d2zdForwarddAlpha(strike);

    double xz = this->xz(strike, z);
    double dxzdAlpha = this->dxzdAlpha(strike);
    double dxzdForward = this->dxzdForward(strike);
    double d2xzdForwarddAlpha = this->d2xzdForwarddAlpha(strike);

    double d2zxzdForwarddAlpha;
    if(fabs(forward - strike) > moneynessDiscontinuityThreshold){
        d2zxzdForwarddAlpha = (d2zdForwarddAlpha * xz + dzdForward * dxzdAlpha - dzdAlpha * dxzdForward - z * d2xzdForwarddAlpha) * (xz * xz);
        d2zxzdForwarddAlpha -= (dzdForward * xz - z * dxzdForward) * (2 * dxzdAlpha * xz);
        d2zxzdForwarddAlpha /= (xz * xz * xz * xz);
    }
    else{
        d2zxzdForwarddAlpha = - rho / 2 * d2zdForwarddAlpha + (2 - 3 * rho * rho) / 12 * 2 * (d2zdForwarddAlpha * z + dzdForward * dzdAlpha);
    }

    return d2zxzdForwarddAlpha;
}

double Sabr::d2zxzdForwarddBeta(double strike){
    double z = this->z(strike);
    double dzdBeta = this->dzdBeta(strike);
    double dzdForward = this->dzdForward(strike);
    double d2zdForwarddBeta = this->d2zdForwarddBeta(strike);

    double xz = this->xz(strike, z);
    double dxzdBeta = this->dxzdBeta(strike);
    double dxzdForward = this->dxzdForward(strike);
    double d2xzdForwarddBeta = this->d2xzdForwarddBeta(strike);

    double d2zxzdForwarddBeta;
    if(fabs(forward - strike) > moneynessDiscontinuityThreshold){
        d2zxzdForwarddBeta = ((d2zdForwarddBeta * xz + dzdForward * dxzdBeta - dzdBeta * dxzdForward - z * d2xzdForwarddBeta) * xz * xz - (dzdForward * xz - z * dxzdForward) * 2 * dxzdBeta * xz) / (xz * xz * xz * xz);
    }
    else{
        d2zxzdForwarddBeta = - rho / 2 * d2xzdForwarddBeta + (2 - 3 * rho * rho) / 12 * 2 * (d2xzdForwarddBeta * z + dzdForward * dzdBeta);
    }

    return d2zxzdForwarddBeta;
}

double Sabr::d2zxzdForwarddRho(double strike){
    double z = this->z(strike);
    double dzdRho = this->dzdRho(strike);
    double dzdForward = this->dzdForward(strike);
    double d2zdForwarddRho = this->d2zdForwarddRho(strike);

    double xz = this->xz(strike, z);
    double dxzdRho = this->dxzdRho(strike);
    double dxzdForward = this->dxzdForward(strike);
    double d2xzdForwarddRho = this->d2xzdForwarddRho(strike);

    double d2zxzdForwarddRho;
    if(fabs(forward - strike) > moneynessDiscontinuityThreshold){
        d2zxzdForwarddRho = ((d2zdForwarddRho * xz + dzdForward * dxzdRho - dzdRho * dxzdForward - z * d2xzdForwarddRho) * xz * xz - (dzdForward * xz - z * dxzdForward) * 2 * dxzdRho * xz) / (xz * xz * xz * xz);
    }
    else{
        d2zxzdForwarddRho = - rho / 2 * d2zdForwarddRho + (2 - 3 * rho * rho) / 12 * 2 * (d2zdForwarddRho * z + dzdForward * dzdRho);
    }

    return d2zxzdForwarddRho;
}

double Sabr::d2zxzdForwarddVovol(double strike){
    double z = this->z(strike);
    double dzdVovol = this->dzdVovol(strike);
    double dzdForward = this->dzdForward(strike);
    double d2zdForwarddVovol = this->d2zdForwarddVovol(strike);

    double xz = this->xz(strike, z);
    double dxzdVovol = this->dxzdVovol(strike);
    double dxzdForward = this->dxzdForward(strike);
    double d2xzdForwarddVovol = this->d2xzdForwarddVovol(strike);

    double d2zxzdForwarddVovol;
    if(fabs(forward - strike) > moneynessDiscontinuityThreshold){
        d2zxzdForwarddVovol = ((d2zdForwarddVovol * xz + dzdForward * dxzdVovol - dzdVovol * dxzdForward - z * d2xzdForwarddVovol) * xz * xz - (dzdForward * xz - z * dxzdForward) * 2 * dxzdVovol * xz) / (xz * xz * xz * xz);
    }
    else{
        d2zxzdForwarddVovol = - rho / 2 * d2zdForwarddVovol + (2 - 3 * rho * rho) / 12 * 2 * (d2zdForwarddVovol * z + dzdForward * dzdVovol);
    }

    return d2zxzdForwarddVovol;
}

double Sabr::d2zxzdForwarddForward(double strike){
    double z = this->z(strike);
    double dzdForward = this->dzdForward(strike);
    double d2zdForwarddForward = this->d2zdForwarddForward(strike);

    double xz = this->xz(strike, z);
    double dxzdForward = this->dxzdForward(strike);
    double d2xzdForwarddForward = this->d2xzdForwarddForward(strike);

    double d2zxzdForwarddForward;
    if(fabs(forward - strike) > moneynessDiscontinuityThreshold){
        d2zxzdForwarddForward = ((d2zdForwarddForward * xz + dzdForward * dxzdForward - dzdForward * dxzdForward - z * d2xzdForwarddForward) * xz * xz - (dzdForward * xz - z * dxzdForward) * 2 * dxzdForward * xz) / (xz * xz * xz * xz);
    }
    else{
        d2zxzdForwarddForward = - rho / 2 * d2zdForwarddForward + (2 - 3 * rho * rho) / 12 * 2 * (d2zdForwarddForward * z + dzdForward * dzdForward);
    }

    return d2zxzdForwarddForward;
}

double Sabr::Rf(double strike){
    double Rf = (1 - beta) * (1 - beta) / 24 * alpha * alpha / pow(forward * strike, 1 - beta);
    Rf += rho * beta * vovol * alpha / 4 / pow(forward * strike, (1 - beta) / 2);
    Rf += (2 - 3 * rho * rho) / 24 * vovol * vovol;
    return Rf;
}

double Sabr::dRfdAlpha(double strike){
    double dRfdAlpha = (1 - beta) * (1 - beta) / 12 * alpha / pow(forward * strike, 1 - beta);
    dRfdAlpha += rho * beta * vovol / 4 / pow(forward * strike, (1 - beta) / 2);
    return dRfdAlpha;
}

double Sabr::d2RfdAlphadAlpha(double strike){
    double d2RfdAlphadAlpha = (1 - beta) * (1 - beta) / 12 / pow(forward * strike, 1 - beta);
    return d2RfdAlphadAlpha;
}

double Sabr::dRfdBeta(double strike){
    double dRfdBeta = ((- 2 + 2 * beta) + (1 - beta) * (1 - beta) * log(forward * strike)) / pow(forward * strike, 1 - beta) / 24 * alpha * alpha;
    dRfdBeta += (rho * vovol * alpha + rho * beta * vovol * alpha * log(forward * strike) / 2) / pow(forward * strike, (1 - beta) / 2) / 4;
    return dRfdBeta;
}

double Sabr::d2RfdBetadAlpha(double strike){
    double d2RfdBetadAlpha = ((- 2 + 2 * beta) * pow(forward * strike, 1 - beta) + (1 - beta) * (1 - beta) * log(forward * strike) * pow(forward * strike, 1 - beta)) / pow(forward * strike, 2 - 2 * beta) / 24 * 2 * alpha;
    d2RfdBetadAlpha += (rho * vovol * pow(forward * strike, (1 - beta) / 2) + rho * beta * vovol * log(forward * strike) / 2 * pow(forward * strike, (1 - beta) / 2)) / 4 / pow(forward * strike, 1 - beta);
    return d2RfdBetadAlpha;
}

double Sabr::d2RfdBetadBeta(double strike){
    double part1 = (2 - 2 * (1 - beta) * log(forward * strike)) * pow(forward * strike, beta - 1);
    part1 += ((- 2 + 2 * beta) + (1 - beta) * (1 - beta) * log(forward * strike)) * log(forward * strike) * pow(forward * strike, beta - 1);
    double part2 = (rho * vovol * alpha * log(forward * strike) / 2) * pow(forward * strike, (beta - 1) / 2);
    part2 += (rho * vovol * alpha + rho * beta * vovol * alpha * log(forward * strike) / 2) * 0.5 * log(forward * strike) * pow(forward * strike, (beta - 1) / 2);
    double d2RfdBetadBeta = part1 / 24 * alpha * alpha + part2 / 4;
    return d2RfdBetadBeta;
}

double Sabr::d2RfdBetadRho(double strike){
    return this->d2RfdRhodBeta(strike);
}

double Sabr::d2RfdAlphadBeta(double strike){
    return this->d2RfdBetadAlpha(strike);
}

double Sabr::d2RfdAlphadRho(double strike){
    return this->d2RfdRhodAlpha(strike);
}

double Sabr::dRfdRho(double strike){
    double dRfdRho = beta * vovol * alpha / 4 / pow(forward * strike, (1 - beta) / 2);
    dRfdRho += - rho / 4 * vovol * vovol;
    return dRfdRho;
}

double Sabr::d2RfdRhodAlpha(double strike){
    double d2RfdRhodAlpha = beta * vovol / 4 / pow(forward * strike, (1 - beta) / 2);
    return d2RfdRhodAlpha;
}

double Sabr::d2RfdRhodBeta(double strike){
    double d2RfdRhodBeta = vovol * alpha / 4 * pow(forward * strike, (beta - 1) / 2);
    d2RfdRhodBeta += beta * vovol * alpha / 4 * 0.5 * log(forward * strike) * pow(forward * strike, (beta - 1) / 2);
    return d2RfdRhodBeta;
}

double Sabr::d2RfdRhodRho(double strike){
    double dRfdRho = - 1 / 4 * vovol * vovol;
    return dRfdRho;
}

double Sabr::dRfdVovol(double strike){
    double dRfdVovol = rho * beta * alpha / 4 * pow(forward * strike, (beta - 1) / 2);
    dRfdVovol += (2 - 3 * rho * rho) / 12 * vovol;
    return dRfdVovol;
}

double Sabr::d2RfdVovoldAlpha(double strike){
    double d2RfdVovoldAlpha = rho * beta / 4 * pow(forward * strike, (beta - 1) / 2);
    return d2RfdVovoldAlpha;
}

double Sabr::d2RfdVovoldBeta(double strike){
    double d2RfdVovoldBeta = rho * alpha / 4 * pow(forward * strike, (beta - 1) / 2);
    d2RfdVovoldBeta += rho * beta * alpha / 4 * 0.5 * log(forward * strike) *  pow(forward * strike, (beta - 1) / 2);
    return d2RfdVovoldBeta;
}

double Sabr::d2RfdVovoldRho(double strike){
    double d2RfdVovoldRho = beta * alpha / 4 * pow(forward * strike, (beta - 1) / 2);
    d2RfdVovoldRho += (- 6 * rho) / 12 * vovol;
    return d2RfdVovoldRho;
}

double Sabr::d2RfdVovoldVovol(double strike){
    double d2RfdVovoldVovol = (2 - 3 * rho * rho) / 12;
    return d2RfdVovoldVovol;
}

double Sabr::dRfdForward(double strike){
    double dRfdForward = (1 - beta) * (1 - beta) / 24 * alpha * alpha * (beta - 1) * pow(forward * strike, beta - 2) * strike;
    dRfdForward += rho * beta * vovol * alpha / 4 * (beta - 1) / 2 * pow(forward * strike, (beta - 3) / 2) * strike;
    return dRfdForward;
}

double Sabr::d2RfdForwarddAlpha(double strike){
    double d2RfdForwarddAlpha = (1 - beta) * (1 - beta) / 24 * 2 * alpha * (beta - 1) * pow(forward * strike, beta - 2) * strike;
    d2RfdForwarddAlpha += rho * beta * vovol / 4 * (beta - 1) / 2 * pow(forward * strike, (beta - 3) / 2) * strike;
    return d2RfdForwarddAlpha;
}

double Sabr::d2RfdForwarddBeta(double strike){
    double d2RfdForwarddBeta =  -2 * (1 - beta) / 24 * alpha * alpha * (beta - 1) * pow(forward * strike, beta - 2) * strike;
    d2RfdForwarddBeta += (1 - beta) * (1 - beta) / 24 * alpha * alpha * (1) * pow(forward * strike, beta - 2) * strike;
    d2RfdForwarddBeta += (1 - beta) * (1 - beta) / 24 * alpha * alpha * (beta - 1) * log(forward * strike) * pow(forward * strike, beta - 2) * strike;
    d2RfdForwarddBeta += rho * vovol * alpha / 4 * (beta - 1) / 2 * pow(forward * strike, (beta - 3) / 2) * strike;
    d2RfdForwarddBeta += rho * beta * vovol * alpha / 4 / 2 * pow(forward * strike, (beta - 3) / 2) * strike;
    d2RfdForwarddBeta += rho * beta * vovol * alpha / 4 * (beta - 1) / 2 * 0.5 * log(forward * strike) * pow(forward * strike, (beta - 3) / 2) * strike;
    return d2RfdForwarddBeta;
}

double Sabr::d2RfdForwarddRho(double strike){
    double d2RfdForwarddRho = beta * vovol * alpha / 4 * (beta - 1) / 2 * pow(forward * strike, (beta - 3) / 2) * strike;
    return d2RfdForwarddRho;
}

double Sabr::d2RfdForwarddVovol(double strike){
    double d2RfdForwarddVovol = rho * beta * alpha / 4 * (beta - 1) / 2 * pow(forward * strike, (beta - 3) / 2) * strike;
    return d2RfdForwarddVovol;
}

double Sabr::d2RfdForwarddForward(double strike){
    double d2RfdForwarddForward = (1 - beta) * (1 - beta) / 24 * alpha * alpha * (beta - 1) * strike * (beta - 2) * pow(forward * strike, beta - 3) * strike;
    d2RfdForwarddForward += rho * beta * vovol * alpha / 4 * (beta - 1) / 2 * strike * (beta - 3) / 2 * pow(forward * strike, (beta - 5) / 2) * strike;
    return d2RfdForwarddForward;
}
