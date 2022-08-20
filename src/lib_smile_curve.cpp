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

// Getters & setters
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
    
    double sigma = alpha / W * zxz * (1 + dcf * this->Rf(strike));
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

    double dSdBeta = - alpha * dWdBeta / (W * W) * zxz * (1 + dcf * Rf) + alpha / W * dzxzdBeta * (1 + dcf * Rf) + alpha / W * zxz * dcf * this->dRfdBeta(strike);
    return dSdBeta;
}

double Sabr::dSdRho(double strike){
    double dcf = static_cast<double>(maturity - calDate) / 365.0;
    double W = this->W(strike);
    double zxz = this->zxz(strike);
    double dzxzdRho = this->dzxzdRho(strike);
    double dRfdRho = this->dRfdRho(strike);

    double dSdRho = alpha / W * (dzxzdRho * (1 + dcf * this->Rf(strike)) + zxz * dcf * dRfdRho);
    return dSdRho;
}

double Sabr::dSdVovol(double strike){
    double dcf = static_cast<double>(maturity - calDate) / 365.0;
    double W = this->W(strike);
    double zxz = this->zxz(strike);
    double dzxzdVovol = this->dzxzdVovol(strike);
    double dRfdVovol = this->dRfdVovol(strike);

    double dSdVovol = alpha / W * (dzxzdVovol * (1 + dcf * this->Rf(strike)) + zxz * dcf * dRfdVovol);
    return dSdVovol;
}

double Sabr::dSdForward(double strike){
    double dcf = static_cast<double>(maturity - calDate) / 365.0;
    double W = this->W(strike);
    double dWdForward = this->dWdForward(strike);
    double zxz = this->zxz(strike);
    double dzxzdForward = this->dzxzdForward(strike);
    double Rf = this->Rf(strike);

    double dSdForward = -alpha * dWdForward / W / W * zxz * (1 + dcf * Rf);
    dSdForward += alpha / W * dzxzdForward * (1 + dcf * Rf);
    dSdForward += alpha / W * zxz * dcf * this->dRfdForward(strike);

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

double Sabr::dWdForward(double strike){
    return 0;
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

double Sabr::z(double strike){
    double z = vovol / alpha * pow(forward * strike, (1 - beta) / 2) * log(forward / strike);
    return z;
}

double Sabr::dzdAlpha(double strike){
    double dzdAlpha =  - 1 / alpha * this->z(strike);
    return dzdAlpha;
}

double Sabr::dzdBeta(double strike){
    double dzdBeta = this->z(strike) * log(forward * strike) * (-0.5);
    return dzdBeta;
}

double Sabr::dzdVovol(double strike){
    double dzdVovol = 1 / alpha * pow(forward * strike, (1 - beta) / 2) * log(forward / strike);
    return dzdVovol;
}

double Sabr::xz(double strike){
    double z = this->z(strike);
    double xz = log((sqrt(1 - 2 * rho * z + z * z) + z - rho) / (1 - rho));
    return xz;
}

double Sabr::dxzdz(double strike){
    double z = this->z(strike);
    double dxzdz = (1 + (z - rho) / sqrt(1 - 2 * rho * z + z * z) ) / (sqrt(1 - 2 * rho * z + z * z) + z - rho);
    return dxzdz;
}

double Sabr::dxzdAlpha(double strike){
    double dxzdAlpha = this->dxzdz(strike) * this->dzdAlpha(strike);
    return dxzdAlpha;
}

double Sabr::dxzdBeta(double strike){
    double dxzBbeta = this->dxzdz(strike) * this->dzdBeta(strike);
    return dxzBbeta;
}

double Sabr::dxzdRho(double strike){
    double z = this->z(strike);
    double dxzdRho = (((-1 + -z / sqrt(1 - 2 * rho * z + z * z)) * (1 - rho) + (sqrt(1 - 2 * rho * z + z * z) + z - rho)) / ((1 - rho) * (1 - rho))) * (1 - rho) / (sqrt(1 - 2 * rho * z + z * z) + z - rho);
    return dxzdRho;
}

double Sabr::dxzdVovol(double strike){
    double dxzdVovol = this->dxzdz(strike) * this->dzdVovol(strike);
    return dxzdVovol;
}

double Sabr::zxz(double strike){
    double z = this->z(strike);

    double zxz;
    if(fabs(forward - strike) > 0.00000000000001){
        zxz = z / this->xz(strike);
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
    double xz = this->xz(strike);
    double dxzdAlpha = this->dxzdAlpha(strike);

    if(fabs(forward - strike) > 0.00000000000001){
        dzxzdAlpha = (dzdAlpha * xz - dxzdAlpha * z) / (xz * xz);
    }
    else{
        dzxzdAlpha = - rho / 2 * dzdAlpha + (2 - 3 * rho * rho) / 6 * z * dzdAlpha;
    }

    return dzxzdAlpha;
}

double Sabr::dzxzdBeta(double strike){
    double dzxzdBeta;

    double z = this->z(strike);
    double dzdBeta = this->dzdBeta(strike);

    double xz = this->xz(strike);
    double dxzdBeta = this->dxzdBeta(strike);

    if(fabs(forward - strike) > 0.00000000000001){
        dzxzdBeta = (dzdBeta * xz - dxzdBeta * z) / (xz * xz);
    }
    else{
        dzxzdBeta =  - rho / 2 * dzdBeta + (2 - 3 * rho * rho) / 6 * z * dzdBeta;
    }

    return dzxzdBeta;
}

double Sabr::dzxzdRho(double strike){
    double dzxzdRho;

    double z = this->z(strike);

    double xz = this->xz(strike);
    double dxzdRho = this->dxzdRho(strike);

    if(fabs(forward - strike) > 0.00000000000001){
        dzxzdRho = - dxzdRho * z / (xz * xz);
    }
    else{
        dzxzdRho =  - z / 2 + (2 - 6 * rho) / 12 * z * z;
    }

    return dzxzdRho;
}

double Sabr::dzxzdVovol(double strike){
    double dzxzdVovol;

    double z = this->z(strike);
    double dzdVovol = this->dzdVovol(strike);

    double xz = this->xz(strike);
    double dxzdVovol = this->dxzdVovol(strike);

    if(fabs(forward - strike) > 0.00000000000001){
        dzxzdVovol = (dzdVovol * xz - dxzdVovol * z) / (xz * xz);
    }
    else{
        dzxzdVovol =  - rho / 2 * dzdVovol + (2 - 3 * rho * rho) / 6 * z * dzdVovol;
    }

    return dzxzdVovol;
}

double Sabr::dzxzdForward(double strike){
    return 0;
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

double Sabr::dRfdBeta(double strike){
    double dRfdBeta = ((- 2 + 2 * beta) * pow(forward * strike, 1 - beta) + (1 - beta) * (1 - beta) * log(forward * strike) * pow(forward * strike, 1 - beta)) / pow(forward * strike, 2 - 2 * beta) / 24 * alpha * alpha;
    dRfdBeta += (rho * vovol * alpha * pow(forward * strike, (1 - beta) / 2) + rho * beta * vovol * alpha * log(forward * strike) / 2 * pow(forward * strike, (1 - beta) / 2)) / 4 / pow(forward * strike, 1 - beta);
    return dRfdBeta;
}

double Sabr::dRfdRho(double strike){
    double dRfdRho = beta * vovol * alpha / 4 / pow(forward * strike, (1 - beta) / 2);
    dRfdRho += - rho / 4 * vovol * vovol;
    return dRfdRho;
}

double Sabr::dRfdVovol(double strike){
    double dRfdVovol = rho * beta * alpha / 4 / pow(forward * strike, (1 - beta) / 2);
    dRfdVovol += (2 - 3 * rho * rho) / 12 * vovol;
    return dRfdVovol;
}

double Sabr::dRfdForward(double strike){
    return 0;
}
