#ifndef SMILE_CURVE_H
#define SMILE_CURVE_H

#define _USE_MATH_DEFINES

using namespace std;

double N(double x);
double Nprime(double x);

class Sabr
{
public:
    // Constructor and destructor
    Sabr();
    Sabr(double alpha, double beta, double rho, double vovol, double forward, int t0, int maturity);
    ~Sabr();

    // Getters and Setters
    void setCalDate(int t0);
    void setAlpha(double alpha);
    void setBeta(double beta);
    void setRho(double rho);
    void setVovol(double vovol);
    void setForward(double forward);
    void setMaturity(int maturity);

    // Volatility
    double volatility(double strike);   // Formula A.69c in "Managing Smile Risk" by Hagan, Kumar, Lesniewski & Woodward.
    double dSdAlpha(double strike);
    double dSdBeta(double strike);
    double dSdRho(double strike);
    double dSdVovol(double strike);
    double dSdForward(double strike);

protected:
    // Intermediates
    double W(double strike);
    double dWdBeta(double strike);
    double dWdForward(double strike);
    double Lf(double strike);
    double dLfdBeta(double strike);

    double z(double strike);
    double dzdAlpha(double strike);
    double dzdBeta(double strike);
    double dzdVovol(double strike);

    double xz(double strike);
    double dxzdz(double strike);
    double dxzdAlpha(double strike);
    double dxzdBeta(double strike);
    double dxzdRho(double strike);
    double dxzdVovol(double strike);

    double zxz(double strike);
    double dzxzdAlpha(double strike);
    double dzxzdBeta(double strike);
    double dzxzdRho(double strike);
    double dzxzdVovol(double strike);
    double dzxzdForward(double strike);

    double Rf(double strike);
    double dRfdAlpha(double strike);
    double dRfdBeta(double strike);
    double dRfdRho(double strike);
    double dRfdVovol(double strike);
    double dRfdForward(double strike);

    // Attributes
    double alpha = 0.001;   // control of the atm volatility
    double beta = 0.5;      // rotation of the smile
    double rho = -0.1;      // rotation of the smile
    double vovol = 0.1;     // convexity of the smile
    double forward = 0.01;  // forward of the underlying
    int maturity = 365;     // option maturity
    int calDate = 0;        // market data date
};

#endif // SMILE_CURVE_H
