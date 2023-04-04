#ifndef SMILE_CURVE_H
#define SMILE_CURVE_H

#define _USE_MATH_DEFINES

using namespace std;

double N(double x);
double Nprime(double x);

#define moneynessDiscontinuityThreshold 0.0000000000001
#define precisionMachine 0.0000000000001
#define numberOfParameters 4

class Sabr
{
public:
    // Constructor and destructor
    explicit Sabr();
    explicit Sabr(double alpha, double beta, double rho, double vovol, double forward, int t0, int maturity);
    virtual ~Sabr();

    // Getters
    int getCalDate() const;
    double getAlpha() const;
    double getBeta() const;
    double getRho() const;
    double getVovol() const;
    double getForward() const;
    double getMaturity() const;

    // Setters
    void setCalDate(int t0);
    void setAlpha(double alpha);
    void setBeta(double beta);
    void setRho(double rho);
    void setVovol(double vovol);
    void setForward(double forward);
    void setMaturity(int maturity);

    // Volatility
    double volatility(double strike);   // Formula A.69c in "Managing Smile Risk" by Hagan, Kumar, Lesniewski & Woodward.

    // First order derivatives
    double dSdAlpha(double strike);
    double dSdBeta(double strike);
    double dSdRho(double strike);
    double dSdVovol(double strike);
    double dSdForward(double strike);

    // Second order derivatives
    double d2SdAlphadAlpha(double strike);
    double d2SdBetadAlpha(double strike);
    double d2SdRhodAlpha(double strike);
    double d2SdVovoldAlpha(double strike);
    double d2SdForwarddAlpha(double strike);

    double d2SdAlphadBeta(double strike);
    double d2SdBetadBeta(double strike);
    double d2SdRhodBeta(double strike);
    double d2SdVovoldBeta(double strike);
    double d2SdForwarddBeta(double strike);     // 13% relative error sensi VS variation, but I cannot find any error in the formulas

    double d2SdAlphadRho(double strike);
    double d2SdBetadRho(double strike);
    double d2SdRhodRho(double strike);
    double d2SdVovoldRho(double strike);        // 20% relative error sensi VS variation, but I cannot find any error in the formulas
    double d2SdForwarddRho(double strike);

    double d2SdAlphadVovol(double strike);
    double d2SdBetadVovol(double strike);
    double d2SdRhodVovol(double strike);        // 20% relative error sensi VS variation, but I cannot find any error in the formulas
    double d2SdVovoldVovol(double strike);
    double d2SdForwarddVovol(double strike);

    double d2SdAlphadForward(double strike);
    double d2SdBetadForward(double strike);     // 13% relative error sensi VS variation, but I cannot find any error in the formulas
    double d2SdRhodForward(double strike);
    double d2SdVovoldForward(double strike);
    double d2SdForwarddForward(double strike);  // 2% relative error sensi VS variation, but I cannot find any error in the formulas

public:
    // Intermediates
    double W(double strike);
    double dWdBeta(double strike);
    double d2WdBetadBeta(double strike);
    double dWdForward(double strike);
    double d2WdForwarddBeta(double strike);
    double d2WdForwarddForward(double strike);

    double Lf(double strike);
    double dLfdBeta(double strike);
    double d2LfdBetadBeta(double strike);
    double dLfdForward(double strike);
    double d2LfdForwarddBeta(double strike);
    double d2LfdForwarddForward(double strike);

    double z(double strike);
    double dzdAlpha(double strike);
    double d2zdAlphadAlpha(double strike);
    double d2zdAlphadBeta(double strike);
    double d2zdAlphadRho(double strike);
    double dzdBeta(double strike);
    double d2zdBetadAlpha(double strike);
    double d2zdBetadBeta(double strike);
    double d2zdBetadRho(double strike);
    double dzdRho(double strike);
    double d2zdRhodAlpha(double strike);
    double d2zdRhodBeta(double strike);
    double d2zdRhodRho(double strike);
    double dzdVovol(double strike);
    double d2zdVovoldAlpha(double strike);
    double d2zdVovoldBeta(double strike);
    double d2zdVovoldRho(double strike);
    double d2zdVovoldVovol(double strike);
    double dzdForward(double strike);
    double d2zdForwarddAlpha(double strike);
    double d2zdForwarddBeta(double strike);
    double d2zdForwarddRho(double strike);
    double d2zdForwarddVovol(double strike);
    double d2zdForwarddForward(double strike);

    double xz(double strike, double z);
    double dxzdz(double strike, double z);
    double d2xzdzdz(double strike, double z);
    double d2xzdzdRho(double strike, double z);
    double dxzdAlpha(double strike);
    double d2xzdAlphadAlpha(double strike);
    double d2xzdAlphadBeta(double strike);
    double d2xzdAlphadRho(double strike);
    double dxzdBeta(double strike);
    double d2xzdBetadAlpha(double strike);
    double d2xzdBetadBeta(double strike);
    double d2xzdBetadRho(double strike);
    double dxzdRho(double strike);
    double d2xzdRhodAlpha(double strike);
    double d2xzdRhodBeta(double strike);
    double d2xzdRhodRho(double strike);
    double d2xzdRhodz(double strike, double z);
    double dxzdVovol(double strike);
    double d2xzdVovoldAlpha(double strike);
    double d2xzdVovoldBeta(double strike);
    double d2xzdVovoldRho(double strike);
    double d2xzdVovoldVovol(double strike);
    double dxzdForward(double strike);
    double d2xzdForwarddAlpha(double strike);
    double d2xzdForwarddBeta(double strike);
    double d2xzdForwarddRho(double strike);
    double d2xzdForwarddVovol(double strike);
    double d2xzdForwarddForward(double strike);

    double zxz(double strike);
    double dzxzdAlpha(double strike);
    double d2zxzdAlphadAlpha(double strike);
    double d2zxzdAlphadBeta(double strike);
    double d2zxzdAlphadRho(double strike);
    double dzxzdBeta(double strike);
    double d2zxzdBetadAlpha(double strike);
    double d2zxzdBetadBeta(double strike);
    double d2zxzdBetadRho(double strike);
    double dzxzdRho(double strike);
    double d2zxzdRhodAlpha(double strike);
    double d2zxzdRhodBeta(double strike);
    double d2zxzdRhodRho(double strike);
    double dzxzdVovol(double strike);
    double d2zxzdVovoldAlpha(double strike);
    double d2zxzdVovoldBeta(double strike);
    double d2zxzdVovoldRho(double strike);
    double d2zxzdVovoldVovol(double strike);
    double dzxzdForward(double strike);
    double d2zxzdForwarddAlpha(double strike);
    double d2zxzdForwarddBeta(double strike);
    double d2zxzdForwarddRho(double strike);
    double d2zxzdForwarddVovol(double strike);
    double d2zxzdForwarddForward(double strike);

    double Rf(double strike);
    double dRfdAlpha(double strike);
    double d2RfdAlphadAlpha(double strike);
    double d2RfdAlphadBeta(double strike);
    double d2RfdAlphadRho(double strike);
    double dRfdBeta(double strike);
    double d2RfdBetadAlpha(double strike);
    double d2RfdBetadBeta(double strike);
    double d2RfdBetadRho(double strike);
    double dRfdRho(double strike);
    double d2RfdRhodAlpha(double strike);
    double d2RfdRhodBeta(double strike);
    double d2RfdRhodRho(double strike);
    double dRfdVovol(double strike);
    double d2RfdVovoldAlpha(double strike);
    double d2RfdVovoldBeta(double strike);
    double d2RfdVovoldRho(double strike);
    double d2RfdVovoldVovol(double strike);
    double dRfdForward(double strike);
    double d2RfdForwarddAlpha(double strike);
    double d2RfdForwarddBeta(double strike);
    double d2RfdForwarddRho(double strike);
    double d2RfdForwarddVovol(double strike);
    double d2RfdForwarddForward(double strike);

    // Attributes
    double alpha = 0.001;   // control of the atm volatility
    double beta = 0.5;      // rotation of the smile
    double rho = -0.1;      // rotation of the smile
    double vovol = 0.1;     // convexity of the smile
    double forward = 0.01;  // forward of the underlying
    int maturity = 365;     // option maturity
    int calDate = 0;        // market data date
};

typedef double(Sabr::*computer)(double);
typedef double(Sabr::*computerIntermediate)(double, double);
typedef double(Sabr::*paramGetter)() const;
typedef void (Sabr::*paramSetter)(double);

#endif // SMILE_CURVE_H
