#include "lib_calibration.h"

Calibrator::Calibrator(){}

Calibrator::Calibrator(Sabr* myS, vector<double>* s, vector<double>* tV, double tMP, int lI){
    this->mySabr = myS;
    this->strikes = s;
    this->targetVolatilities = tV;
    this->targetMatchingPrecision = tMP;
    this->limitOfIterations = lI;

    this->numberOfPoints = s->size();
}

Calibrator::~Calibrator(){}

void Calibrator::calibrate(){
    gradient increment = {0};
    int iteration = 0;
    double matchingPrecision = this->normOfGradientOfError();

    while(matchingPrecision > targetMatchingPrecision && iteration < limitOfIterations){
        this->increment(&increment);
        mySabr->setAlpha(mySabr->getAlpha() + increment[0]);
        mySabr->setBeta(mySabr->getBeta() + increment[1]);
        mySabr->setRho(mySabr->getRho() + increment[2]);
        mySabr->setVovol(mySabr->getVovol() + increment[3]);

        iteration++;
        matchingPrecision = this->normOfGradientOfError();
    }
}

void Calibrator::fillStructure(vector<double> *s, int calibrationDate, int optionMaturity){
    this->strikes = s;
    this->mySabr->setCalDate(calibrationDate);
    this->mySabr->setMaturity(optionMaturity);
    this->numberOfPoints = s->size();
}

void Calibrator::fillMarketData(double forward, vector<double> *tV){
    this->targetVolatilities = tV;
    this->mySabr->setForward(forward);
}

void Calibrator::increment(gradient* gradIncrement){
    gradient grad;
    this->gradientOfError(&grad);
    hessian invHess;
    this->inverseHessianOfError(&invHess);

    *gradIncrement[0] = invHess[0][0] * grad[0] + invHess[0][1] * grad[1] + invHess[0][2] * grad[2] + invHess[0][3] * grad[3];
    *gradIncrement[1] = invHess[1][0] * grad[0] + invHess[1][1] * grad[1] + invHess[1][2] * grad[2] + invHess[1][3] * grad[3];
    *gradIncrement[2] = invHess[2][0] * grad[0] + invHess[2][1] * grad[1] + invHess[2][2] * grad[2] + invHess[2][3] * grad[3];
    *gradIncrement[3] = invHess[3][0] * grad[0] + invHess[3][1] * grad[1] + invHess[3][2] * grad[2] + invHess[3][3] * grad[3];

    return;
}

double Calibrator::computeError(){
    double error = 0;
    double ecart;

    for(int i = 0; i<numberOfPoints; i++){
        ecart = mySabr->volatility(strikes->at(i)) - targetVolatilities->at(i);
        error += ecart * ecart;
    }

    return error;
}

void Calibrator::gradientOfError(gradient* gradientError){
    double ecart;

    (*gradientError)[0] = 0;
    (*gradientError)[1] = 0;
    (*gradientError)[2] = 0;
    (*gradientError)[3] = 0;

    for(int i = 0; i<numberOfPoints; i++){
        ecart = mySabr->volatility(strikes->at(i)) - targetVolatilities->at(i);

        (*gradientError)[0] += 2 * mySabr->dSdAlpha(strikes->at(i)) * ecart;
        (*gradientError)[1] += 2 * mySabr->dSdBeta(strikes->at(i)) * ecart;
        (*gradientError)[2] += 2 * mySabr->dSdRho(strikes->at(i)) * ecart;
        (*gradientError)[3] += 2 * mySabr->dSdVovol(strikes->at(i)) * ecart;
    }

    return;
}

double Calibrator::normOfGradientOfError(){
    gradient grad;
    this->gradientOfError(&grad);
    return grad[0] * grad[0] + grad[1] * grad[1] + grad[2] * grad[2] + grad[3] * grad[3];
}

void Calibrator::hessianOfError(hessian* hessianError){
    double ecart;

    (*hessianError)[0][0] = 0;
    (*hessianError)[0][1] = 0;
    (*hessianError)[0][2] = 0;
    (*hessianError)[0][3] = 0;

    (*hessianError)[1][0] = 0;
    (*hessianError)[1][1] = 0;
    (*hessianError)[1][2] = 0;
    (*hessianError)[1][3] = 0;

    (*hessianError)[2][0] = 0;
    (*hessianError)[2][1] = 0;
    (*hessianError)[2][2] = 0;
    (*hessianError)[2][3] = 0;

    (*hessianError)[3][0] = 0;
    (*hessianError)[3][1] = 0;
    (*hessianError)[3][2] = 0;
    (*hessianError)[3][3] = 0;

    for(int i = 0; i<numberOfPoints; i++){
        ecart = mySabr->volatility(strikes->at(i)) - targetVolatilities->at(i);

        (*hessianError)[0][0] += 2 * (mySabr->d2SdAlphadAlpha(strikes->at(i)) * ecart + mySabr->dSdAlpha(strikes->at(i)) * mySabr->dSdAlpha(strikes->at(i)));
        (*hessianError)[0][1] += 2 * (mySabr->d2SdBetadAlpha(strikes->at(i)) * ecart + mySabr->dSdBeta(strikes->at(i)) * mySabr->dSdAlpha(strikes->at(i)));
        (*hessianError)[0][2] += 2 * (mySabr->d2SdRhodAlpha(strikes->at(i)) * ecart + mySabr->dSdRho(strikes->at(i)) * mySabr->dSdAlpha(strikes->at(i)));
        (*hessianError)[0][3] += 2 * (mySabr->d2SdVovoldAlpha(strikes->at(i)) * ecart + mySabr->dSdVovol(strikes->at(i)) * mySabr->dSdAlpha(strikes->at(i)));

        (*hessianError)[1][0] += 2 * (mySabr->d2SdAlphadBeta(strikes->at(i)) * ecart + mySabr->dSdAlpha(strikes->at(i)) * mySabr->dSdBeta(strikes->at(i)));
        (*hessianError)[1][1] += 2 * (mySabr->d2SdBetadBeta(strikes->at(i)) * ecart + mySabr->dSdBeta(strikes->at(i)) * mySabr->dSdBeta(strikes->at(i)));
        (*hessianError)[1][2] += 2 * (mySabr->d2SdRhodBeta(strikes->at(i)) * ecart + mySabr->dSdRho(strikes->at(i)) * mySabr->dSdBeta(strikes->at(i)));
        (*hessianError)[1][3] += 2 * (mySabr->d2SdVovoldBeta(strikes->at(i)) * ecart + mySabr->dSdVovol(strikes->at(i)) * mySabr->dSdBeta(strikes->at(i)));

        (*hessianError)[2][0] += 2 * (mySabr->d2SdAlphadRho(strikes->at(i)) * ecart + mySabr->dSdAlpha(strikes->at(i)) * mySabr->dSdRho(strikes->at(i)));
        (*hessianError)[2][1] += 2 * (mySabr->d2SdBetadRho(strikes->at(i)) * ecart + mySabr->dSdBeta(strikes->at(i)) * mySabr->dSdRho(strikes->at(i)));
        (*hessianError)[2][2] += 2 * (mySabr->d2SdRhodRho(strikes->at(i)) * ecart + mySabr->dSdRho(strikes->at(i)) * mySabr->dSdRho(strikes->at(i)));
        (*hessianError)[2][3] += 2 * (mySabr->d2SdVovoldRho(strikes->at(i)) * ecart + mySabr->dSdVovol(strikes->at(i)) * mySabr->dSdRho(strikes->at(i)));

        (*hessianError)[3][0] += 2 * (mySabr->d2SdAlphadVovol(strikes->at(i)) * ecart + mySabr->dSdAlpha(strikes->at(i)) * mySabr->dSdVovol(strikes->at(i)));
        (*hessianError)[3][1] += 2 * (mySabr->d2SdBetadVovol(strikes->at(i)) * ecart + mySabr->dSdBeta(strikes->at(i)) * mySabr->dSdVovol(strikes->at(i)));
        (*hessianError)[3][2] += 2 * (mySabr->d2SdRhodVovol(strikes->at(i)) * ecart + mySabr->dSdRho(strikes->at(i)) * mySabr->dSdVovol(strikes->at(i)));
        (*hessianError)[3][3] += 2 * (mySabr->d2SdVovoldVovol(strikes->at(i)) * ecart + mySabr->dSdVovol(strikes->at(i)) * mySabr->dSdVovol(strikes->at(i)));
    }

    return;
}

void Calibrator::inverseHessianOfError(hessian* inverseHessian){
    double buffer, dbuff;
    hessian hess;
    this->hessianOfError(&hess);

    (*inverseHessian)[0][0] = 1;
    (*inverseHessian)[0][1] = 0;
    (*inverseHessian)[0][2] = 0;
    (*inverseHessian)[0][3] = 0;

    (*inverseHessian)[1][0] = 0;
    (*inverseHessian)[1][1] = 1;
    (*inverseHessian)[1][2] = 0;
    (*inverseHessian)[1][3] = 0;

    (*inverseHessian)[2][0] = 0;
    (*inverseHessian)[2][1] = 0;
    (*inverseHessian)[2][2] = 1;
    (*inverseHessian)[2][3] = 0;

    (*inverseHessian)[3][0] = 0;
    (*inverseHessian)[3][1] = 0;
    (*inverseHessian)[3][2] = 0;
    (*inverseHessian)[3][3] = 1;

    int i;
    for (int j = 0; j < numberOfParameters; j++){
        i = j;
        while (hess[i][j] == 0 && i <= numberOfParameters) {
            i++;
        }
        if (i > numberOfParameters){
            std::cout << std::endl << "Error in inverseHessianOfError: the matrix is not inversible." << std::endl;
        }

        else{
            buffer = hess[j][0];
            hess[j][0] = hess[i][0];
            hess[i][0] = buffer;
            buffer = hess[j][1];
            hess[j][1] = hess[i][1];
            hess[i][1] = buffer;
            buffer = hess[j][2];
            hess[j][2] = hess[i][2];
            hess[i][2] = buffer;
            buffer = hess[j][3];
            hess[j][3] = hess[i][3];
            hess[i][3] = buffer;

            buffer = (*inverseHessian)[j][0];
            (*inverseHessian)[j][0] = (*inverseHessian)[i][0];
            (*inverseHessian)[i][0] = buffer;
            buffer = (*inverseHessian)[j][1];
            (*inverseHessian)[j][1] = (*inverseHessian)[i][1];
            (*inverseHessian)[i][1] = buffer;
            buffer = (*inverseHessian)[j][2];
            (*inverseHessian)[j][2] = (*inverseHessian)[i][2];
            (*inverseHessian)[i][2] = buffer;
            buffer = (*inverseHessian)[j][3];
            (*inverseHessian)[j][3] = (*inverseHessian)[i][3];
            (*inverseHessian)[i][3] = buffer;
        }

        dbuff = hess[j][j];

        buffer = hess[j][0];
        hess[j][0] = buffer + (1 - dbuff) / dbuff * buffer;
        buffer = hess[j][1];
        hess[j][1] = buffer + (1 - dbuff) / dbuff * buffer;
        buffer = hess[j][2];
        hess[j][2] = buffer + (1 - dbuff) / dbuff * buffer;
        buffer = hess[j][3];
        hess[j][3] = buffer + (1 - dbuff) / dbuff * buffer;

        buffer = (*inverseHessian)[j][0];
        (*inverseHessian)[j][0] = buffer + (1 - dbuff) / dbuff * buffer;
        buffer = (*inverseHessian)[j][1];
        (*inverseHessian)[j][1] = buffer + (1 - dbuff) / dbuff * buffer;
        buffer = (*inverseHessian)[j][2];
        (*inverseHessian)[j][2] = buffer + (1 - dbuff) / dbuff * buffer;
        buffer = (*inverseHessian)[j][3];
        (*inverseHessian)[j][3] = buffer + (1 - dbuff) / dbuff * buffer;

        for (int k = 0; k < numberOfParameters; k++){
            if (k != j){
                dbuff = -hess[k][j];
                hess[k][0] = hess[k][0] + dbuff * hess[j][0];
                hess[k][1] = hess[k][1] + dbuff * hess[j][1];
                hess[k][2] = hess[k][2] + dbuff * hess[j][2];
                hess[k][3] = hess[k][3] + dbuff * hess[j][3];

                (*inverseHessian)[k][0] = (*inverseHessian)[k][0] + dbuff * (*inverseHessian)[j][0];
                (*inverseHessian)[k][1] = (*inverseHessian)[k][1] + dbuff * (*inverseHessian)[j][1];
                (*inverseHessian)[k][2] = (*inverseHessian)[k][2] + dbuff * (*inverseHessian)[j][2];
                (*inverseHessian)[k][3] = (*inverseHessian)[k][3] + dbuff * (*inverseHessian)[j][3];
            }
        }
    }

    return;
}
