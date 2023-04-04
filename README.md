# SABR Volatility

## Financial purpose
Library computing the volatility and its first order and second order derivatives based on the article  : "Managing Smile Risk" by Hagan, Kumar, Lesniewski & Woodward - and in particular the formula A.69c.

The derivative with regard to the forward are useful to compute the sensitivities of options coming from the forward used in the interpolation of the volatility.

Leveraging on the second order derivatives, a library of calibration permits to calibrate SABR parameters (alpha, beta, rho and vovol), based on an array of smile volatilities altrogether with an array of there corresponding strikes.

## Content
This repository contains :
- the source code of the two libraries, available under :
> src/

- the unitary tests validating the formulas of the first order and second order derivatives, available under :
> tst/

## Roadmap
The calibration algorithm is not tested yet. In a near future, it shall be :
 - unit tested,
 - backtested with market data representing various typology of market state.

To make sure the fitness of the calibrated smile is better around the ATM than on the tails, a weighted version of the function optimized during the calibration shall be implemented.

To avoid interference between the parameters beta and rho during the calibration, calibration modes with fixed beta or fixed rho shall be implemented.

A GUI to test scenarios and validate the calibration is planned.

Both interpolation and calibration functions can be encapsulated in an API.