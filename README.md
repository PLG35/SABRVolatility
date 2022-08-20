# SABR Volatility
Library computing the volatility and its first order derivatives based on the article "Managing Smile Risk" by Hagan, Kumar, Lesniewski & Woodward. And in particular the formula A.69c.

## Content
This repository contains :
- the source code of the library, available under :
> src/

- the unitary tests validating the formulas of the first order derivatives, available under :
> tst/

## Roadmap
The first order derivative with regard to the forward is not ready yet. It would be useful to compute the DV01 of an swaption coming from the forward used in the interpolation of the volatility.

The second order derivatives are not computed yet. They would be useful to build a calibration tool that compute the parameters alpha, beta, rho and vovol, for a given forward and maturity, to best match market volatilities given as input.

A GUI to test scenarios and validate the calibration is under build.