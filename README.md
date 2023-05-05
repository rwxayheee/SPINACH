# SPINACH

This is a numerical interpolation-extrapolation-integration tool I wrote to analyze alchemical TI simulations. 

I made this object (cubic) for error propagation of integrals that are resulted from cubic spline interpolation. The W_extra attribute of this object gives the weight (contribution) of discrete f(x) in the integral. 

## Dependency
Python3, numpy.

## Method Explained
This is a small set of numerical integration tools I made for the analysis of alchemical thermodynamic integration (TI) simulations. The method is originally intended for simulations with dual-topology, single-step alchemical transformations that are achieved with soft-core potential. The numerical methods, although generally not as commonly used as the perturbation-based methods (e,g. BAR, see more at: doi:10.1007/s10822-015-9840-9), are more suitable for the above types of TI simulations because the outcome, DV/DL(L) (the derivative of potential energy function with respect to the coupling parameter), is more likely a smooth integrand.

To numerically compute the integral of DV/DL(L), literally, the free energy difference, from a set of discrete DV/DL(L) data, the following steps are taken as the general procedure: 
1. Interpolation
2. Extrapolation (if needed)
3. Integration
All are doable with respective functions in scipy. However, to compute the propagated errors within the integral, the integral needs to be written as a weighted sum of the DV/DL(L) (doi/10.1021/ct2003995).

The motive of this effort is to obtain the weights for error propagation. Here are some useful features of the weights: 
1. The weights are determined exclusively by the spacing of lambda. 
2. The formula of weights depends on the interpolation and extrapolation methods used. 
3. The weights by intervals are additive. 
4. The breakdown of weights (by order of the terms in polynomial interpolators) are additive. 

Segmental interpolation and extrapolation by polynomials have exact solutions to weights. The implementation of linear and cubic splines in this tool is mainly inspired by: 
The appendix of doi/10.1021/ct2003995l; 
https://github.com/MobleyLab/alchemical-analysis/blob/507e157b1a9658aa21c7bee100695ce66f2ad71d/alchemical_analysis/alchemical_analysis.py; and
Applied numerical methods with MATLAB® for engineers and scientists / Steven C. Chapra, Berger Chair in Computing and Engineering, Tufts University. 

In the interpolation step, it generally follows Chapra’s procedure to find the coefficients, with a tweak to this matrix equation (Chapra’s Eq 18.27):

TH C = TFD

The above equation is intended to solve coefficients c from finite differences. TH is a tridiagonal, NXN matrix that is made of the width of intervals, h. C is an NX1 column vector containing the unknown coefficients c. TFD is an NX1 column vector that is made of the finite differences. With a change of basis, the coefficients can be solved in the space of f(x). 

TH C = COB TF

The change-of-basis matrix COB, is NXN and it can be written with respect to the interpolators. If the natural cubic spline is used. TF is an NX1 column vector that is made of f(x). At this point, the matrix equation can be rearranged based on the invertibility of TH: 

C = TH-1 COB TF

The multiplication of the inverse of TH and COB gives the projection of coefficients c in the space of f(x). 

The coefficients a are set to be equal to the respective f(x) based on the continuity conditions. Coefficients b and d can be solved from coefficients c with the equality and smoothness conditions. This implementation solves b in a propragative manner. The coefficient in the first interval, b1, is computed first with Chapra’s Eq 18.21. The coefficient in the next interval is computed with Chapra’s Eq 18.23. The coefficients d are computed with Chapra’s Eq 18.18. All are solved in the space of f(x). 

In the interpolated segments, the segment-wise weights of f(x) can now be obtained by substituting the coefficients in the integral of each segment. The weights of f(x) of the integral can be obtained by summing up the segment-wise weights of f(x). If extrapolation is used, then the weights shall be updated accordingly. 

The cubic function uses natural cubic spline as the interpolator and it currently supports natural and linear extrapolation. The linear extrapolation essentially ignores the quartic components in the integral (or equivalently, the d is forced to be 0) of the extrapolated segments. The extrapolated segments also have c=0 in accordance to the natural end condition. 

Full doc with equations at: 
https://docs.google.com/document/d/1cdzjvWCLHGfLKSw0Qeun3R_36pQrMPxskibRasMrRYI/edit?usp=sharing
