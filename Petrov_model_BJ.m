%% Fitting with Petrov & Schwille Equation (2025-02-14 by Yohan Lee)
% Reference: Petrov and Schwille. Biophys. J. 94, 2008, pages L41-L43.
close all
clear variables
clc
%% Constants
k = 1.38065 * 10^(-23) ; % Boltzmann constant [=] J/K
T = 293 ; % Temperature [=] K (= 20 degC)
eta_fluid = 0.001 ; % Water viscosity [=] Pa.s
gamma = 0.577215 ; % Euler constant
%% Petrov & Schwille Equation Constants
b1 = 2.74819 ;
b2 = 0.51465 ; % Corrected by the authors.Update. Biophys. J. 103, 2012, pg 375
c1 = 0.73761 ;
c2 = 0.52119 ;
%% Experimental Data (PUT YOUR DATA)
R = [0.75 ; 1; 2; 3; 4] ; % Condensate radius [=] um
D_exp = [0.292; 0.203; 0.108; 0.077; 0.051] ; % Diffusivity [=] um2/s
m = length (R) ;
%% Membrane viscosity candidates
n = 10000 ; % Number of candidates for membrane viscosity. Set reasonably large.
eta_memb = linspace(10^(-11), 10^(-7), n)' ; % Membrane viscosity [=] Pa.s.m.
%% Matrix for fitting results
eps = ones(m, n) ; % Epsilon (Reduced radius) = (2*eta_fluid*R / eta_memb). Dimensionless. 
for i = 1 : n
    for j = 1 : m
        eps(j, i) = 10^(-6) * (2 * eta_fluid * R (j)) / eta_memb (i) ; % Multiply 10^(-6) due to the Radius unit being um.
    end
end

D_fit = ones(m, n) ; % Matrix for D_fit(R) by Petrov Equation for each eta_memb.
for i = 1 : n
    for j = 1 : m
        D_fit (j, i) = (10^12 * k * T) / (4 * pi * eta_memb(i)) * (log(2 / eps(j, i)) - gamma + 4 * eps(j, i) / pi - ((eps(j, i)^2) / 2) * log(2 / eps(j, i))) / (1 - (eps(j, i)^3 / pi) * log(2 / eps(j, i)) + (c1 * eps(j, i)^b1) / (1 + c2 * eps(j, i)^b2)) ; % [=] um2/s
    end
end
%% Least square method to find the best-fit
% Squared Errors (an error being the difference between the measured and fitted value)
Sq_error = ones(m, n) ;
for i = 1 : n
    for j = 1 : m
        Sq_error (j, i) = (D_exp(j) - D_fit(j, i))^2 ; % Squared residuals
    end
end

% Sum of squared errors (SSE)
SSE = ones(n, 1) ;
for i = 1 : n
    SSE (i) = sum(Sq_error(:, i)) ; % Sum of the squared residuals
end

% Find the membrane viscosity that gives the smallest sum of the squared residuals.
bestfit_idx = find(SSE == min(SSE)) ;
eta_memb_bestfit = eta_memb(bestfit_idx) ; % Best-fit membrane viscosity

% Best-fit result
D_bestfit = D_fit(:, bestfit_idx) ;
%% Plot
logR = log10(R) ;
logD_exp = log10(D_exp) ;
logD_bestfit = log10(D_bestfit) ;

figure()
hold on
plot(logR, logD_exp, 'r--o')
plot(logR, logD_bestfit, 'b--o')
legend('Exp', 'Fitted')
xlabel ('log10_R (um)')
ylabel ('log10_D (um2/s)')
%% Statistics
yresid =  D_exp - D_bestfit ;
sq_yresid = yresid.^2 ; % Squared residuals (SR)
SSresid = sum(yresid.^2) ; % Sum of squared residuals (SSR)
SStotal = (m - 1) * var(D_exp) ;
R_squared = 1 - SSresid/SStotal ;
p = 1 ; % Number of predictors in the model (membrane viscosity)
RSE = sqrt(SSresid) / (m - p - 1) ; % Residual standard errors having the same units of the dependent variables.

%% Print results
fprintf('Best fitted membrane viscosity is %d Pa.s.m\n', eta_memb_bestfit) 
fprintf('R^2 = %f\n', R_squared)
fprintf('RSE = %f um2/s\n', RSE)