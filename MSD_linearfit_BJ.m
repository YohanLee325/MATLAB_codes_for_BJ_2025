%% MSD over time-lag linear fitting (MSD = B*Tau) (2024. 8. 1. by Yohan Lee)
% This code is for linearly fitting each MSD curve by the least square method.
% From the slope of each fitted line, the diffusivity (D) is calculated.
% Slope = 4 * D. (MSD = 4D*Tau)
% Function file included: "linearfit.m"
close all
clear variables
clc
%% Variables
filename0 = 'MSD_example_BJ_uncoupled_bottom' ; % !!MAKE SURE IT IS THE RIGHT FILE!!
MSD = readmatrix(filename0,'Sheet', 1) ;
tau = readmatrix(filename0,'Sheet', 2) ;
n_tau = length(tau) ;
n_track = length(MSD(1, :)) ; % Number of tracks
%% Fitting MSD = b*Tau (zero-intercept)
D = ones(n_track, 1) ; % Diffusivity
x = tau ;
y = MSD ;

for i = 1 : n_track
    D(i) = linearfit(x, y(:, i)) / 4 ; % Function 'linearfit' is used.
end
D_avg = mean(D) ;
D_se = std(D) / sqrt(n_track) ; % Standard error
%% Final results summary
fprintf('Mean Diffusivity : %f um2/s\n', D_avg)