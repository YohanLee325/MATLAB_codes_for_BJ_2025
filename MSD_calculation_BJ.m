%% Mean squared displacement (MSD) calculation (2024-07-23 coded by Yohan Lee)
% For plotting Figure 1e.
close all
clear variables
clc
%% Position data
filename0 = 'XY_position_example_uncoupled_bottom.xlsx' ; % MAKE SURE IT IS THE RIGHT FILE
X = readmatrix(filename0, 'Sheet', 'X') ; % X position
Y = readmatrix(filename0, 'Sheet', 'Y') ; % Y position
tau = readmatrix(filename0, 'Sheet', 'tau') ; % Tau (time lag). Zero tau not included.
n_time = length(X(:, 1)) ; % Number of time points
n_track = length(X(1, :)) ; % Number of tracks
%% Empty MSD matrix generation
msd_X = zeros((n_time - 1), n_track) ; % Empty MSD matrix
msd_Y = zeros((n_time - 1), n_track) ; % Empty MSD matrix
%% MSD calculation
square_disp = zeros (length(tau) , length(tau));

idx_x = n_time - 1 ;
for k = 1 : n_track
    idx_x = n_time - 1 ; % Reset idx_x
    for i = 1 : n_time - 1
        for j = 1 : idx_x
        square_disp(j, i) = (X(j, k) - X((j + i), k))^2 ;
        end
        msd_X (i, k) = sum(square_disp(: , i)) / idx_x; % MSD_X
        idx_x = idx_x - 1 ; 
    end
end

idx_y = n_time - 1 ;
for k = 1 : n_track
    idx_y = n_time - 1 ; % Reset idx_y
    for i = 1 : n_time - 1
        for j = 1 : idx_y
        square_disp(j, i) = (Y(j, k) - Y((j + i), k))^2 ;
        end
        msd_Y (i, k) = sum(square_disp(: , i)) / idx_y; % MSD_Y
        idx_y = idx_y - 1 ; 
    end
end

msd_sum = msd_X + msd_Y ; % FINAL MSD FOR EACH TRACK
%% Statistics
msd_mean = ones(length(tau), 1);
msd_std = ones(length(tau), 1); % standard deviation
msd_se = ones(length(tau), 1); % standard error

for i = 1 : length(tau)
    msd_mean(i) = mean(msd_sum(i, :));
end

for i = 1 : length(tau)
    msd_std(i) = std(msd_sum(i, :));
end

for i = 1 : length(tau)
    msd_se(i) = msd_std(i) / sqrt(n_track) ;
end
%% log_tau vs. log_MSD linear fitting to obtain the diffusion exponent.
% x = log10_tau, y = log10_MSD.
% Fitting to "y = a*x + b", a: diffusion exponent.
x = log10(tau) ;
y = log10(msd_mean) ;
[p, S] = polyfit(x, y, 1) ;

log_msd_se = ones(length(tau), 1) ; % Standard errors for log10-MSD 
for i = 1 : length(tau)
    log_msd_se(i) = (1/log(10))*(msd_se(i) / msd_mean(i)); % Error propagation.
end

f = polyval(p, x)' ;
figure(1)
plot(x, y, 'o')
errorbar(x, y, log_msd_se)
hold on
plot(x, f, '--')
xlabel ('log10-tau (s)')
ylabel ('log10-MSD (micron^2)')
legend('data', 'fitted')
%% Print results
fprintf('Diffusion exponent is %f\n', p(1)) 
fprintf('R^2 for linear fitting = %f\n', S.rsquared)