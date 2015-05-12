% EST_PAR: Estimate signal parameters (mean, amplitude, phase and SNR) 
% using trigonometric regression
% Usage: [mu,A,theta,SNR,err] = est_par(signal)
% INPUT: signal should be a row vector
% OUTPUTS:
% theta (phase-angle) is reported in degrees
% err=1: A is NaN, err=2: theta is NaN, err=3: SNR is NaN
% err=0 if (A,theta,SNR) are all valid numbers

function [mu,A,theta,SNR,err] = est_par(signal)

% ------Use trigonometric regression------
% Create regressors
N = length(signal);
I = (0:1:N-1)';
X = [ones(N,1) sin(2*pi*(1/3)*I) cos(2*pi*(1/3)*I)]; % works!!

alpha=0.05;
[b,bint,R,Rint,Stats] = REGRESS(signal',X,alpha); % Put data in a column vector

% ----Estimate mean----
mu = mean(signal);

% ----Estimate magnitude----
A = sqrt((b(2)^2) + (b(3)^2));

% ----Estimate phase using regression----
% This is correct since the coeff of b(3) is a sin(theta) and that of b(2) is a cos(theta)
% going by the sin(2.pi.f.t+theta) model fit
theta = (180/pi)*atan2(b(3),b(2)); 

% % ----Estimate phase using polar plot----
% [mag,phase]=cumm_mag_phase(signal); theta=(180/pi)*phase(end);

% Estimate SNR
SNR = 10*log10((A^2)/2) - 10*log10(mean((R-mean(R)*ones(size(R))).^2));

err=0;
if isnan(A)
    err=1;
elseif isnan(theta)
    err=2;
elseif isnan(SNR)
    err=3;
end
