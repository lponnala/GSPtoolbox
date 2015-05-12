% UNBANOVA : "Unbalanced ANOVA"
% This function performs unbalanced ANOVA
% USAGE: P = unbanova(S)
% INPUTS: 
% S = (a x 2) matrix of summary statistics, rows representing groups
% S(k,:) = [samplesize_k, samplemean_k, samplevariance_k]
% OUTPUTS: 
% p-value for the test of 
% (H0: mu1 = mu2 = mu3 = ... = mua) Vs. (H1: all means are not equal)

function P = unbanova(S)

a = size(S,1);
N = S(:,1); M = S(:,2); V = S(:,3); 

MSE = ((N'*V)-sum(V))/(sum(N)-a); dfE = sum(N)-a;
OM = (N'*M)/sum(N);
MST = (N'*((M-OM*ones(size(M))).^2))/(a-1); dfT = a-1;

F = MST/MSE;
P = 1-fcdf(F,dfT,dfE);
