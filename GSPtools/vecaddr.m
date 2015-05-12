% VECADDR
% This function adds a set of vectors given in polar form
% All angles in *radians* (inputs and output), hence the R in VECADDR
% Usage: [M, ph] = vecaddr(M1, ph1, M2, ph2)

function [M, ph] = vecaddr(M1, ph1, M2, ph2)

avec = (M1*cos(ph1)+M2*cos(ph2))+j*(M1*sin(ph1)+M2*sin(ph2));
% M = abs(avec); ph = angle(avec);
% M = abs(avec); ph = phase(avec);

M = abs(avec); ph = atan2(imag(avec),real(avec)); 
