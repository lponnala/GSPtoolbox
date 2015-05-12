% VECADDD
% This function adds two vectors given in polar form
% All angles in *degrees* (inputs and output), hence the D in VECADDD
% Usage: [M, ph] = vecaddd(M1, ph1, M2, ph2)

function [M, ph] = vecaddd(M1, ph1, M2, ph2)

ph1 = ph1*(pi/180);
ph2 = ph2*(pi/180);

avec = (M1*cos(ph1)+M2*cos(ph2))+j*(M1*sin(ph1)+M2*sin(ph2));
% M = abs(avec); ph = angle(avec);
% M = abs(avec); ph = phase(avec);

M = abs(avec); ph = (180/pi)*atan2(imag(avec),real(avec));