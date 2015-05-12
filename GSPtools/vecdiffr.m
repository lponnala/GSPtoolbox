% VECDIFFR
% This function calculates the difference of two vectors given in polar
% form
% All angles in *radians* (inputs and output), hence the R in VECDIFFR
% Usage: [M, ph] = vecdiffr(M1, ph1, M2, ph2)

function [M, ph] = vecdiffr(M1, ph1, M2, ph2)

dvec = (M1*cos(ph1)-M2*cos(ph2))+j*(M1*sin(ph1)-M2*sin(ph2));
% M = abs(dvec); ph = angle(dvec);
% M = abs(dvec) ph = phase(dvec);
M = abs(dvec); ph = atan2(imag(dvec),real(dvec));
