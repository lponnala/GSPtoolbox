% VECDIFFD
% This function calculates the difference of two vectors given in polar
% form
% All angles in *degrees* (inputs and output), hence the D in VECDIFFD
% Usage: [M, ph] = vecdiffd(M1, ph1, M2, ph2)

function [M, ph] = vecdiffd(M1, ph1, M2, ph2)

ph1 = ph1*(pi/180);
ph2 = ph2*(pi/180);

dvec = (M1*cos(ph1)-M2*cos(ph2))+j*(M1*sin(ph1)-M2*sin(ph2));
% M = abs(dvec); ph = angle(dvec);
% M = abs(dvec) ph = phase(dvec); 

M = abs(dvec); ph = (180/pi)*atan2(imag(dvec),real(dvec));
