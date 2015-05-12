% CUMM_MAG_PHASE : "Cummulative magnitude and phase"
% This function calculates the cummulative magnitude and phase of the input
% signal
%
% USAGE:
% [Mag, Phase, Err] = cumm_mag_phase(x)
% x is a row vector that contains the free energy values
% Phase is in radians

function [Mag, Phase, Err] = cumm_mag_phase(x)

if rem(length(x),3)~=0
    fprintf('\ncumm_mag_phase: Signal rounded off to a codon multiple');
    x = x(1:(length(x)-rem(length(x),3)));
end

Mag = zeros(1,length(x)/3);
Phase = zeros(1,length(x)/3);
Err = zeros(2,length(x)/3);

for j=1:length(Mag)
    [Mag(1,j), Phase(1,j), Err(1:2,j)] = calc_mag_phase(x(1:(j*3)),0);    
end
