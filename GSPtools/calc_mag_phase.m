% CALC_MAG_PHASE : "Calculate magnitude, phase"
% This function calculates the magnitude and phase of a given length of
% signal using sine-wave interpolation
% The length of the signal must be a codon multiple
% 
% INPUTS: 
% x = signal upto a certain number of codons
% avg_choice = (1 : energy values will be averaged during magnitude calculation) (any other value : no averaging)
% 
% USAGE:
% [A,theta,Err] = calc_mag_phase(x, avg_choice);
% theta = phase angle in radians
% Err(1)=1: Magnitude negative, Err(2)=2: Equations not satisfied

function [A,theta,Err] = calc_mag_phase(x, avg_choice);

% Initialization
A = 0; theta = 0; Err = [0 0]';

if rem(length(x),3)~=0
    display('Signal rounded off to a codon multiple');
    x = x(1:(length(x)-rem(length(x),3)));
end

% Initialize the values of the memory registers
M = zeros(1,3);      

% Sum the third values to fill up the memory registers                        
for j = 1:3:length(x)
    for k = 1:1:3
        M(1,k) = M(1,k) + x(1,j+k-1); 
    end
end    

% The averaging calculation works only for the case when length(x) is a multiple of 3    
if avg_choice==1
    M = (M-mean(M))/(length(x)/3);            									
else
    M = (M-mean(M)); 
end

% Calculate magnitude and phase using formulae
if (M(1) ~= 0)
    theta = atan2( M(1)*sqrt(3),(M(1) + 2*M(2)) ); 
    A = M(1)/sin(theta);
else
    if((M(3)+M(2))~=0)
        theta = atan2( (M(3)+M(2))*sqrt(3),(M(3)-M(2)) );
        A = -(M(3)+M(2))/sin(theta); 
    end
end

% Check if magnitude is negative
if A<0
    Err(1)=1; % fprintf(1,'Magnitude negative!! A = %f\n',A);    
end

% Check if solution satisfies all 3 equations
eps = 0.0001;
if abs(M(1)-(A*sin(theta)))>eps
    Err(2)=2.1; % fprintf(1,'Equation 1 not satisfied\n');
end
if abs(M(2)-(A*sin(theta+(2*pi/3))))>eps
    Err(2)=2.2; % fprintf(1,'Equation 2 not satisfied\n');
end
if abs(M(3)-(A*sin(theta+(4*pi/3))))>eps
    Err(2)=2.3; % fprintf(1,'Equation 3 not satisfied\n');
end
