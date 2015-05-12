% CALCMPX: "Calculate magnitude, phase, displacement"
% Calculates the cummulative magnitude, phase and the
% displacement vector based on the frameshift model
% The wait-time (number of loops) at the stop codon can be passed as a
% parameter to this function. 
% Starts calculating wait-time from codon 2, since codon 1 is locked into
% the P-site at the beginning of translation. 
% 
% USAGE: [Mag,Phase,x] = calcmpx(seq,signal,phi_sp,Names,TAV,C1,C2,Nstop)
% seq = string containing the RNA sequence in lower case
% signal = row vector containing the free energy signal
% phi_sp = species-specific phase angle, typically the mean of all the
% individual signal phase angles
% Names,TAV = pre-calculated tRNA availability
% C1 = sets the size of the incremental displacement
% C2 = estimate of initial displacement
% Nstop = number of loops at the stop codon
% 
% Note: This function does not check if the signal length is a multiple of 3

function [Mag,Phase,x] = calcmpx(seq,signal,phi_sp,Names,TAV,C1,C2,Nstop)

% Make sure sequence and signal are of the same length!
if length(seq)~=length(signal)
    error('Sequence and signal are NOT of the same length!');
end

% ------------------------------------------------------------
% CALCULATE CUMMULATIVE MAGNITUDE AND PHASE
% ------------------------------------------------------------    
[Mag, Phase, Err] = cumm_mag_phase(signal);
I1 = find(Err(1,:));
if ~isempty(I1)
    fprintf(1,'\nMagnitude negative at %d indices:',length(I1));
end
I2 = find(Err(2,:)); 
if ~isempty(I2)
    fprintf(1,'\nEquations not satisfied at %d indices:',length(I2));
end
fprintf(1,'\nNumber of codons: %d',length(Mag));
numcodons = length(Mag);

% ------------------------------------------------------------
% CALCULATE DIFFERENTIAL VECTORS
% ------------------------------------------------------------    
L = 3; % input('Enter number of points: ');
if rem(L,2)==0
    error('L has to be an odd number');
end
P = 1; % input('Enter order of polynomial: ');
for k=1:numcodons
    index = [min(max(1,k-(L-1)/2),numcodons-L+1):1:min(max(1,k-(L-1)/2),numcodons-L+1)+L-1];
    polyMag = polyfit(1:L,Mag(index),P);
    polyPhase = polyfit(1:L,Phase(index),P);
    
    dA_dc(1,k) = polyval(polyder(polyMag),k);
    dphi_dc(1,k) = polyval(polyder(polyPhase),k);
    
    D = exp(j*Phase(1,k))*(dA_dc(1,k) + j*(Mag(1,k)*dphi_dc(1,k)));
    
    Dvec(k,1) = abs(D); 
    % Dvec(k,1) = sqrt(dA_dc(1,k)^2+(Mag(1,k)*dphi_dc(1,k))^2);
    % Dvec(k,2) = phase(D); 
    % --- for those who don't have the Signal Processing Toolbox ---
    Dvec(k,2) = Phase(1,k) + atan2(Mag(1,k)*dphi_dc(1,k),dA_dc(1,k));  
end

% ------------------------------------------------------------
% CALCULATE DISPLACEMENT 
% ------------------------------------------------------------    
Nloop = []; 
currRF=0; % current reading frame

x = 0; x(1,2) = C2;

for k=2:numcodons-1    
    if abs(x(1,k))<1
        codon=seq(3*(k-1)+1:3*(k-1)+3); Nloop(k)=nloopcalc(codon,0,1,Names,TAV,Nstop);
    elseif x(1,k)<-1
        codon=seq(3*(k-1):3*(k-1)+2); Nloop(k)=nloopcalc(codon,0,1,Names,TAV,Nstop);
    elseif x(1,k)>1
        codon=seq(3*(k-1)+2:3*(k-1)+4); Nloop(k)=nloopcalc(codon,0,1,Names,TAV,Nstop);
    end    
  
%     if k==26
%         fprintf('\n\n---- k = %d, codon = %s ----',k,codon);
%         fprintf('\nDiffVecMag = %f, DiffVecPhase = %f, Nloop(k) = %d',Dvec(k,1),(180/pi)*Dvec(k,2),Nloop(k));
%     end
    
    phi_signal(1,k) = Dvec(k,2); x_temp = x(1,k); 
    for wt=1:Nloop(k)
        phi_dx = ((pi/3)*x_temp)-phi_sp;
        dx = -C1*Dvec(k,1)*sin(phi_signal(1,k) + phi_dx); % Correct
%         if k==26
%             fprintf('\nArg of sin() = %f',(180/pi)*(phi_signal(1,k) + phi_dx));
%         end
        x_temp = dx + x_temp;
    end
    x(1,k+1) = x_temp;             
end
