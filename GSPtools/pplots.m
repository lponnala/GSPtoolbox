% This function produces polar plots of specified signals and returns the
% final cummulative magnitude and phase for each gene
% USAGE: 
% [M,theta] = pplots(SIGFILE,I,preset,postset,display)
% INPUTS:
% SIGFILE = text-file containing the signals, one per line
% I = vector containing the indices of signals to be examined
% preset, postset refer to the signals
% (display=1):polar plots will be shown, (display=0):polar plots will not be shown
% OUTPUTS:
% [M,theta] = set of final values of cummulative magnitude and phase (angle in radians)
% NOTES:
% This function does not check if the signal length is a multiple of 3

function [M,Theta] = pplots(SIGFILE,I,preset,postset,display)

% Input files
fid = fopen(SIGFILE); % File containing the signals
if ~fid
    fprintf('\nUnable to open file containing signals!\n');
end

linenum = 0; % line counter for signal-file
index=1;

while 1
    tline = fgetl(fid);
    if ~ischar(tline), break, end
    linenum=linenum+1;    
    if linenum==I(index)        
        fprintf('\nProcessing signal %d of %d',index,length(I));
        signalfull = str2num(tline);
        signal = signalfull(preset+1:end-postset); 
        [mag,phase] = cumm_mag_phase(signal); M(index,1)=mag(end); Theta(index,1)=phase(end);
        if display==1
            polar(phase,mag); title(strcat('Polar plot',' :',int2str(index),'of ',int2str(length(I))));
            pause;
        end
        if index<length(I), index=index+1; else, return, end
    end    
end
fclose(fid);