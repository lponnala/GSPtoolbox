% PICKGG: Helps pick good genes based on the periodogram. The logic
% is that a sine-wave can be fit only if there is significant activity at
% f=1/3. Thus, this program picks out good genes based on signal criteria. 
% Procedure:
% 1. Check if signal length is equal to the sequence length (1 signal point for every base)
% 2. Check if signal length is a multiple of 3 (ignore if not: frameshifter!)
% 3. Check if signal length is greater than the minimum length
% 4. Do periodogram test - check if strong 1/3 harmonic
% Good gene = gene that has no frameshifts, is longer than the specified
% length and passes the detection test
% 
% USAGE: 
% GG = pickgg(Data,SIGFILE,preset,postset,minL,alpha,display)
% INPUTS:
% preset, postset refer to the signals
% minL = signals greater than or equal to this length are chosen for analysis
% alpha = significance level of periodogram test
% (display=1): plot periodogram
% OUTPUTS:
% GG = [linenum, length(signal), mu, A, theta, SNR]
% 
% Calls other functions: DETECTF, EST_PAR

function GG = pickgg(Data,SIGFILE,preset,postset,minL,alpha,display)

% Input files
fid = fopen(SIGFILE); % File containing the signals
if ~fid
    fprintf('\nUnable to open file containing signals!\n');
end

MAX_SEQ_LEN = 15000; % maximum sequence length
linenum = 0; % line counter for signal-file
selectcount = 0; % counts selected signals

% Create an empty matrix to store indices and parameters of good genes
GG = []; % [linenum,length(signal),mu,A,theta,SNR]
sgf = []; % selected genes that failed
FS = []; % genes for which length is not a multiple of 3, i.e. frameshifters

while 1
    tline = fgetl(fid);
    if ~ischar(tline), break, end
    
    linenum=linenum+1;    
    fprintf('\nProcessing line %d of %d',linenum,size(Data,1)); 
    signalfull = str2num(tline);
    signal = signalfull(preset+1:end-postset); 
    
    % Make sure signal length = sequence length
    seqlen = Data{linenum,5}-Data{linenum,4}+1;
    if (length(signal)~=seqlen)
        error('  Signal length NOT equal to Sequence length');
    end
    
    % Ignore frameshift genes
    % if rem(length(signal),3)~=0, continue, end; 
    if rem(length(signal),3)~=0
        FS = [FS; linenum]; continue;
    end
    
    if length(signal)>=minL
        selectcount = selectcount+1;
        % fprintf('\nSelect count = %d \t Gene index = %d',selectcount,linenum);
        % Do periodogram test
        H = detectf(signal,alpha,display);
        if H==1 % freq=1/3 component exists
            [mu,A,theta,SNR,err] = est_par(signal); % Calculate parameters
            if err==0
                GG = [GG; [linenum, length(signal), mu, A, theta, SNR]];
            end
        else
            sgf = [sgf; linenum]; % store in the list of "selected but failed" genes
        end        
    end
end
fclose(fid);

fprintf('\n\nNumber of input genes = %d',size(Data,1));
fprintf('\nNumber of genes selected = %d, passed = %d, failed = %d',selectcount,size(GG,1),length(sgf));
fprintf('\nNumber of genes whose length is not a codon-multiple = %d.',length(FS));
fprintf('\tTheir names are: \n'); 
for i=1:length(FS)
    fprintf('%s\n',Data{FS(i),3});
end
% dispopt = input('\nDo you want to see indices of selected genes that failed the detection test? (1=yes,0=no): ');
% if dispopt==1
%     fprintf('\nLinennumbers of genes that were selected but failed the detection test: \n');
%     disp(sgf);
% end