% CHECKWRSIGS: "Check written signals"
% This function checks if the signals written-out to the text-files
% are correct in length. Using the starts and stops listed in the
% information matrix, it first computes the lengths of all the genes. It
% then compares these lengths with the lengths of the signals read-in
% from the text-files. 
% 
% Usage:
% I = checkwrsigs(Data,SIGFILE,sigpreset,sigpostset)
% Inputs:
% Data = information matrix
% SIGFILE = file containing signals, one gene per line
% sigpreset/sigpostset = preset/postset for the signals in SIGFILE
% Output:
% I = indices of signals whose lengths don't match up

function I = checkwrsigs(Data,SIGFILE,sigpreset,sigpostset)

% Calculate lengths of signals from the information matrix
for i=1:size(Data,1)
    L(i,1)=Data{i,5}-Data{i,4}+1;
end

% Calculate lengths of signals from the text-file
fp=fopen(SIGFILE,'r');
i=1;
while 1
    sigline=fgetl(fp); 
    if ~ischar(sigline), break, end
    sigfull=str2num(sigline); signal=sigfull(sigpreset+1:end-sigpostset);
    Lw(i,1)=length(signal);
    i=i+1;
end
fclose(fp);

if length(L)~=length(Lw)
    error('\nThe number of signals written-out is different from the number of genes in the information matrix!');
end

% Find indexes where there are differences
I=find(L~=Lw);
