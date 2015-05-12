% CHECKWRSEQS: "Check written sequences"
% This function checks if the sequences written-out to the text-files
% are correct in length. Using the starts and stops listed in the
% information matrix, it first computes the lengths of all the genes. It
% then compares these lengths with the lengths of the sequences read-in
% from the text-files. Also checks if the number of sequences annotated in
% the information matrix is equal to the number of sequences in the
% text-files. 
% 
% Usage:
% I = checkwrseqs(Data,SEQFILE,seqpreset,seqpostset)
% Inputs:
% Data = information matrix
% SEQFILE = file containing sequences, one gene per line
% seqpreset/seqpostset = preset/postset for the sequences in SEQFILE
% Output:
% I = indices of sequences whose lengths don't match up

function I = checkwrseqs(Data,SEQFILE,seqpreset,seqpostset)

% Calculate lengths of sequences from the information matrix
for i=1:size(Data,1)
    L(i,1)=Data{i,5}-Data{i,4}+1;
end

% Calculate lengths of sequences from the text-file
fp=fopen(SEQFILE,'r');
i=1;
while 1
    seqline=fgetl(fp);
    if ~ischar(seqline), break, end
    seq=seqline(seqpreset+1:end-seqpostset);
    Lw(i,1)=length(seq);
    i=i+1;
end
fclose(fp);

if length(L)~=length(Lw)
    error('\nThe number of sequences written-out is different from the number of genes in the information matrix!');
end

% Find indexes where there are differences
I=find(L~=Lw);
