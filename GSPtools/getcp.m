% GETCP : "Get codon-bias profile"
% This function gets the codon bias profile from a text-file containing
% sequences, one gene per line
% Returns the counts and codon-triplets (in lower-case RNA format)
% 
% Usage: [Counts,Names] = getcp(SEQFILE,Case,Type,seqpreset,seqpostset)
% Inputs:
% SEQFILE = text-file containing sequences, one gene per line
% seqpreset/seqpostset = preset/postset in the sequences
% Outputs:
% Counts = 64-vector containing the number of times each codon appears
% Names = 64-cell array containing codon triplets in lower-case RNA format
% Note:
% There is a one-to-one match between Counts and Names

function [Counts,Names] = getcp(SEQFILE,Case,Type,seqpreset,seqpostset)

FSind=[]; % indices of sequences NOT having an integer-number of codons
Codons=''; % character array containing all codons in the text-file
linenum=0;

fp=fopen(SEQFILE,'r');
if ~fp
    error('Unable to open SEQFILE!');
end

while 1
    seqline=fgetl(fp);
    if ~ischar(seqline), break, end
    linenum=linenum+1; fprintf('\nProcessing line: %d',linenum);
    S=seqline(seqpreset+1:end-seqpostset);
    if rem(length(S),3)~=0 % Skip sequence
        FSind=[FSind; linenum];
    else 
        % First, make sure the sequence is in lower-case RNA format
        if (Case~=0)|(Type~=1)
            S=num2char(char2num(S,Case,Type),0,1);
        end        
        % Chop the sequence into codons (triplets)        
        C=''; % codons in the current sequence
        for j=1:(length(S)/3)
            C=strvcat(C,S(3*(j-1)+1:3*(j-1)+3));
        end
    end
    Codons=strvcat(Codons,C); % append codons from the current sequence to the overall codon list    
end
fclose(fp);
fprintf('\nNumber of sequences NOT having an integer number of codons = %d',length(FSind));
fprintf('\nThese sequences have been excluded from the calculation of codon-counts');

Counts = zeros(1,64);
Names = {'uuu' 'uuc' 'uua' 'uug' 'cuu' 'cuc' 'cua' 'cug' 'auu' 'auc' 'aua' 'aug' 'guu' 'guc' 'gua' 'gug' 'ucu' 'ucc' 'uca' 'ucg' 'agu' 'agc' 'ccu' 'ccc' 'cca' 'ccg' 'acu' 'acc' 'aca' 'acg' 'gcu' 'gcc' 'gca' 'gcg' 'uau' 'uac' 'uaa' 'uag' 'uga' 'cau' 'cac' 'caa' 'cag' 'aau' 'aac' 'aaa' 'aag' 'gau' 'gac' 'gaa' 'gag' 'ugu' 'ugc' 'ugg' 'cgu' 'cgc' 'cga' 'cgg' 'aga' 'agg' 'ggu' 'ggc' 'gga' 'ggg'};
for i=1:length(Names)
    % fprintf('\nProcessing codon %d',i);    
    Counts(i)=length(strmatch(Names{i},Codons,'exact'));
end
