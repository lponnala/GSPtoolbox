% CODON2AACID
% This function returns the name of the amino acid which the input codon
% corresponds to. The amino acid name will be returned in one-letter or
% three-letter abbreviation, depending on the input 'Choice'. 
% Based on the Reverse Codon Table at: http://en.wikipedia.org/wiki/Codon
% This table is also saved as D:/Research/BioTools/MATLAB/ReverseCodonTable.xls
% 
% USAGE:
% y = codon2aacid(codon, Case, Type, Choice)
% INPUTS:
% codon = character string of specified Case and Type
% Case=1: capital letters, Case=0: small letters
% Type=1: RNA sequence, Type=0: DNA sequence
% Choice=3: three-letter abbreviation, Choice=1: one-letter abbreviation
% 
% Note: Start codon is not annotated separately, it is the same as Met

function y = codon2aacid(codon, Case, Type, Choice)

if length(codon)~=3
    error('Codon must be a triplet!');
end

% Convert codon to upper-case RNA format, if needed
if (Case~=1)|(Type~=1)
    codon= num2char(char2num(codon,Case,Type),1,1);
end

% Reverse Codon Table
C = cell(21,3);
C{1,1}= 'Ala'; C{1,2}= 'A'; C{1,3}= ['GCU'; 'GCC'; 'GCA'; 'GCG'];
C{2,1}= 'Arg'; C{2,2}= 'R'; C{2,3}= ['CGU'; 'CGC'; 'CGA'; 'CGG'; 'AGA'; 'AGG'];
C{3,1}= 'Asn'; C{3,2}= 'N'; C{3,3}= ['AAU'; 'AAC'];
C{4,1}= 'Asp'; C{4,2}= 'D'; C{4,3}= ['GAU'; 'GAC'];
C{5,1}= 'Cys'; C{5,2}= 'C'; C{5,3}= ['UGU'; 'UGC'];
C{6,1}= 'Gln'; C{6,2}= 'Q'; C{6,3}= ['CAA'; 'CAG'];
C{7,1}= 'Glu'; C{7,2}= 'E'; C{7,3}= ['GAA'; 'GAG'];
C{8,1}= 'Gly'; C{8,2}= 'G'; C{8,3}= ['GGU'; 'GGC'; 'GGA'; 'GGG'];
C{9,1}= 'His'; C{9,2}= 'H'; C{9,3}= ['CAU'; 'CAC'];
C{10,1}= 'Ile'; C{10,2}= 'I'; C{10,3}= ['AUU'; 'AUC'; 'AUA'];
C{11,1}= 'Leu'; C{11,2}= 'L'; C{11,3}= ['UUA'; 'UUG'; 'CUU'; 'CUC'; 'CUA'; 'CUG'];
C{12,1}= 'Lys'; C{12,2}= 'K'; C{12,3}= ['AAA'; 'AAG'];
C{13,1}= 'Met'; C{13,2}= 'M'; C{13,3}= ['AUG'];
C{14,1}= 'Phe'; C{14,2}= 'F'; C{14,3}= ['UUU'; 'UUC'];
C{15,1}= 'Pro'; C{15,2}= 'P'; C{15,3}= ['CCU'; 'CCC'; 'CCA'; 'CCG'];
C{16,1}= 'Ser'; C{16,2}= 'S'; C{16,3}= ['UCU'; 'UCC'; 'UCA'; 'UCG'; 'AGU'; 'AGC'];
C{17,1}= 'Thr'; C{17,2}= 'T'; C{17,3}= ['ACU'; 'ACC'; 'ACA'; 'ACG'];
C{18,1}= 'Trp'; C{18,2}= 'W'; C{18,3}= ['UGG'];
C{19,1}= 'Tyr'; C{19,2}= 'Y'; C{19,3}= ['UAU'; 'UAC'];
C{20,1}= 'Val'; C{20,2}= 'V'; C{20,3}= ['GUU'; 'GUC'; 'GUA'; 'GUG'];
C{21,1}= 'Stop'; C{21,2}= '*'; C{21,3}= ['UAG'; 'UGA'; 'UAA']; 
% As per NCBI BLAST-search input format, * = translation stop
% See http://www.ncbi.nlm.nih.gov/blast/html/search.html

y='';
for i=1:size(C,1)
    I=strmatch(codon,C{i,3});
    if ~isempty(I)
        if Choice==3
            y=C{i,1};
        elseif Choice==1
            y=C{i,2};
        end
    end
end


% % ----OLD CODE----
% codon = char2num(codon, Case, Type);
% 
% % Create the coding dictionary
% Dict = zeros(21,18);
% Dict(1,1:6) = [4 4 4 4 4 3]; % Phe
% Dict(2,1:18) = [4 4 1 4 4 2 3 4 4 3 4 3 3 4 1 3 4 2]; % Leu
% Dict(3,1:9) = [1 4 4 1 4 3 1 4 1]; % Iso
% Dict(4,1:3) = [1 4 2]; % Met
% Dict(5,1:12) = [2 4 4 2 4 3 2 4 1 2 4 2]; % Val
% Dict(6,1:18) = [4 3 4 4 3 3 4 3 1 4 3 2 1 2 4 1 2 3]; % Ser
% Dict(7,1:12) = [3 3 4 3 3 3 3 3 1 3 3 2]; % Pro
% Dict(8,1:12) = [1 3 4 1 3 3 1 3 1 1 3 2]; % Thr        
% Dict(9,1:12) = [2 3 4 2 3 3 2 3 1 2 3 2]; % Ala
% Dict(10,1:6) = [4 1 4 4 1 3]; % Tyr
% Dict(11,1:9) = [4 1 1 4 1 2 4 2 1]; % Stop
% Dict(12,1:6) = [3 1 4 3 1 3]; % His
% Dict(13,1:6) = [3 1 1 3 1 2]; % Glu
% Dict(14,1:6) = [1 1 4 1 1 3]; % Asp
% Dict(15,1:6) = [1 1 1 1 1 2]; % Lys
% Dict(16,1:6) = [2 1 4 2 1 3]; % AspA
% Dict(17,1:6) = [2 1 1 2 1 2]; % GluA
% Dict(18,1:6) = [4 2 4 4 2 3]; % Cys
% Dict(19,1:3) = [4 2 2]; % Trp
% Dict(20,1:18) = [3 2 4 3 2 3 3 2 1 3 2 2 1 2 1 1 2 2]; % Arg
% Dict(21,1:12) = [2 2 4 2 2 3 2 2 1 2 2 2]; % Gly
% 
% % Create the list of names
% aacidnames = {'Phe' 'Leu' 'Iso' 'Met' 'Val' 'Ser' 'Pro' 'Thr' 'Ala' 'Tyr' 'Stop' 'His' 'Glu' 'Asp' 'Lys' 'AspA' 'GluA' 'Cys' 'Trp' 'Arg' 'Gly'};
% 
% % Look for codon in the dictionary
% for row = 1:size(Dict,1)
%     for col = 1:3:size(Dict,2)
%         if sum(codon==Dict(row,col:col+2))==3
%             y = aacidnames(row);
%         end
%     end
% end
