% LISTCHARBYCODON
% This function lists a character sequence by codon and name of the
% amino acid for which the triplet codes
% (Case = 1 : capital letters ; Case = 0 : small letters)
% (Type = 1 : RNA sequence ; Type = 0 : DNA sequence)
% (Display = 1: display output; Display = 0: hide output)
% 
% USAGE: y = listcharbycodon(x, Case, Type, Display)
% To list a small-character RNA sequence by codon and amino acid name, use
% listcharbycodon(sequence, 0, 1);

function y = listcharbycodon(x, Case, Type, Display)

if rem(length(x),3)~=0
    error('Signal length is not of codon multiple!');
end

numcodons = length(x)/3;

for i=1:numcodons    
    codon = x(3*(i-1)+1:3*(i-1)+3);
    aacidname = codon2aacid(codon, Case, Type, 3); 
    S{i} = strcat(int2str(i),'_______________',codon,'______',aacidname);
end

y = S;

if Display==1
    disp('-------------------------------------------------------');
    disp('CodonNumber    Codon    Name');
    for i=1:numcodons
        disp(y{i});
    end
    disp('-------------------------------------------------------');
end