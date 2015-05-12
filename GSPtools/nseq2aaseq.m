% NSEQ2AASEQ
% Convert nucleotide sequence to amino-acid sequence
% Based on the Reverse Codon Table at: http://en.wikipedia.org/wiki/Codon
% 
% Usage: AA = nseq2aaseq(S,Case,Type,Choice)
% S = nucleotide sequence of specified Case and Type
% Choice=3: three-letter abbreviations, Choice=1: one-letter abbreviations
% AA = cell-array of amino acids

function AA = nseq2aaseq(S,Case,Type,Choice)

if rem(length(S),3)~=0
    error('Nucleotide sequence does not contain an integer number of codons!');
end

C=''; % character-array of codons in the sequence S
for j=1:(length(S)/3)
    C=strvcat(C,S(3*(j-1)+1:3*(j-1)+3));
end

AA=cell(1,length(S)/3); % cell-array of amino acids
for k=1:(length(S)/3)
    AA{1,k}=codon2aacid(C(k,:),Case,Type,Choice);
end

if Choice==1
    AA=char(AA)';
elseif Choice==3
    AA=char(AA);
end