% REVCOMP
% This function returns the reverse-complement of a character sequence. The
% reverse-complement will be of the same case and type as the input
% sequence. 
% 
% USAGE: y = revcomp(seq,Case,Type)
% INPUT:
% seq = sequence of characters
% Case=0: lower-case, Case=1: upper-case
% Type=0: DNA, Type=1: RNA

function y = revcomp(seq,Case,Type)

numseq=char2num(seq,Case,Type); % original (input) sequence in numeric form

numcseq=zeros(size(numseq)); % complement sequence in numeric form
numcseq(find(numseq==1))=4;
numcseq(find(numseq==4))=1;
numcseq(find(numseq==2))=3;
numcseq(find(numseq==3))=2;

numrcseq=fliplr(numcseq); % reverse-complement sequence in numeric form
y=num2char(numrcseq,Case,Type); % reverse-complement sequence in character form
