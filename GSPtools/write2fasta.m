% WRITE2FASTA: Writes a string to a file in FASTA format
% A file in FASTA format starts with a header >.... followed by lines
% containing the sequence, maximum of L characters per line
% USAGE: write2fasta(fname,S,header,L)
% INPUTS:
% fname = name of the FASTA file (with the .fasta extension)
% S = sequence
% header = line in the header (without the >)
% L = number of bases to a line

function write2fasta(fname,S,header,L)

fid=fopen(fname,'w');
if ~fid
    error('Unable to open file to write FASTA sequence\n');
end

H=strcat('>',header); fprintf(fid,H,'%c'); fprintf(fid,'\n');
% Print to file, L characters per line
for i=1:floor(length(S)/L)
    seq=S(L*(i-1)+1:L*(i)); fprintf(fid,seq,'%c'); fprintf(fid,'\n');
end
seq=S(L*(i)+1:end); fprintf(fid,seq,'%c'); fprintf(fid,'\n');
fclose(fid);
