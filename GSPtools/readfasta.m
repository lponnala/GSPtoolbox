% READFASTA: Reads sequence from a file in FASTA format, and returns the
% sequence as a string
% INPUTS:
% fname = string containing the full filename (including .fasta extension)
% OUTPUTS:
% S = sequence string

function S = readfasta(fname)

fid=fopen(fname,'r');
if ~fid
    error('Could not locate FASTA file');
end

S = ''; % create empty string
while 1
    tline = fgetl(fid);
    if ~ischar(tline), break, end    
    % If the first character on the line is a '>', ignore the line
    if tline(1)=='>', continue, end    
    S = strcat(S,tline);    
end
fclose(fid);