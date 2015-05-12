% GETTAILS
% This function extracts the last L bases of the 16S rRNA sequences and
% returns them as a character array in 5'-3' format. The function string as
% annotated in the structural RNA table needs to be provided as input. 
% USAGE:
% Tails = gettails(func,SRNATABLE,S,L)
% INPUTS:
% func = string containing the function
% SRNATABLE = Excel95 table containing the structural RNAs
% S = whole-genome sequence of that organism
% L = length of tail
% OUTPUT:
% Tails = char array containing tails in 5'-3' format

function Tails = gettails(func,SRNATABLE,S,L)

[A, B] = XLSREAD(SRNATABLE); 
AStartCol = 1; AStopCol = 2; % From XLSREAD, check if these are correct
BFuncCol = 1; BStrandCol = 4; % From XLSREAD, check if these are correct

if size(A,1)~=size(B,1)
    error('A,B matrices NOT of same size');
end

Tails=[];
for i=1:size(B,1)
    funcstr = B{i,BFuncCol}; % extract function
    if strcmp(funcstr,func)
        if B{i,BStrandCol}=='+'
            Seq=S(A(i,AStartCol):A(i,AStopCol));
            Tails=strvcat(Tails,Seq(end-L+1:end));
        elseif B{i,BStrandCol}=='-'
            Seq=revcomp(S(A(i,AStartCol):A(i,AStopCol)),0,0);
            Tails=strvcat(Tails,Seq(end-L+1:end));
        end
    end
end
