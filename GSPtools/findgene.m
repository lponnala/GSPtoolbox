% FINDGENE
% This function searches the information matrix and returns the index (row
% number) of a specified gene. All other information pertaining to the 
% gene may be extracted from that row of the information matrix. 
% 
% Usage: [y,err] = findgene(genename,Data,display)
% y = matching index, i.e. row number
% (err=0): found uniquely, (err=1): multiple matches, (err=2): not found
% genename = string containing the gene name
% Data = information matrix
% (display=1): display information; (display=0): suppress display

function [y,err] = findgene(genename,Data,display)

% Extract all gene-names
Names=[];
for i=1:size(Data,1)
    Names=strvcat(Names,Data{i,3});
end

% Search for matching index
I = strmatch(genename,Names,'exact'); % look for exact match
if length(I)==1
    y = I; err = 0;
    if display==1
        fprintf('\nName = %s,\tStrand = %s,\tStart = %d,\tStop = %d,\nFunction = %s\n',Data{I,3},Data{I,2},Data{I,4},Data{I,5},Data{I,1});
    end
elseif length(I)>1
    y = I; err = 1;
    if display==1
        fprintf('\nMultiple matches among Certains!\n');
    end
elseif length(I)==0 
    y = []; err = 2;
    if display==1
        fprintf('\nNOT FOUND!!\n');
    end
end