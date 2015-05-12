% FINDGENES
% This function looks in an information matrix to find genes supplied in a
% text-file (one gene-name to a line). It returns a vector of indices of
% genes found in the information matrix. 
% 
% Usage: I = findgenes(FILE,Data)
% FILE = text-file containing one gene-name per line (could be blank lines
% too, as would happen if the genes were ranked and the line-numbers were
% given to ranks, some ranks may be missing... see LinkGenes.txt)
% Data = information matrix

function I = findgenes(FILE,Data)

% Extract names of all genes in the information matrix
Names=[];
for i=1:size(Data,1)
    Names=strvcat(Names,Data{i,3});
end

% Read-in the contents of the text-file
Genes=textread(FILE,'%s')'; 
I=[];
% Find genes one by one
for i=1:length(Genes)
    genename=Genes{i}; 
    K = strmatch(genename,Names,'exact'); % look for exact match
    if length(K)==1
        I=[I; K]; 
    end       
end
fprintf('\n\nOut of %d genes in the text-file, %d were found in the information matrix',length(Genes),length(I));
