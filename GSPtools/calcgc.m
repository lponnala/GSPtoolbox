% CALCGC : "Calculate (G+C) content"
% This function calculates the (G+C) content of genes
% Usage: GC = calcgc(Data,S)
% Data = information matrix 
% S = whole-genome sequence of the organism

function GC = calcgc(Data,S)

for i=1:size(Data,1)
    start=Data{i,4}; stop=Data{i,5};    
    gene=S(start:stop);         
    GC(i,1) = length(find((gene=='g')|(gene=='c')))/length(gene);
end
