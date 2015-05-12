% GETANNO
% This function extracts the annotation from the protein coding table and
% returns the information matrices for certain and hypothetical sequences
% Usage: [CertData,HypData] = getanno(PCTABLE)
% Inputs:
% PCTABLE = Excel95 file containing the protein coding table
% Outputs:
% Datafile format: [Function,Strand,Name,Start,Stop]

function [CertData,HypData] = getanno(PCTABLE)

[A, B] = XLSREAD(PCTABLE); 
BFuncCol = 1; BStrandCol = 4; BNameCol = 8; % From XLSREAD, check if these are correct
AStartCol = 1; AStopCol = 2; % From XLSREAD, check if these are correct

% Separately store information for certains and hypotheticals
if size(A,1)~=size(B,1)
    error('A,B matrices NOT of same size');
end

HypCount=1; CertCount=1;
for i=1:size(A,1)
    HypFlag=0; % set hyp flag to zero
    funcstr = B{i,BFuncCol}; % extract function
    % Test if hypothetical
    if (~isempty(regexp(funcstr,'hypothetical')))|(~isempty(regexp(funcstr,'putative')))|(~isempty(regexp(funcstr,'possible')))|(~isempty(regexp(funcstr,'probable')))
        HypFlag=1;
    end
    
    if HypFlag==1
        HypData{HypCount,1}=B{i,BFuncCol}; HypData{HypCount,2}=B{i,BStrandCol}; HypData{HypCount,3}=B{i,BNameCol}; 
        HypData{HypCount,4}=A(i,AStartCol); HypData{HypCount,5}=A(i,AStopCol); 
        HypCount=HypCount+1;
    elseif HypFlag==0
        CertData{CertCount,1}=B{i,BFuncCol}; CertData{CertCount,2}=B{i,BStrandCol}; CertData{CertCount,3}=B{i,BNameCol}; 
        CertData{CertCount,4}=A(i,AStartCol); CertData{CertCount,5}=A(i,AStopCol); 
        CertCount=CertCount+1;
    end
end
if (HypCount-1)~=size(HypData,1)
    error('HypCount index does not match HypData size');
end
if (CertCount-1)~=size(CertData,1)
    error('CertCount index does not match CertData size');
end
if (size(HypData,1)+size(CertData,1))~=size(A,1)
    error('Sizes of HypData and CertData dont add up!');
end
