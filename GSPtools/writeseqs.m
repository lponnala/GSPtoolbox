% WRITESEQS
% This function writes out sequences one gene per line. Use GETSEQSIG1 
% to pull-up a sequence and its signal quickly. 
% 
% Usage: [Ncert,Nhyp]=writeseqs(CertData,CERTFILE,HypData,HYPFILE,preset,postset,S)
% INPUTS
% CertData = information matrix for the verified sequences
% CERTFILE = name of text file to which verified sequences will be written
% HypData = information matrix for the hypothetical sequences
% HYPFILE = name of text file to which hypothetical sequences will be written
% preset = number of bases before the first base of the start codon
% postset = number of bases after the last base of the stop codon
% S = whole-genome sequence of that organism
% OUTPUTS
% Ncert = number of verified sequences written out
% Nhyp = number of hypothetical sequences written out

function [Ncert,Nhyp]=writeseqs(CertData,CERTFILE,HypData,HYPFILE,preset,postset,S)

% Certain sequences
countCert=0;
fout=fopen(CERTFILE,'w');
for i=1:size(CertData,1)
    fprintf('\nProcessing line %d of certains',i);
    start=CertData{i,4}; stop=CertData{i,5};    
    if CertData{i,2}=='+'
        if ((start-preset)>0)&((stop+postset)<=length(S))             
            gene=S(start-preset:stop+postset);         
            fprintf(fout,'%s',gene); fprintf(fout,'\n');
            countCert=countCert+1;
        else % ignore if index is out-of-bounds
            fprintf('\nIndex out-of-bounds for forward-strand verified gene named %s',CertData{i,3});
        end
    elseif CertData{i,2}=='-'
        if ((start-postset)>0)&((stop+preset)<=length(S)) 
            gene=revcomp(S(start-postset:stop+preset),0,0);            
            fprintf(fout,'%s',gene); fprintf(fout,'\n');
            countCert=countCert+1;
        else % ignore if index is out-of-bounds            
            fprintf('\nIndex out-of-bounds for reverse-strand verified gene named %s',CertData{i,3});
        end
    end
end
fclose(fout);
Ncert=countCert;

% Hypothetical sequences
countHyp=0;
fout=fopen(HYPFILE,'w');
for i=1:size(HypData,1)
    fprintf('\nProcessing line %d of hypotheticals',i);
    start=HypData{i,4}; stop=HypData{i,5};    
    if HypData{i,2}=='+'
        if ((start-preset)>0)&((stop+postset)<=length(S))             
            gene=S(start-preset:stop+postset);         
            fprintf(fout,'%s',gene); fprintf(fout,'\n');
            countHyp=countHyp+1;
        else % ignore if index is out-of-bounds
            fprintf('\nIndex out-of-bounds for forward-strand hypothetical gene named %s',HypData{i,3});
        end
    elseif HypData{i,2}=='-'
        if ((start-postset)>0)&((stop+preset)<=length(S))             
            gene=revcomp(S(start-postset:stop+preset),0,0);                        
            fprintf(fout,'%s',gene); fprintf(fout,'\n');
            countHyp=countHyp+1;
        else % ignore if index is out-of-bounds
            fprintf('\nIndex out-of-bounds for reverse-strand hypothetical gene named %s',HypData{i,3});
        end
    end
end
fclose(fout);
Nhyp=countHyp;
