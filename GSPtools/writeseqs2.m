% WRITESEQS2
% This function writes out sequences one gene per line. Use GETSEQSIG1 
% to pull-up a sequence and its signal quickly. 
% Ignores genes having non-permissible characters (such as N, R, etc)
% and updates the information matrices. 
% 
% Usage: [Ncert,Nhyp]=writeseqs2(TempCertData,CERTFILE,TempHypData,HYPFILE,preset,postset,S)
% INPUTS
% TempCertData = (initial) information matrix for the verified sequences
% CERTFILE = name of text file to which verified sequences will be written
% TempHypData = (initial) information matrix for the hypothetical sequences
% HYPFILE = name of text file to which hypothetical sequences will be written
% preset = number of bases before the first base of the start codon
% postset = number of bases after the last base of the stop codon
% S = whole-genome sequence of that organism
% OUTPUTS
% Ncert = number of verified sequences written out
% Nhyp = number of hypothetical sequences written out

function [Ncert,Nhyp]=writeseqs2(TempCertData,CERTFILE,TempHypData,HYPFILE,preset,postset,S)

% Certain sequences
countCert=0; CertInd=[];
fout=fopen(CERTFILE,'w');
for i=1:size(TempCertData,1)
    fprintf('\nProcessing row %d of certains information matrix',i);
    start=TempCertData{i,4}; stop=TempCertData{i,5};    
    if TempCertData{i,2}=='+'
        if ((start-preset)>0)&((stop+postset)<=length(S))             
            seq=S(start-preset:stop+postset);      
            % Make sure all characters in the sequence are legitimate,
            % ignore sequence if any of the characters are non-(a,t,c,g)
            L = length(find(seq=='a'))+length(find(seq=='t'))+length(find(seq=='c'))+length(find(seq=='g'));
            if L==length(seq)
                gene=seq;
                fprintf(fout,'%s',gene); fprintf(fout,'\n');
                countCert=countCert+1; CertInd=[CertInd; i];
            else
                fprintf('\nGene ignored for containing illegitimate nucleotides');
            end
        else % ignore if index is out-of-bounds
            fprintf('\nIndex out-of-bounds for forward-strand verified gene named %s',TempCertData{i,3});
        end
    elseif TempCertData{i,2}=='-'
        if ((start-postset)>0)&((stop+preset)<=length(S)) 
            seq=S(start-postset:stop+preset);                        
            % Make sure all characters in the sequence are legitimate,
            % ignore sequence if any of the characters are non-(a,t,c,g)
            L = length(find(seq=='a'))+length(find(seq=='t'))+length(find(seq=='c'))+length(find(seq=='g'));            
            if L==length(seq)
                gene=revcomp(seq,0,0);
                fprintf(fout,'%s',gene); fprintf(fout,'\n');
                countCert=countCert+1; CertInd=[CertInd; i];
            else
                fprintf('\nGene ignored for containing illegitimate nucleotides');
            end            
        else % ignore if index is out-of-bounds            
            fprintf('\nIndex out-of-bounds for reverse-strand verified gene named %s',TempCertData{i,3});
        end
    end
end
fclose(fout);
Ncert=countCert;
CertData=TempCertData(CertInd,:); save CertDataFile.mat CertData;


% Hypothetical sequences
countHyp=0; HypInd=[];
fout=fopen(HYPFILE,'w');
for i=1:size(TempHypData,1)
    fprintf('\nProcessing row %d of hypotheticals information matrix',i);
    start=TempHypData{i,4}; stop=TempHypData{i,5};    
    if TempHypData{i,2}=='+'
        if ((start-preset)>0)&((stop+postset)<=length(S))             
            seq=S(start-preset:stop+postset);         
            % Make sure all characters in the sequence are legitimate,
            % ignore sequence if any of the characters are non-(a,t,c,g)
            L = length(find(seq=='a'))+length(find(seq=='t'))+length(find(seq=='c'))+length(find(seq=='g'));
            if L==length(seq)
                gene=seq;
                fprintf(fout,'%s',gene); fprintf(fout,'\n');
                countHyp=countHyp+1; HypInd=[HypInd; i];
            else
                fprintf('\nGene ignored for containing illegitimate nucleotides');
            end                                    
        else % ignore if index is out-of-bounds
            fprintf('\nIndex out-of-bounds for forward-strand hypothetical gene named %s',TempHypData{i,3});
        end
    elseif TempHypData{i,2}=='-'
        if ((start-postset)>0)&((stop+preset)<=length(S))             
            seq=S(start-postset:stop+preset);                                                
            % Make sure all characters in the sequence are legitimate,
            % ignore sequence if any of the characters are non-(a,t,c,g)
            L = length(find(seq=='a'))+length(find(seq=='t'))+length(find(seq=='c'))+length(find(seq=='g'));            
            if L==length(seq)
                gene=revcomp(seq,0,0); 
                fprintf(fout,'%s',gene); fprintf(fout,'\n');
                countHyp=countHyp+1; HypInd=[HypInd; i];                                
            else
                fprintf('\nGene ignored for containing illegitimate nucleotides');
            end                                                
        else % ignore if index is out-of-bounds
            fprintf('\nIndex out-of-bounds for reverse-strand hypothetical gene named %s',TempHypData{i,3});
        end
    end
end
fclose(fout);
Nhyp=countHyp;
HypData=TempHypData(HypInd,:); save HypDataFile.mat HypData;
