% WRITEINDEXFILE
% This function writes the index file to be used with J. Starmer's
% free-energy signal calculation software
% Usage: writeindexfile(CertData,CERTFILE,HypData,HYPFILE)

function writeindexfile(CertData,CERTFILE,HypData,HYPFILE)

% Certain sequences
fout=fopen(CERTFILE,'w');
for i=1:size(CertData,1)
    start=CertData{i,4}; stop=CertData{i,5};    
    if CertData{i,2}=='+'
        fprintf(fout,'%d..%d\n',start,stop);
    elseif CertData{i,2}=='-'
        fprintf(fout,'complement(%d..%d)\n',start,stop);
    end
end
fclose(fout);

% Hypothetical sequences
fout=fopen(HYPFILE,'w');
for i=1:size(HypData,1)
    start=HypData{i,4}; stop=HypData{i,5};    
    if HypData{i,2}=='+'
        fprintf(fout,'%d..%d\n',start,stop);
    elseif HypData{i,2}=='-'
        fprintf(fout,'complement(%d..%d)\n',start,stop);
    end
end
fclose(fout);
