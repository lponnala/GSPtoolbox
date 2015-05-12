% GETSEQSIG1 : "Get sequence and signal"
% Sequence and signal of the specified gene are extracted from the
% text-files supplied. Sequence is returned in lower-case RNA format, ready
% for use with the frameshift-model program. 
% Use this function repeatedly in a loop to extract sequences and signals
% for a bunch of genes. It is best to extract the sequence and signal for
% one gene at a time. Here's why:
% If multiple sequences and signals need to be extracted with just one
% function call, then the function will need to return cell arrays of the
% sequences and signals - it is better to avoid this complication!
% 
% USAGE:
% [Seq,Signal] = getseqsig1(genename,Data,SEQFILE,seqpreset,seqpostset,SIGFILE,sigpreset,sigpostset)
% Inputs:
% genename = string containing the name of the gene, e.g. 'prfB'
% Data = information matrix
% SEQFILE = text-file containing sequences, one gene per line
% (seqpreset/seqpostset) = sequence preset/postset
% SIGFILE = text-file containing signals, one gene per line
% (sigpreset/sigpostset) = signal preset/postset
% Outputs:
% Seq = string in lower-case RNA format
% Signal = vector containing free energies
% Note:
% SEQFILE is assumed to contain sequences in lower-case DNA format. 

function [Seq,Signal] = getseqsig1(genename,Data,SEQFILE,seqpreset,seqpostset,SIGFILE,sigpreset,sigpostset)

% Get index of the specified gene in the information matrix
[I,err] = findgene(genename,Data,0);

if err==0
    fseq=fopen(SEQFILE,'r'); if ~fseq, error('Unable to open sequence-file!'), end
    fsig=fopen(SIGFILE,'r'); if ~fseq, error('Unable to open signal-file!'), end    
    linenum=0;
    while 1
        seqline = fgetl(fseq); 
        sigline = fgetl(fsig);
        if ~ischar(seqline), break, end
        linenum=linenum+1;        
        if linenum==I        
            % Extract sequence
            Seq=seqline(seqpreset+1:end-seqpostset);
            Seq=num2char(char2num(Seq,0,0),0,1); 
            % Extract signal
            sigfull=str2num(sigline); 
            Signal=sigfull(sigpreset+1:end-sigpostset);
        end
    end    
    fclose(fseq); fclose(fsig);
elseif err==1
    error('Multiple matches found for the specified gene %s!',genename);
elseif err==2
    error('Specified gene %s NOT FOUND!',genename);
end
