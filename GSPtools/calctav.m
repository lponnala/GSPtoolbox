% CALCTAV : "Calculate tRNA availability"
% This function calculates the codon abundance ratios and the availability
% of tRNAs in a species. 
% The availability of the species tRNAs, is calculated by taking
% into consideration that several codons can have the same tRNA. The codon
% abundance ratios are lumped together, based on the table online at:
% http://www.sci.sdsu.edu/~smaloy/MicrobialGenetics/topics/rev-sup/wobble.html
% 
% USAGE:
% calctav(SEQFILE,Case,Type,seqpreset,seqpostset,NAMESFILE,CARFILE,TAVFILE,TRNAFILE)
% SEQFILE = text file containing sequences, one gene per line
% Case,Type = case, type of sequences in the text-file
% seqpreset = number of bases before the start, common for all sequences
% NAMESFILE = .mat file to which the codon triplets (RNA) will be written
% CARFILE = .mat file to which the frequencies (abundance ratios) of codons will be written
% TAVFILE = .mat file to which the estimated tRNA concentrations (see above) will be written
% TRNAFILE = text file to which the codon-name and tRNA availability will be written

function calctav(SEQFILE,Case,Type,seqpreset,seqpostset,NAMESFILE,CARFILE,TAVFILE,TRNAFILE)

% ------Calculate CAR and save------
[Counts,Names] = getcp(SEQFILE,Case,Type,seqpreset,seqpostset);
CAR = Counts/sum(Counts); % normalize the counts across all codons
save(NAMESFILE,'Names'); save(CARFILE,'CAR');

% ------Calculate TAV and save------
TAV = zeros(size(CAR));

TAV(31)=CAR(31)+CAR(32); TAV(32)=CAR(31)+CAR(32); 
TAV(33)=CAR(33)+CAR(34); TAV(34)=CAR(33)+CAR(34); 
TAV(59)=CAR(59)+CAR(60); TAV(60)=CAR(59)+CAR(60); 
TAV(57)=CAR(57)+CAR(58); TAV(58)=CAR(57)+CAR(58); 
TAV(55)=CAR(55)+CAR(56); TAV(56)=CAR(55)+CAR(56); 
TAV(48)=CAR(48)+CAR(49); TAV(49)=CAR(48)+CAR(49); 
TAV(44)=CAR(44)+CAR(45); TAV(45)=CAR(44)+CAR(45); 
TAV(52)=CAR(52)+CAR(53); TAV(53)=CAR(52)+CAR(53); 
TAV(50)=CAR(50)+CAR(51); TAV(51)=CAR(50)+CAR(51); 
TAV(42)=CAR(42)+CAR(43); TAV(43)=CAR(42)+CAR(43); 
TAV(63)=CAR(63)+CAR(64); TAV(64)=CAR(63)+CAR(64); 
TAV(61)=CAR(61)+CAR(62); TAV(62)=CAR(61)+CAR(62); 
TAV(40)=CAR(40)+CAR(41); TAV(41)=CAR(40)+CAR(41); 
TAV(11)=CAR(11);
TAV(9)=CAR(9)+CAR(10); TAV(10)=CAR(9)+CAR(10); 
TAV(3)=CAR(3)+CAR(4); TAV(4)=CAR(3)+CAR(4); 
TAV(7)=CAR(7)+CAR(8); TAV(8)=CAR(7)+CAR(8); 
TAV(5)=CAR(5)+CAR(6); TAV(6)=CAR(5)+CAR(6); 
TAV(46)=CAR(46)+CAR(47); TAV(47)=CAR(46)+CAR(47); 
TAV(12)=CAR(12);
TAV(1)=CAR(1)+CAR(2); TAV(2)=CAR(1)+CAR(2); 
TAV(25)=CAR(25)+CAR(26); TAV(26)=CAR(25)+CAR(26); 
TAV(23)=CAR(23)+CAR(24); TAV(24)=CAR(23)+CAR(24); 
TAV(21)=CAR(21)+CAR(22); TAV(22)=CAR(21)+CAR(22); 
TAV(19)=CAR(19)+CAR(20); TAV(20)=CAR(19)+CAR(20); 
TAV(17)=CAR(17)+CAR(18); TAV(18)=CAR(17)+CAR(18); 
TAV(29)=CAR(29)+CAR(30); TAV(30)=CAR(29)+CAR(30); 
TAV(27)=CAR(27);
TAV(54)=CAR(54);
TAV(35)=CAR(35)+CAR(36); TAV(36)=CAR(35)+CAR(36); 
TAV(15)=CAR(15)+CAR(16); TAV(16)=CAR(15)+CAR(16); 
TAV(13)=CAR(13)+CAR(14); TAV(14)=CAR(13)+CAR(14);
% Special cases
TAV(28)=CAR(28); % Thr codon has no anticodon specified

save(TAVFILE,'TAV');

% Write out file containing estimated tRNA concentrations
fid=fopen(TRNAFILE,'w');
for i=1:length(Names)
    fprintf(fid,'%s\t%f\n',Names{i},TAV(i));
end
fclose(fid);
