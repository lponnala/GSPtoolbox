% SEQDIFF2
% This function lines-up two character sequences against each other 
% and finds the locations where they differ. Results are written to the
% screen and depicted in a stem-plot. 
% NOTES:
% 1. Make sure both sequences are of the same type (DNA/RNA), 
%    to make the results meaningful in the genetic context! 
% 2. The sequences can be of different cases (upper or lower), they are 
%    converted to a uniform case upon entry into the function. 
% 
% USAGE: y=seqdiff2(seq1,seq2)
% y=0 if seq1=seq2, y=1 otherwise
% seq1 = char string
% seq2 = char string
% 
% Possible extensions:
% 1. Throw an error if both sequences are NOT of the same type (DNA or RNA)
% 2. Comparison over codons (triplets of characters)
% 3. Compare more than 2 sequences

function y=seqdiff2(seq1,seq2)

% Convert both sequences to UPPER case
seq1=upper(seq1); seq2=upper(seq2);

% Make sure both sequences are of same type
% ---- PUT OFF FOR LATER ----

% Compare and print/display results
if length(seq1)==length(seq2)
    % Case: Two Sequences of equal length
    fprintf('\nThe two sequences are of EQUAL length');
    L = length(seq1); D = zeros(1,L);
    I = find(seq1~=seq2); D(I) = 1;
    if isempty(I)
        y=0;
        fprintf(1,'\nThe two sequences are IDENTICAL');
    else
        y=1;
        fprintf('\nThe two sequences DIFFER at %d indices', length(I)); 
        choice=input('\nDo you want to see the locations? (0=no,1=yes): ');
        if choice==1
            disp(I)
            stem(1:L,D); title('Diff plot');
            xlabel('Nucleotide position'); ylabel('');
            AXIS([1 L 0 2])
            set(gca,'XTick',1:1:L); 
            set(gca,'YTick',0:1:3); set(gca,'YTickLabel','');
        end
    end
else
    % Case: Two Sequences of unequal length    
    fprintf('\nThe two sequences are NOT of EQUAL length');
    if input('\nCompare over length of shorter sequence? (0=no,1=yes): ')==1
        L = min(length(seq1),length(seq2)); D = zeros(1,L);    
        fprintf('Comparing over the length of the shorter sequence = %d', L);
        I = find(seq1(1:L)~=seq2(1:L)); D(I) = 1;
        if isempty(I)
            y=0;
            fprintf(1,'\nThe two sequences are IDENTICAL');
        else
            y=1;
            fprintf('\nThe two sequences DIFFER at %d indices', length(I));
            choice=input('\nDo you want to see the locations? (0=no,1=yes): ');
            if choice==1
                disp(I);
                stem(1:L,D); title('Diff plot');
                xlabel('Nucleotide position'); ylabel('');
                AXIS([1 L 0 2])
                set(gca,'XTick',1:1:L);     
                set(gca,'YTick',0:1:3); set(gca,'YTickLabel','');
            end
        end
    end
end
