% CHAR2NUM : "Convert from characters to numbers"
% This function maps the characters of the input sequence to the
% corresponding genetic alphabet representation in GF(5)
% A = 1, G = 2, C = 3, T/U = 4
% Usage: y = char2num(x, Case, Type);
% 
% (Case = 1 : capital letters ; Case = 0 : small letters)
% (Type = 1 : RNA sequence ; Type = 0 : DNA sequence)

function y = char2num(x, Case, Type)

y = zeros(size(x));

if Type==0
    if Case==0
        y(find(x=='a')) = 1;
        y(find(x=='c')) = 3;
        y(find(x=='g')) = 2;
        y(find(x=='t')) = 4;
    end
    if Case==1
        y(find(x=='A')) = 1;
        y(find(x=='C')) = 3;
        y(find(x=='G')) = 2;
        y(find(x=='T')) = 4;
    end
end

if Type==1
    if Case==0
        y(find(x=='a')) = 1;
        y(find(x=='c')) = 3;
        y(find(x=='g')) = 2;
        y(find(x=='u')) = 4;
    end
    if Case==1
        y(find(x=='A')) = 1;
        y(find(x=='C')) = 3;
        y(find(x=='G')) = 2;
        y(find(x=='U')) = 4;
    end
end


I = find(y==0);
if length(I)~=0
    fprintf('\nCheck if all the characters at the following positions in the sequence belong to the specified input alphabet: \n'); 
    disp(I);
    error('Mapping of bases done incorrectly!');
end
