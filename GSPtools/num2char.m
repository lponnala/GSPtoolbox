% NUM2CHAR
% This function maps a numeric (GF(5)) sequence to the
% corresponding genetic alphabet 
% 1 = A, 2 = G, 3 = C, 4 = T/U
% Usage: y = num2char(x, Case, Type);
% 
% (Case = 1 : capital letters ; Case = 0 : small letters)
% (Type = 1 : RNA sequence ; Type = 0 : DNA sequence)

function y = num2char(x, Case, Type);

if prod(x)==0
    error('Input DNA sequence is invalid! It has zeros in it!');
end

y = cell(size(x));

if Type==0
    if Case==0
        y(find(x==1)) = {'a'};
        y(find(x==3)) = {'c'};
        y(find(x==2)) = {'g'};
        y(find(x==4)) = {'t'};
    end
    if Case==1
        y(find(x==1)) = {'A'};
        y(find(x==3)) = {'C'};
        y(find(x==2)) = {'G'};
        y(find(x==4)) = {'T'};
    end
end
    
if Type==1
    if Case==0
        y(find(x==1)) = {'a'};
        y(find(x==3)) = {'c'};
        y(find(x==2)) = {'g'};
        y(find(x==4)) = {'u'};
    end
    if Case==1
        y(find(x==1)) = {'A'};
        y(find(x==3)) = {'C'};
        y(find(x==2)) = {'G'};
        y(find(x==4)) = {'U'};
    end
end

y = char(y)';
