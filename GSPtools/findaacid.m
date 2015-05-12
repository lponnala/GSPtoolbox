% FINDAACID
% This function returns the 3-letter abbreviation for a one-letter amino acid

function y = findaacid(x, Case)

if Case==0
    switch x
        case 'f'
            y = 'Phe';
        case 'l'
            y = 'Leu';
        case 'i'
            y = 'Iso';
        case 'm'
            y = 'Met';
        case 'v'
            y = 'Val';
        case 's'
            y = 'Ser';
        case 'p'
            y = 'Pro';
        case 't'
            y = 'Thr';
        case 'a'
            y = 'Ala';
        case 'y'
            y = 'Tyr';
        case 'h'
            y = 'His';
        case 'q'
            y = 'Glu';
        case 'n'
            y = 'Asp';
        case 'k'
            y = 'Lys';
        case 'd'
            y = 'AspA';
        case 'e'
            y = 'GluA';
        case 'c'
            y = 'Cys';
        case 'w'
            y = 'Trp';
        case 'r'
            y = 'Arg';
        case 'g'
            y = 'Gly';
        otherwise
            error('amino acid not found!');
    end
end

if Case==1
    switch x
        case 'F'
            y = 'Phe';
        case 'L'
            y = 'Leu';
        case 'I'
            y = 'Iso';
        case 'M'
            y = 'Met';
        case 'V'
            y = 'Val';
        case 'S'
            y = 'Ser';
        case 'P'
            y = 'Pro';
        case 'T'
            y = 'Thr';
        case 'A'
            y = 'Ala';
        case 'Y'
            y = 'Tyr';
        case 'H'
            y = 'His';
        case 'Q'
            y = 'Glu';
        case 'N'
            y = 'Asp';
        case 'K'
            y = 'Lys';
        case 'D'
            y = 'AspA';
        case 'E'
            y = 'GluA';
        case 'C'
            y = 'Cys';
        case 'W'
            y = 'Trp';
        case 'R'
            y = 'Arg';
        case 'G'
            y = 'Gly';
        otherwise
            error('amino acid not found!');
    end
end
