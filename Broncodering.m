%%% 3.1 SOURCE CODING


%% 3.1.1 relative frequencies of the macrosymbols in the image

% import image
[image,Y] = imread('bobmarley.bmp');

% macrosymbols:
%   ms1 = 0000 = 00     ms2 = 0001 = 00     ms3 = 0010 = 00     ms4 = 0011 = 00
%                00                  01                  10                  11   
%   
%   ms5 = 0100 = 01     ms6 = 0101 = 01     ms7 = 0110 = 01     ms8 = 0111 = 01
%                00                  01                  10                  11   
%   
%   ms9 = 1000 = 10     ms10 = 1001 = 10    ms11 = 1010 = 10    ms12 = 1011 = 10
%                00                   01                  10                  11   
%   
%   ms13 = 1100 = 11    ms14 = 1101 = 11    ms15 = 1110 = 11    ms16 = 1111 = 11
%                 00                  01                  10                  11    
%
% alphabet
alphabet = {'ms1', 'ms2', 'ms3', 'ms4', 'ms5', 'ms6', 'ms7', 'ms8', 'ms9', 'ms10', 'ms11', 'ms12', 'ms13', 'ms14', 'ms15', 'ms16'};
N_ms = numel(alphabet);

% dimensions of the image
[M,N]=size(image);

% counter holds the number of times a macrosymbol is in the image
counter=zeros(1,16);

% counting the macrosymbols
for m=1:2:M
    for n=1:2:N
        x=image(m,n)*2^3+image(m,n+1)*2^2+image(m+1,n)*2+image(m+1,n+1);
        counter(x+1)=counter(x+1)+1;
    end
end

% relative frequencies of the macrosymbols
rel_freq=counter/sum(counter);
% round to 4 digits after comma
rel_freq_round = round(rel_freq*10000)/10000;
relative_frequencies = cell(N_ms,2);
relative_frequencies(:,1)=alphabet';
for i = 1:length(relative_frequencies)
   relative_frequencies{i,2} =  num2str(rel_freq_round(i));
end

relative_frequencies

%% 3.1.2 Huffmancode 

codebook = Source_Coding.create_codebook(alphabet, rel_freq);

% codewords in bits as a string
codebook_Huffman = cell(N_ms,2);
codebook_Huffman(:,1) = alphabet';
codebook_Huffman(:,2) = codebook';

codebook_Huffman        % a matrix with in collum 1 the macrosymbols and in collum 2 the codewords
codewords=codebook_Huffman(:,end);



%% 3.1.3 Entropy

% see Entropy.m



%% 3.1.4 E(n)

% see Entropy.m



%% 3.1.5 Lower limit of E(n)

% see Verslag.docx



%% 3.1.6 Canonical Huffmancode image

% see Verslag.docx



%% 3.1.7 Canonical Huffmancode example

% see Verlsag.docx



%% 3.1.8 Canonical Huffmancode image

% lengths of the codewords
lengths = cellfun(@length,codebook);

canonical_codebook = Source_Coding.create_canonical_codebook(alphabet, lengths);

% codewords in bits as a string
canonical_Huffman = canonical_codebook';
for i=1:length(canonical_Huffman)
    canonical_Huffman{i}= num2str(canonical_Huffman{i});
end
canonical_Huffman = [alphabet' canonical_Huffman]










