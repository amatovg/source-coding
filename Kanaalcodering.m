%% 3.2 Channel Coding

%% 3.2.1 Cyclic Hamming Code
    
% See Verslag.docx for the checkpolynomial, the error-correcting and 
% -detecting performance of this code.
%
% The generator- and checkmatrix can be created by using the functions
% Generator_Matrix() and Check_Matrix() respectively, while the systematic 
% forms can be created by using Systematic_Generator_Matrix() and 
% Systematic_Check_Matrix(), which all can be found in Channel_Coding.m.
% The matrices themselves can be found in Verslag.docx

% % decomment code:
% generator = [1 0 0 1 1];
% check = [1 1 1 1 0 1 0 1 1 0 0 1];
% G = Channel_Coding.Generator_Matrix(generator)
% Gsys = Channel_Coding.Systematic_Generator_Matrix(generator)
% H = Channel_Coding.Check_Matrix(check)
% Hsys = Channel_Coding.Systematic_Check_Matrix(generator)

%% 3.2.2 Syndrome Table

% The syndrome table can be created by using Syndrome_Table.m. 

% % decomment code:
% [syndromes,cosets]=Syndrome_Table(Hsys)

% The probability of a decoding-error can be found in Verslag.docx.

%% 3.2.3 Hamming Encoding and Decoding + Hamming Code Simulations

% The Ham_encode() and Ham_decode() functions can be found in 
% Channel_Coding.m. 
% 
% For the simulation, see decodeErrorHamming.m and Verslag.docx.

% % decomment code:
% decodeErrorHamming



%% 3.2.4 Product Encoding

% For the function Prod_encode(), see Channel_Coding.m.


%% 3.2.5 Product Decoding

% For the function Prod_decode(), see Channel_Coding.m.

%% 3.2.6 Product Code Simulations

% For the simulation, see decodeErrorProduct.m and Verslag.docx.

% % decomment code:
%decodeErrorProduct


%% 3.2.7 Sending the Images (uncompressed)



% See Verslag.docx

[im,Y] = imread('bobmarley.bmp');
bitseq = Image_to_Bits(im);
[height,width] = size(im);

% Case 1: Uncoded transmission

fig = figure;
endbits = Channel_Data(bitseq,0.01);
endim = Bits_to_Image(endbits,height,width);
imshow(endim,Y);
imwrite(endim,Y,'Kanaalcodering_i.bmp')
unco = endbits;

% Case 2: Hamming-coded transmission

fig = figure;
encoded_bits = Channel_Coding.Ham_encode(bitseq);
channel_bits = Channel_Data(encoded_bits,0.01);
decoded_bits = Channel_Coding.Ham_decode(channel_bits);
endbits = decoded_bits(1:(height*width));
endim = Bits_to_Image(endbits,height,width);
imshow(endim,Y);
imwrite(endim,Y,'Kanaalcodering_ii.bmp')
ham = endbits;

% Case 3: Product-coded transmission

fig = figure;
encoded_bits = Channel_Coding.Prod_encode(bitseq);
channel_bits = Channel_Data(encoded_bits,0.01);
decoded_bits = Channel_Coding.Prod_decode(channel_bits);
endbits = decoded_bits(1:(height*width));
endim = Bits_to_Image(endbits,height,width);
imshow(endim,Y);
imwrite(endim,Y,'Kanaalcodering_iii.bmp')

prod = endbits;

% aantal fouten
N_errors_uncoded = sum(bitxor(bitseq,unco))
N_errors_hamming = sum(bitxor(bitseq,ham))
N_errors_product = sum(bitxor(bitseq,prod))




