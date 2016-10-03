BER=1e-3;
EbN0=qfuncinv(BER)^2/2;
T=10^-6;
Ns=8;
frequency = 3*10^6;
Lf=5*Ns;
alpha=1;

k=11;
n=15;
m_B = 1;
m_Q = 2;
sigma_uncoded_B=sqrt(Ns/EbN0/2/m_B/T);
sigma_uncoded_Q = sqrt(Ns/EbN0/2/m_Q/T);
sigma_coded_B = sqrt(Ns*n/EbN0/2/m_B/T/k);
sigma_coded_Q = sqrt(Ns*n/EbN0/2/m_Q/T/k);

vlag=0;

[image,Y] = imread('bobmarley.bmp');
bitstring = Image_to_Bits(image);
[height,width] = size(image);

% comprimeringoverhead

alphabet = {'A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P'};
[M,N]=size(image);
counter=zeros(1,16);
symbols = [];
for m=1:2:M
    for n=1:2:N
        x=image(m,n)*2^3+image(m,n+1)*2^2+image(m+1,n)*2+image(m+1,n+1);
        counter(x+1)=counter(x+1)+1;
        symbols = strcat(symbols,alphabet{x+1});
    end
end
rel_freq=counter/sum(counter);
codewords = Source_Coding.create_codebook(alphabet, rel_freq);


%% QPSK

constellation='QPSK';



%%% zonder codering of comprimering %%%



a=PHY.mapper(bitstring, constellation);
s = PHY.modulate(a,T,Ns,frequency,alpha,Lf);
r = PHY.channel_AWGN(s, Ns, sigma_uncoded_Q);
rdemod = PHY.demodulate(r,T,Ns,frequency,alpha,Lf);
y=PHY.downSample(rdemod,Ns,Lf);
a_estim = PHY.harddecisions(y, constellation);
bitstring2 = PHY.demapper(a_estim, constellation);

endim = Bits_to_Image(bitstring2,height,width);
imshow(endim,Y);
imwrite(endim,Y,'Complete_QPSK_i.bmp')

diff1=sum(bitxor(bitstring(1:size(bitstring2,2)),bitstring2))




%% zonder codering met comprimering %%%



bits_c = Source_Coding.Huffman_encode(symbols, alphabet, codewords);

if(mod(size(bits_c,2),2)==1)
    bits_c=[bits_c 0];
    vlag=1;
end

a=PHY.mapper(bits_c, constellation); 
s = PHY.modulate(a,T,Ns,frequency,alpha,Lf);
r = PHY.channel_AWGN(s, Ns, sigma_uncoded_Q);
rdemod = PHY.demodulate(r,T,Ns,frequency,alpha,Lf);
y=PHY.downSample(rdemod,Ns,Lf);
a_estim = PHY.harddecisions(y, constellation);
bits_c2 = PHY.demapper(a_estim, constellation);

if(vlag)
   bits_c2=bits_c2(1:size(bits_c2,2)-1);
   vlag = 0;
end

post_symbols = Source_Coding.Huffman_decode(bits_c2, alphabet, codewords);
bitstring2 = symbolsToBits(post_symbols);

endim = Bits_to_Image(bitstring2,height,width);
imshow(endim,Y);
imwrite(endim,Y,'Complete_QPSK_ii.bmp')

diff2=sum(bitxor(bitstring(1:size(bitstring2,2)),bitstring2))




%%% met codering met comprimering %%%



bits_c = Source_Coding.Huffman_encode(symbols, alphabet, codewords);

bits_c_c = Channel_Coding.Prod_encode(bits_c);

if(mod(size(bits_c_c,2),2)==1)
    bits_c_c=[bits_c_c 0];
    vlag=1;
end

a=PHY.mapper(bits_c_c, constellation); 
s = PHY.modulate(a,T,Ns,frequency,alpha,Lf);
r = PHY.channel_AWGN(s, Ns, sigma_uncoded_Q);
rdemod = PHY.demodulate(r,T,Ns,frequency,alpha,Lf);
y=PHY.downSample(rdemod,Ns,Lf);
a_estim = PHY.harddecisions(y, constellation);
bits_c_c2 = PHY.demapper(a_estim, constellation);

if(vlag)
   bits_c_c2=bits_c_c2(1:size(bits_c_c2,2)-1);
   vlag = 0;
end

bits_c2 = Channel_Coding.Prod_decode(bits_c_c2);
bits_c2 = bits_c2(1:size(bits_c,2));

post_symbols = Source_Coding.Huffman_decode(bits_c2, alphabet, codewords);
bitstring2 = symbolsToBits(post_symbols);

endim = Bits_to_Image(bitstring2,height,width);
imshow(endim,Y);
imwrite(endim,Y,'Complete_QPSK_iii.bmp')

diff3=sum(bitxor(bitstring(1:size(bitstring2,2)),bitstring2))
