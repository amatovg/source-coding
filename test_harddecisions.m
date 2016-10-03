 bitstring=random_bitstring(3000);
 constellation='D-4PSK';
 a=PHY.mapper(bitstring, constellation);
 T=10^-6;
 Ns=8;
 frequency = 3*10^6;
 Lf=5*Ns;
 alpha=1;
 s = PHY.modulate(a,T,Ns,frequency,alpha,Lf);
 sigma=0;
 r = PHY.channel_AWGN(s, Ns, sigma);
 rdemod = PHY.demodulate(r,T,Ns,frequency,alpha,Lf);
 y=PHY.downSample(rdemod,Ns,Lf);
 a_estim = PHY.harddecisions(y, constellation);
 bitstring2 = PHY.demapper(a_estim, constellation);
 bitstring;
 bitstring2;
 diff=sum(bitxor(bitstring,bitstring2))