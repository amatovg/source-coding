% Bepaal de BER voor beide constellaties (BPSK en QPSK) aan de hand van
% simulaties voor verschillende waarden van Eb/N0 (zie appendix A1.5 en
% A1.7) en maak een figuur. Beschouw enkel het interval Eb/N0 waarvoor
% 10^-1>BER>10^-4
    Es = 1; % zelf gekozen voor dit project om Es=1 te gebruiken
    M_BPSK = 2; % BPSK : 2 constellatiepunten
    M_QPSK = 4; % QPSK : 4 constellatiepunten
    m_BPSK = 1; % m = log2(M)
    m_QPSK = 2; % m = log2(M)
    T = 10^-6;
    frequency = 3*10^6;
    Ns = 8;
    Lf = 5*Ns;
    alpha = 1;
    k = 11; % lengte informatiewoorden : 11 bits
    n = 15; % lengte codewoorden bij kanaalcodering
%Grensen bepalen : 
    % Eb/No in dB
    E_N_dB = [-3:7];
    size_EN = 11;
    % Eb/No omgezet (niet in dB)
    E_N = zeros(1,size_EN);
    for i = 1:11
        E_N(i) = 10^(E_N_dB(i)/10);
    end
% sigma² = (No*Ns)/(2*T*Es) = (No*Ns)/(2*Eb*m*T) kanaal zonder codering
% sigma² = (No*Ns)/(2*T*Es) = (No*Ns*n)/(2*Eb*m*T*k) kanaal met codering
    N_E = ones(size(E_N))./E_N; % Eb/No -> No/Eb
    sigma_BPSK = (N_E*Ns/2/m_BPSK/T).^(1/2);  
    sigma_BPSK_coded = (N_E*Ns*n/2/T/k/m_BPSK).^(1/2);
    sigma_QPSK = (N_E*Ns/2/m_QPSK/T).^(1/2);
    sigma_QPSK_coded = (N_E*Ns*n/2/m_QPSK/T/k).^(1/2);
    sigma_DQPSK = sigma_QPSK;
    sigma_DQPSK_coded = sigma_QPSK_coded;
% Variables
    BER_BPSK = zeros(1,size_EN);
    BER_QPSK = zeros(1,size_EN);
    BER_DQPSK = zeros(1,size_EN);
    BER_DQPSK_phasenoise = zeros(1,size_EN);
    BER_c_BPSK = zeros(1,size_EN);
    BER_c_QPSK = zeros(1,size_EN);
    BER_c_DQPSK = zeros(1,size_EN);
    BER_c_DQPSK_phasenoise = zeros(1,size_EN);

for i = 1:size_EN
    K = 100; % aantal te versturen symbolen
    % ONGECODEERD KANAAL
        % BER = Ne/Nbits met Ne = foute bits in foute symbolen, Nbits = m*K =
        % totaal van verstuurde bits 
            Nbits_BPSK = K*m_BPSK;
            Nbits_QPSK = K*m_QPSK;
            Nbits_DQPSK = K*m_QPSK;
        % Kanaal : ongecodeerd
            b_BPSK = random_bitstring(K*m_BPSK);
            b_QPSK = random_bitstring(K*m_QPSK);
            b_DQPSK = random_bitstring(K*m_QPSK);
            
            a_BPSK = PHY.mapper(b_BPSK,'BPSK');
            a_QPSK = PHY.mapper(b_QPSK,'QPSK');
            a_DQPSK = PHY.mapper(b_DQPSK,'D-4PSK');
            
            s_BPSK = PHY.modulate(a_BPSK,T,Ns,frequency,alpha,Lf);
            s_QPSK = PHY.modulate(a_QPSK,T,Ns,frequency,alpha,Lf);
            s_DQPSK = PHY.modulate(a_DQPSK,T,Ns,frequency,alpha,Lf);
            
            r_BPSK = PHY.channel_AWGN(s_BPSK,Ns,sigma_BPSK(i));
            r_QPSK = PHY.channel_AWGN(s_QPSK,Ns,sigma_QPSK(i));
            r_DQPSK = PHY.channel_AWGN(s_DQPSK,Ns,sigma_DQPSK(i));
            r_DQPSK_phasenoise = PHY.channel_phaseNoise(s_DQPSK, Ns, sigma_DQPSK(i), 0.05);
            
            y_BPSK = PHY.demodulate(r_BPSK,T,Ns,frequency,alpha,Lf);
            y_QPSK = PHY.demodulate(r_QPSK,T,Ns,frequency,alpha,Lf);
            y_DQPSK = PHY.demodulate(r_DQPSK,T,Ns,frequency,alpha,Lf);
            y_DQPSK_phasenoise = PHY.demodulate(r_DQPSK_phasenoise,T,Ns,frequency,alpha,Lf);
            
            y_down_BPSK = PHY.downSample(y_BPSK,Ns,Lf);
            y_down_QPSK = PHY.downSample(y_QPSK,Ns,Lf);
            y_down_DQPSK = PHY.downSample(y_DQPSK,Ns,Lf);
            y_down_DQPSK_phasenoise = PHY.downSample(y_DQPSK_phasenoise,Ns,Lf);
            
            a_d_BPSK = PHY.harddecisions(y_down_BPSK,'BPSK');
            a_d_QPSK = PHY.harddecisions(y_down_QPSK,'QPSK');
            a_d_DQPSK = PHY.harddecisions(y_down_DQPSK,'D-4PSK');
            a_d_DQPSK_phasenoise = PHY.harddecisions(y_down_DQPSK_phasenoise,'D-4PSK');
            
            b_d_BPSK = PHY.demapper(a_d_BPSK,'BPSK');
            b_d_QPSK = PHY.demapper(a_d_QPSK,'QPSK');
            b_d_DQPSK = PHY.demapper(a_d_DQPSK,'D-4PSK');
            b_d_DQPSK_phasenoise = PHY.demapper(a_d_DQPSK_phasenoise,'D-4PSK');

            Ne_BPSK = sum(bitxor(b_BPSK,b_d_BPSK));
            Ne_QPSK = sum(bitxor(b_QPSK,b_d_QPSK));
            Ne_DQPSK = sum(bitxor(b_DQPSK,b_d_DQPSK));
            Ne_DQPSK_phasenoise = sum(bitxor(b_DQPSK,b_d_DQPSK_phasenoise));

    % GECODEERD KANAAL : voor BPSK en QPSK
            k=11;
            Nbits_c_BPSK = K*k;
            Nbits_c_QPSK = K*k;
            Nbits_c_DQPSK = K*k;
            % Kanaal
            b_c_BPSK = random_bitstring(K*k); % Te versturen bitsequentie
            b_c_QPSK = random_bitstring(K*k);
            b_c_DQPSK = random_bitstring(K*k);

            c_BPSK = Channel_Coding.Ham_encode(b_c_BPSK); % Geencodeerde bitsequentie
            c_QPSK = Channel_Coding.Ham_encode(b_c_QPSK);
            c_DQPSK = Channel_Coding.Ham_encode(b_c_DQPSK);

            a_c_BPSK = PHY.mapper(c_BPSK,'BPSK'); % Datasymbolen
            a_c_QPSK = PHY.mapper(c_QPSK,'QPSK');
            a_c_DQPSK = PHY.mapper(c_DQPSK,'D-4PSK');

            s_c_BPSK = PHY.modulate(a_c_BPSK,T,Ns,frequency,alpha,Lf); % Modulatie
            s_c_QPSK = PHY.modulate(a_c_QPSK,T,Ns,frequency,alpha,Lf);
            s_c_DQPSK = PHY.modulate(a_c_DQPSK,T,Ns,frequency,alpha,Lf);

            r_c_BPSK = PHY.channel_AWGN(s_c_BPSK,Ns,sigma_BPSK_coded(i)); % Kanaal
            r_c_QPSK = PHY.channel_AWGN(s_c_QPSK,Ns,sigma_QPSK_coded(i));
            r_c_DQPSK = PHY.channel_AWGN(s_c_DQPSK,Ns,sigma_DQPSK_coded(i));
            r_c_DQPSK_phasenoise = PHY.channel_phaseNoise(s_c_DQPSK, Ns, sigma_DQPSK_coded(i), 0.05);

            y_c_BPSK = PHY.demodulate(r_c_BPSK,T,Ns,frequency,alpha,Lf); % demodulatie
            y_c_QPSK = PHY.demodulate(r_c_QPSK,T,Ns,frequency,alpha,Lf);
            y_c_DQPSK = PHY.demodulate(r_c_DQPSK,T,Ns,frequency,alpha,Lf);
            y_c_DQPSK_phasenoise = PHY.demodulate(r_c_DQPSK_phasenoise,T,Ns,frequency,alpha,Lf);

            y_c_down_BPSK = PHY.downSample(y_c_BPSK,Ns,Lf); % Downsampling
            y_c_down_QPSK = PHY.downSample(y_c_QPSK,Ns,Lf);
            y_c_down_DQPSK = PHY.downSample(y_c_DQPSK,Ns,Lf);
            y_c_down_DQPSK_phasenoise = PHY.downSample(y_c_DQPSK_phasenoise,Ns,Lf);

            a_c_d_BPSK = PHY.harddecisions(y_c_down_BPSK,'BPSK'); % decisie
            a_c_d_QPSK = PHY.harddecisions(y_c_down_QPSK,'QPSK');
            a_c_d_DQPSK = PHY.harddecisions(y_c_down_DQPSK,'D-4PSK');
            a_c_d_DQPSK_phasenoise = PHY.harddecisions(y_c_down_DQPSK_phasenoise,'D-4PSK');

            c_d_BPSK = PHY.demapper(a_c_d_BPSK,'BPSK');% demapping
            c_d_QPSK = PHY.demapper(a_c_d_QPSK,'QPSK');
            c_d_DQPSK = PHY.demapper(a_c_d_DQPSK,'D-4PSK');
            c_d_DQPSK_phasenoise = PHY.demapper(a_c_d_DQPSK_phasenoise,'D-4PSK');

            b_c_d_BPSK = Channel_Coding.Ham_decode(c_d_BPSK); % decodering
            b_c_d_QPSK = Channel_Coding.Ham_decode(c_d_QPSK);
            b_c_d_DQPSK = Channel_Coding.Ham_decode(c_d_DQPSK);
            b_c_d_DQPSK_phasenoise = Channel_Coding.Ham_decode(c_d_DQPSK_phasenoise);

            Ne_c_BPSK = sum(bitxor(b_c_d_BPSK,b_c_BPSK));
            Ne_c_QPSK = sum(bitxor(b_c_d_QPSK,b_c_QPSK));
            Ne_c_DQPSK = sum(bitxor(b_c_d_DQPSK,b_c_DQPSK));
            Ne_c_DQPSK_phasenoise = sum(bitxor(b_c_d_DQPSK_phasenoise,b_c_DQPSK));
    

    i %Afprinten waar we zitten in de lus (voor de ongeduldige wachters)
    BER_BPSK(1,i) = Ne_BPSK/Nbits_BPSK;
    BER_QPSK(1,i) = Ne_QPSK/Nbits_QPSK;
    BER_DQPSK(1,i) = Ne_DQPSK/Nbits_DQPSK;
    BER_DQPSK(1,i) = Ne_DQPSK_phasenoise/Nbits_DQPSK;
    BER_c_BPSK(1,i) = Ne_c_BPSK/Nbits_c_BPSK;
    BER_c_QPSK(1,i) = Ne_c_QPSK/Nbits_c_QPSK;
    BER_c_DQPSK(1,i) = Ne_c_DQPSK/Nbits_c_DQPSK;
    BER_c_DQPSK_phasenoise(1,i) = Ne_c_DQPSK_phasenoise/Nbits_c_DQPSK;
end

QBER = qfunc(sqrt(E_N*2));

fig = figure;
semilogy(E_N_dB,BER_BPSK,'ko-',E_N_dB,BER_QPSK,'ks-',E_N_dB,BER_DQPSK,'kx-',E_N_dB,BER_DQPSK_phasenoise,'k*-',E_N_dB,QBER,'r');
xlabel('Eb/No (dB)');
ylabel('BER');
title('Schatting BER zonder kanaalcodering');
legend('BPSK : AWGN','QPSK : AWGN','D-4PSK : AWGN','D-4PSK : phasenoise','Theoretische BER','Location','Southwest');
print(fig, '-djpeg', 'BER_zonder_kanaalcodering.jpg');

fig = figure;
semilogy(E_N_dB,BER_c_BPSK,'ko-',E_N_dB,BER_c_QPSK,'ks-',E_N_dB,BER_c_DQPSK,'kx-',E_N_dB,BER_c_DQPSK_phasenoise,'k*-',E_N_dB,QBER,'r');
xlabel('Eb/No (dB)');
ylabel('BER');
title('Schatting BER met kanaalcodering');
legend('BPSK : AWGN','QPSK : AWGN','D-4PSK : AWGN','D-4PSK : phasenoise','Theoretische BER','Location','Southwest');
print(fig, '-djpeg', 'BER_met_kanaalcodering.jpg');