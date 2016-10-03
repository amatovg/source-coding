classdef PHY
    
    methods(Static=true)
        
        %% Functies voor Mapping/Demapping
        function a = mapper(bitstring, constellation)
            % Functie die een bitstring mapped naar complexe symbolen.
            % input:
            % bitstring: vector van ongecodeerde bits met lengte 1xN
            % constellation: 'BPSK', 'QPSK'
            % output: vector met (complexe) symbolen
            
            switch(constellation)
                case 'BPSK'
                    % BPSK mapping here
                    a = (bitstring-1/2)*2;
                case 'QPSK'
                    % QPSK mapping here
                    Mapping = [1 0+1j 0-1j -1]; % Volgorde bepaald door graymapping
                    N = size(bitstring,2);
                    if(mod(N,2)==1) % we komen in de problemen bij oneven bitstrings
                        error('Uneven number of bits');
                    end
                    a = zeros(1,N/2);
                    for i = 1 : N/2
                        a(i) = Mapping(bi2de(bitstring(2*i-1:2*i),'left-msb')+1);
                    end
                case 'D-4PSK'
                    % 3.3.10 : niet-coherente detectie : differentieel
                    % geencodeerde 4-PSK
                    Mapping = [0 pi/2 -pi/2 pi];
                    N = size(bitstring,2);
                    if(mod(N,2)==1) % we komen in de problemen bij oneven bitstrings
                        error('Uneven number of bits');
                    end
                    a = zeros(1,N/2+1); % het eerst verstuurde symbool dient om de faseverschillen te kunnen berekenen.
                    a(1) = 1;
                    arg = 0;
                    for i = 1 : N/2
                        arg_diff = Mapping(bi2de(bitstring(2*i-1:2*i),'left-msb')+1);
                        arg=arg+arg_diff;
                        a(i+1)=exp(1j*arg);
                    end
                otherwise
                    error('Constellation not recognized');
            end
            
        end
        function bitstring = demapper(a, constellation)
            % Functie die complexe symbolen omzet naar bits
            % input:
            % a: vector met (complexe symbolen)
            % constellation: ofwel 'BPSK', 'QAM', 'QPSK', '4PAM', etc
            
            switch(constellation)
                case 'BPSK'
                    bitstring = (a+1)/2;
                case 'QPSK'
                    N_size = size(a);
                    N = N_size(2);
                    bitstring = zeros(1,N*2);
                    for i = 1 : N
                        switch(a(1,i))
                            case 1
                                bitstring(1,2*i-1) = 0;
                                bitstring(1,2*i) = 0;
                            case -1
                                bitstring(1,2*i-1) = 1;
                                bitstring(1,2*i) = 1;
                            case 0+1j
                                bitstring(1,2*i-1) = 0;
                                bitstring(1,2*i) = 1;
                            case 0-1j
                                bitstring(1,2*i-1) = 1;
                                bitstring(1,2*i) = 0;
                            otherwise
                                error ('Datasymbol not recognized');
                        end
                    end
                case 'D-4PSK'
                    % demap differentieel geencodeerde 4-PSK hier
                    N=size(a,2)-1;
                    arg=angle(a(1));
                    bitstring = zeros(1,2*N);
                    C=1e-5;
                    for i = 1 : N
                        arg_diff = angle(a(i+1))-arg;
                        arg=angle(a(i+1));
                        
                        if(abs(arg_diff)<= C || abs(arg_diff-2*pi)<=C)
                            bitstring(2*i-1:2*i)=[0 0];
                        elseif(abs(arg_diff-pi/2)<=C || abs(arg_diff+3*pi/2)<=C)
                            bitstring(2*i-1:2*i)=[0 1];    
                        elseif(abs(arg_diff+pi/2)<=C || abs(arg_diff-3*pi/2)<=C)
                            bitstring(2*i-1:2*i)=[1 0];    
                        elseif(abs(arg_diff-pi)<=C || abs(arg_diff+pi)<=C)
                            bitstring(2*i-1:2*i)=[1 1];    
                        else
                            error ('Datasymbol not recognized');
                        end
                    end                    
                    
                otherwise
                    error('Constellation not recognized');
            end
        end
        function a_estim = harddecisions(x, constellation)
            % Functie die harde desicie toepast op x
            % input:
            % x: vector met ruizige (complexe) symbolen
            % constellation ofwel 'BPSK', '4QAM', 'QPSK', '4PAM', etc
            
            switch(constellation)
                case 'BPSK'
                    
                    a_estim(x<0) = -1;
                    a_estim(x>0) = 1;
                case 'QPSK'
                    
                    a = real(x);
                    b = imag(x);
                    a_estim(a>0  & abs(a)>=abs(b)) = 1+0j;
                    a_estim(a<=0 & abs(a)>abs(b))  = -1+0j;
                    a_estim(b>=0 & abs(a)<abs(b))  = 0+1j;
                    a_estim(b<0  & abs(a)<=abs(b)) = 0-1j;
                    
                case 'D-4PSK'
                    
                    % harde decisie voor differentieel geencodeerde 4-PSK hier
                    
                    a = real(x);
                    b = imag(x);
                    a_estim(a>0  & abs(a)>=abs(b)) = 1+0j;
                    a_estim(a<=0 & abs(a)>abs(b))  = -1+0j;
                    a_estim(b>=0 & abs(a)<abs(b))  = 0+1j;
                    a_estim(b<0  & abs(a)<=abs(b)) = 0-1j;
                    
                otherwise
                    error('Constellation not recognized');
                    
            end
            
            
        end
        
  %% moduleren/demoduleren
        function s = modulate(a,T,Ns,frequency,alpha,Lf)
            % Zet een vector symbolen op puls en moduleert dit signaal op een golf
            % input:
            % a = vector of symbols
            % T = symbol duration = 1 microsec = 10^-6 s
            % Ns = samples per symbol = 8
            % frequency = carrier frequency = 3*10^6 Hz
            % alpha = roll-off factor = 1
            % Lf = pulse duration (in samples) = 5*Ns = 40
            % s = vector containing samples of modulated signal
            K_lengtes = size(a);
            K = K_lengtes(2);
            Es = 1;
            u = zeros(K+1,2*Lf+(K-1)*Ns+1);
            % eerste K rijen vormen telkens een term uit u(t)
            for k = 1:K
                var_t = [-Lf/Ns+(k-1):1/Ns:Lf/Ns+(k-1)]*T; %10Ns+1 monsterwaarden maar verschoven interval
                u(k,1+(k-1)*Ns:2*Lf+1+(k-1)*Ns) = sqrt(Es)*(a(k)*PHY.pulse(var_t-(k-1)*T,T,alpha));
            end
            % laatste rij van U vormt u(t)
            u(K+1,1:2*Lf+(K-1)*Ns+1) = sum(u(1:K,1:2*Lf+(K-1)*Ns+1));
            % s(t)
            t = [-Lf/Ns:1/Ns:Lf/Ns+(K-1)]*T;
            s = sqrt(2)*real(u(K+1,1:2*Lf+(K-1)*Ns+1).*exp(0+1j*2*pi*frequency*t));
            % Plot van symbolen op zenderpulsen
%             fig = figure;
%             plot(t,u,'b.:');
%             xlabel('t');
%             ylabel('u(t)');
%             title('complex basisbandsignaal');
%             print(fig, '-djpeg', 'complex basisbandsignaal.jpg');
            % Plot van zendersignaal
%             plot(t,s);
%             xlabel('t');
%             ylabel('s(t)');
%             title('Zendersignaal');
            % Vergelijking van de fouriergetransformeerden van de originele
            % puls en de afgeknotte puls is te vinden in
            % Square_root_raised_cosine_pulse.m
            
        end        
       function rdemod = demodulate(r,T,Ns,frequency,alpha,Lf)
            % Demoduleert het signaal r en voert het matched filter uit. 
            % input:
            % r: vector of received signal samples
            % T: symbol duration
            % Ns: samples per symbol
            % frequency: carrier frequency
            % alpha: roll-off factor
            % Lf: pulse duration
            %
            % output y: vector of demodulated samples
            
            % Es = 1, kanaal is genormaliseerd
            Es=1;
            lengte = size(r);
            rl = zeros(1,lengte(2));
            t=(-Lf/Ns:1/Ns:Lf/Ns)*T;
            for i = 1:lengte(2);
                rl(i) = r(i)/sqrt(Es);
                rl(i) = sqrt(2)*rl(i)*exp(-1j*2*pi*frequency*(i-1)*(T/Ns));
            end
            rdemod = (T/Ns)*conv(rl,PHY.pulse(t,T,alpha),'same');
        end           
        function y = pulse(t,T,alpha)
            % Functie: De square root raised cosine pulse
            % input:
            % t sample waarden
            % T tijdsinterval van 1 symbool: standaard gelijk aan 1 
            % microseconde, seconde? )
            % alpha: rolloff factor
            % vb van gebruik
            %
            % alpha = 0.5;
            % t = [-5:0.1:5];
            % s = PHY.pulse(t, 1, alpha);
            % plot(t, s)
            % xlabel('tijd t/T');
            % ylabel('y(t)');
            % title(['Square root raised cosine pulse met rollofffactor ' num2str(alpha)]);
                        
            een=(1-alpha)*sinc(t*(1-alpha)/T);
            twee=(alpha)*cos(pi*(t/T-0.25)).*sinc(alpha*t/T-0.25);
            drie=(alpha)*cos(pi*(t/T+0.25)).*sinc(alpha*t/T+0.25);
            y=1/(sqrt(T))*(een+twee+drie);

        end
        
        function y=downSample(x,Ns,Lf)
            % Functie die het gegeven signaal downsampled
            % input:
            % x: vector met samples
            % Ns: samples per symbool
            % Lf: pulse duration
            % output: rijvector met 1 sample per symbool
            
            y = x(Lf+1:Ns:size(x,2)-Lf);
        end   
        
        %% het kanaal
        function r = channel_AWGN(s, Ns, sigma)
            % Functie die een AWGN kanaal simuleert
            % input:
            % s: rijvector met samples van het signaal
            % Ns: samples per symbol
            % sigma: standaard deviatie van de ruis.
            
            %Es = 1, omdat het kanaal genormaliseerd is
            Es=1;
            lengte=size(s);
            r = zeros(1,lengte(2));
            % e^(-j*theta) zitten er nog niet ingerekend.
            for i = 1:(lengte(2))
               r(i) = s(i) + normrnd(0,sigma);
            end
        end
        
        function r = channel_phaseNoise(s, Ns, sigma, faseNoiseStd)
            % Functie die een kanaal simuleert waarbij de
            % draaggolfoscillator van de ontvanger onderhevig is aan
            % fasefluctuaties ('faseruis').
            % input:
            % s: rijvector met samples van het signaal
            % Ns: samples per symbol
            % sigma: standaard deviatie van de ruis.
            % faseNoiseStd: standaard deviatie van het faseincrement voor één symboolperiode
            
            N = length(s);            
            n_theta = zeros(1, N);
            for k=1:N
                n_theta(k) = n_theta(k-1) + faseNoiseStd/sqrt(Ns)*randn();    
            end
		
		% faseruis toegepast op het signaal
		s = s.*exp(1j*n_theta);
            
            % voeg nog steeds AWGN toe
            r = PHY.channel_AWGN(s, Ns, sigma);
        end      
        
        
    end
    
end