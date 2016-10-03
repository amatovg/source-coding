% Variabelen
Ns = 8;
Lf = 5*Ns;
T = 10^-6;
t = [-Lf/Ns:1/Ns:Lf/Ns]*T;


%___________________________________________________________________
% ALPHA = 1

pulse_1 = PHY.pulse([-Lf/Ns:1/Ns:Lf/Ns]*T,T,1);

%DTFT
% y fft van pulse  (fast fourier transform -> DTFT)
pulse_11 = PHY.pulse([-1:0.1:1]*T,T,1);
y_1=fft(pulse_11);

%____________________________________________________________________
% ALPHA = 0


pulse_2 = PHY.pulse([-Lf/Ns:1/Ns:Lf/Ns]*T,T,0);
fig = figure;
plot(t,pulse_1,'k.:',t,pulse_2,'r.:')
xlabel('t');
ylabel('pulse ');
legend('alpha = 1','alpha = 2','Location','Northeast');
title('square-root raised cosine pulse with rolloff');
print(fig, '-djpeg', 'square-root raised cosine pulse.jpg');

%DTFT
% y fft van pulse  (fast fourier transform -> DTFT)
pulse_21 = PHY.pulse([-1:0.1:1]*T,T,0);
y_2=fft(pulse_21);

%__________________________________________________________________
% plot of the dtft
fig = figure;
plot([-1:0.1:1],y_1/T,'k.:',[-1:0.1:1],y_2/T,'r.:');
ylabel('G(f)/T');
xlabel('f');
title('Fouriertransformatie met verschillende roll-off');
legend('alpha = 1','alpha = 0','Location','Southeast');
print(fig, '-djpeg', 'square-root raised cosine pulse FOURIER.jpg');