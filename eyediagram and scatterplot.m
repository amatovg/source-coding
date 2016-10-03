
b=random_bitstring(100);
T=10^-6;
Ns=8;
frequency=3*10^6;
alpha=1;
Lf = 5*Ns;
sigma=200;
x=[1:100];
x2=[1:50];
a=PHY.mapper(b,'BPSK');
a2=PHY.mapper(b,'QPSK');
K=100;
t=[-Lf/Ns:1/Ns:Lf/Ns+(K-1)]*T;
t2=[-Lf/Ns:1/Ns:Lf/Ns+(K/2-1)]*T;


for k = 1:K
                var_t = [-Lf/Ns+(k-1):1/Ns:Lf/Ns+(k-1)]*T; %10Ns+1 monsterwaarden maar verschoven interval
                u(k,1+(k-1)*Ns:2*Lf+1+(k-1)*Ns) = sqrt(1)*(a(k)*PHY.pulse(var_t-(k-1)*T,T,alpha));
            end
            % laatste rij van U vormt u(t)
            u(K+1,1:2*Lf+(K-1)*Ns+1) = sum(u(1:K,1:2*Lf+(K-1)*Ns+1));

fig = figure;
plot(x,a,'b.')
xlabel('symbolen');
ylabel('a');
title('Plot van symbolen BPSK');
print(fig, '-djpeg', 'Mapping a-BPSK.jpg');

fig = figure;
plot(x2,a2,'b.')
xlabel('symbolen');
ylabel('a2');
title('Plot van symbolen QPSK');
print(fig, '-djpeg', 'Mapping a-QPSK.jpg');

s = PHY.modulate(a,T,Ns,frequency,alpha,Lf);
s2= PHY.modulate(a2,T,Ns,frequency,alpha,Lf);

fig = figure;
plot(t,s,'b.:');
xlabel('t');
ylabel('s(t)');
title('Zendersignaal BPSK');
print(fig, '-djpeg', 'draaggolfsignaal s-BPSK.jpg');

fig = figure;
plot(t2,s2,'b.:');
xlabel('t');
ylabel('s(t)');
title('Zendersignaal QPSK');
print(fig, '-djpeg', 'draaggolfsignaal S-QSPK.jpg');

r = PHY.channel_AWGN(s, Ns, sigma);
r2 = PHY.channel_AWGN(s2, Ns, sigma);

fig = figure;
plot(t,r)
xlabel('t');
ylabel('ontvangerpuls r');
title('Ontvanger R BPSK');
print(fig, '-djpeg', 'ontvanger r-BPSK.jpg');

fig = figure;
plot(t2,r2)
xlabel('t');
ylabel('ontvangerpuls r');
title('Ontvanger R QPSK');
print(fig, '-djpeg', 'ontvanger r-QPSK.jpg');

%  mp = PHY.pulse([-Lf/Ns:1/Ns:Lf/Ns]*T,T,1);
%  fig = figure;
%  plot([-Lf/Ns:1/Ns:Lf/Ns]*T,mp,'b.')
%  xlabel('t');
%  ylabel('pulse ');
%  title('square-root raised cosine pulse');
%  print(fig, '-djpeg', 'square-root raised cosine pulse.jpg');

y = PHY.demodulate(r,T,Ns,frequency,alpha,Lf);
y2 = PHY.demodulate(r2,T,Ns,frequency,alpha,Lf);

fig = figure;
plot(t,y,'b.')
xlabel('t');
ylabel('ontvanger y');
title('Ontvanger Symbolen');
print(fig, '-djpeg', 'ontvanger y-BPSK.jpg');

fig = figure;
plot(t2,y2,'b.')
xlabel('t');
ylabel('ontvanger y');
title('Ontvanger Symbolen');
print(fig, '-djpeg', 'ontvanger y-QPSK.jpg');

z = PHY.demodulate(s,T,Ns,frequency,alpha,Lf);
z =  z(Lf-Ns+1:Lf +(K-1)*Ns +1);
z2 = PHY.demodulate(s2,T,Ns,frequency,alpha,Lf);
z2 = z2(Lf-Ns+1:Lf +(K/2-1)*Ns +1);


eyediagram(z,2*Ns,2*T,0) %eyediagram BPSK, signaal zonder ruis
eyediagram(z2,2*Ns,2*T,0) %eyediagram QPSK, signaal zonder ruis

w=PHY.downSample(y,Ns,Lf);
w2=PHY.downSample(y2,Ns,Lf);

scatterplot(w) %scatterplot BPSK, signaal met ruis
scatterplot(w2) %scatterplot QPSK, signaal met ruis