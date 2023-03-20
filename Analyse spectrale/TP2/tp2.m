clear all ;
close all ;
clc ;

%% 2.1 Etude du filtre
%1
A = [1 -0.5 1];
B = [1 0 0 -0.729];

[H,W]=freqz(A,B);

figure;
plot(W,log10(abs(H)));
title("Réponse en fréquence du filtre H");
%Apres avoir tracé la reponnse en frrequence du filtre nous voyons qu'il se
%comporte comme un filtre coupe bande 

%2

[h,n]=impz(A,B);

figure;
plot(n,h);
title("Réponse impulsionnelle h");


%5
load('sig.mat');

P=50;
N=length(e);

Rse=xcorr(s,e);
Re=xcorr(e);
figure;
subplot(2,1,1); plot(1:2047,Re);
title("autocorrelation Re");
subplot(2,1,2); plot(1:2047,Rse);
title("intercorrelation Rse");

%calcul des matrice/vecteur
sigma= toeplitz(Re(N:N+P-1));
c= Rse(N:N+P-1);

%calcul de h estime
hest= inv(sigma)*c;

figure;
plot(1:50,hest);
hold on
plot(1:93,h);
legend("h estime","h");
title("Réponse impulsionnelle h estime");

%8
load('sig.mat');

P=50;
N=length(e);

Rse=xcorr(s,e);
Re=xcorr(e);
figure;
subplot(2,1,1); plot(1:2047,Re);
title("autocorrelation Re");
subplot(2,1,2); plot(1:2047,Rse);
title("intercorrelation Rse");

%calcul des matrice/vecteur
sigma= toeplitz(Re(N:N+P-1));
c= Rse(N:N+P-1);

%calcul de h estime
hest= inv(sigma)*c;

Se=fft(Re);
Sse=fft(Rse);

Hest= Sse./Se;
hest2=ifft(Hest);

figure;
plot(1:N,log10(abs(Hest(1:N))));
title("Réponse en fréquence du filtre h estime");

figure;
plot(1:50,hest2(1:P));
hold on
plot(1:93,h);
legend("h estime","h");
title("Réponse impulsionnelle h estime");


%% 3 Debruitage d'un signal avec référence de bruit 
%3.1) Application biomédicale
%2
load('estim.mat');
N=length(mes);
Rsb= xcorr(mes,ref);
Rb= xcorr(ref);

figure;
subplot(2,1,1); plot(1:7499,Rb);
title("autocorrelation Rb");
subplot(2,1,2); plot(1:7499,Rsb);
title("intercorrelation Rsb");

P=50;
sigma= toeplitz(Rb(N:N+P-1));
c= Rsb(N:N+P-1);

hest3= inv(sigma)*c;

figure;
plot(hest3);
title("Réponse impulsionnelle h estime 3");

%3
p_h= conv(ref,hest3);
e2 = mes - p_h(1:N);
figure;
plot(mes)

figure;
plot(1:N,e2);
hold on
plot(1:N,sig);
legend("signal estime","signal sig");
title("Signal estime");

%3.2) Application audio
%4
[z2 Fe]=audioread('z2.wav');
[b2 Fs]=audioread('b2.wav');

%5
P=300;
N=length(z2);
Rzb2= xcorr(z2,b2);
Rb2= xcorr(b2);

sigma= toeplitz(Rb2(N:N+P-1));
c= Rzb2(N:N+P-1);

hest4= inv(sigma)*c;

figure;
plot(hest4);
title("Réponse impulsionnelle h estime 4");

p_h = conv(b2,hest4);
e4= z2 - p_h(1:N);

figure;
plot(1:N,e4);
title("Signal estime");

%6
sound(e4 ,Fs);

