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


%%

%4.1)Compréhension sur des signaux simulés
%1


% Construction du signal
  Fe = 8000;     % Fréquence d'échantillonnage
  dt =1/Fe; % Pas d'échantillonnage
  N  = 1024;     % Nombre d'échantillons du signal
  t  = (0:N-1)/Fe;     % Vecteur de temps
  f  = (-(N-1)/2:(N-1)/2)*Fe/N ;   % Vecteur de fréquence
  % Choix de la fenêtre
  L   =  10   % longueur de la fenêtre
  l=  0 : L-1 ;   % vecteur de la fenêtre

  %rect
fen = ones(1,L);

%Hanning:
fen= 0.5*(1- cos((2*pi*l)/(L-1)));

%Blackman:
fen=(27/64)- 0.5*cos((2*pi*l)/(L-1)) + (5/64)*cos((4*pi*l)/(L-1));
  
  d = 0:1:N-1;
  x = double(d==40);% Signal à analyser
% Parametres de la STFT
  NFFT = N   % taille des fft
  pas  = 1	   % pas d'échantillonnage en temps de la STFT


  %rect
fen1 = ones(1,L);

%Hanning:
fen2 = 0.5*(1- cos((2*pi*l)/(L-1)));

%Blackman:
fen3 =(27/64)- 0.5*cos((2*pi*l)/(L-1)) + (5/64)*cos((4*pi*l)/(L-1));



% Calcul de la Transformée de Fourier à court terme
  [X,tp,fp] = stft(x,NFFT,Fe,fen3,pas);


% Affichage
  figure(2); clf
% Transformée de Fourier à court terme
  subplot('position',[0.1 0.1 0.67 0.58])
  imagesc(tp,fp,abs(X)); axis xy
  title('Spectrogramme')
  xlabel('Temps (sec)')
  ylabel('fréquence (Hz)')
  Xlim = axis; Xlim=Xlim(1:2);

% Périodogramme
  subplot('position',[0.85 0.1 0.14 0.58])
  N2 = 8*N;
  freq = (0:N2-1)/N2*Fe; TFx = fft(x,N2);
  plot(abs(TFx(1:end/2)).^2/N,freq(1:end/2))
  xlabel('Puissance')
  title('Spectre')

% Représentation temporelle
  subplot('position',[0.1 0.75 0.67 0.2])
  plot(t,x)
  xlabel('Temps (sec)')
  ylabel('Amplitude')
  title('Représentation temporelle du signal')
  ax = axis; ax(1:2) = Xlim; axis(ax)

% Barre des couleurs
  hand = subplot('position',[0.01 0.1 0.02 0.58]);
  colorbar(hand)


%2
Fe=8000;
N=1024;
dt=1/Fe;
t=(0:N-1)/Fe;
f=(-(N-1)/2:(N-1)/2)*Fe/N ;
d = 0:1:N-1;
x1 = double(d==500);
x= cos(2*pi*2000*t) ;
t= 1/Fe*(1:length(x));

% Parametres de la STFT
NFFT = N; % taille des fft
pas = 1; % pas d'echantillonnage en temps de la STFT
% Choix de la fen^etre
L =60;% longueur de la fenetre
l=(0:L-1);

%rect
fen1 = ones(1,L);

%Hanning:
fen2= 0.5*(1- cos((2*pi*l)/(L-1)));

%Blackman:
fen3=(27/64)- 0.5*cos((2*pi*l)/(L-1)) + (5/64)*cos((4*pi*l)/(L-1));

% Calcul de la Transform ?ee de Fourier à court terme
[X,tp,fp]= stft(x,NFFT,Fe,fen1,pas);

% Affichage
figure(2)
% Transformee de Fourier `a court terme
subplot('position',[0.07 0.1 0.67 0.58])
imagesc(tp,fp,abs(X)); axis xy
title('Spectrogramme')
xlabel('Temps (sec)')
ylabel('fréquence (Hz)')
Xlim = axis; Xlim=Xlim(1:2);

% Periodogramme
subplot('position',[0.78 0.1 0.2 0.58])
N2 = 8*N;
freq = (0:N2-1)/N2*Fe; TFx = fft(x,N2);
plot(abs(TFx(1:end/2)).^2/N,freq(1:end/2))
xlabel('Puissance')
title('Spectre')

% Representation temporelle
subplot('position',[0.07 0.75 0.67 0.2])
plot(1/Fe*(1:length(x)),x)
xlabel('Temps (sec)')
ylabel('Amplitude')
title('Représentation temporelle du signal')
ax = axis; ax(1:2) = Xlim; axis(ax)


%3
Fe=8000;
N=1024;
dt=1/Fe;
t=(0:N-1)/Fe;
f=(-(N-1)/2:(N-1)/2)*Fe/N ;
x= [cos(2*pi*100*t(1:N/4)) ,  cos(2*pi*400*t((N/4)+1:end))];
t= 1/Fe*(1:length(x));

% Parametres de la STFT
NFFT = N; % taille des fft
pas = 1; % pas d2 ?echantillonnage en temps de la STFT
% Choix de la fen^etre
L = 100;% longueur de la fen^etre
l=(0:L-1);

%rect
fen = ones(1,L);

%Hanning:
fen= 0.5*(1- cos((2*pi*l)/(L-1)));

%Blackman:
fen=(27/64)- 0.5*cos((2*pi*l)/(L-1)) + (5/64)*cos((4*pi*l)/(L-1));

% Calcul de la Transform ?ee de Fourier à court terme
[X,tp,fp]= stft(x,NFFT,Fe,fen,pas);

% Affichage
figure(2)
% Transform ?ee de Fourier `a court terme
subplot('position',[0.07 0.1 0.67 0.58])
imagesc(tp,fp,abs(X)); axis xy
title('Spectrogramme')
xlabel('Temps (sec)')
ylabel('fréquence (Hz)')
Xlim = axis; Xlim=Xlim(1:2);

% P ?eriodogramme
subplot('position',[0.78 0.1 0.2 0.58])Fe=8000;
N=1024;
dt=1/Fe;
t=(0:N-1)/Fe;
f=(-(N-1)/2:(N-1)/2)*Fe/N ;
x= [cos(2*pi*100*t(1:N/4)) ,  cos(2*pi*400*t((N/4)+1:end))];
t= 1/Fe*(1:length(x));

% Parametres de la STFT
NFFT = N; % taille des fft
pas = 1; % pas d2 ?echantillonnage en temps de la STFT
% Choix de la fen^etre
L = 100;% longueur de la fen^etre
l=(0:L-1);

%rect
fen = ones(1,L);

%Hanning:
fen= 0.5*(1- cos((2*pi*l)/(L-1)));

%Blackman:
fen=(27/64)- 0.5*cos((2*pi*l)/(L-1)) + (5/64)*cos((4*pi*l)/(L-1));

% Calcul de la Transform ?ee de Fourier à court terme
[X,tp,fp]= stft(x,NFFT,Fe,fen,pas);

% Affichage
figure(2)
% Transform ?ee de Fourier `a court terme
subplot('position',[0.07 0.1 0.67 0.58])
imagesc(tp,fp,abs(X)); axis xy
title('Spectrogramme')
xlabel('Temps (sec)')
ylabel('fréquence (Hz)')
Xlim = axis; Xlim=Xlim(1:2);
Fe=8000;
N=1024;
dt=1/Fe;
t=(0:N-1)/Fe;
f=(-(N-1)/2:(N-1)/2)*Fe/N ;
x= [cos(2*pi*100*t(1:N/4)) ,  cos(2*pi*400*t((N/4)+1:end))];
t= 1/Fe*(1:length(x));

% Parametres de la STFT
NFFT = N; % taille des fft
pas = 1; % pas d2 ?echantillonnage en temps de la STFT
% Choix de la fen^etre
L = 100;% longueur de la fen^etre
l=(0:L-1);

%rect
fen = ones(1,L);

%Hanning:
fen= 0.5*(1- cos((2*pi*l)/(L-1)));

%Blackman:
fen=(27/64)- 0.5*cos((2*pi*l)/(L-1)) + (5/64)*cos((4*pi*l)/(L-1));

% Calcul de la Transform ?ee de Fourier à court terme
[X,tp,fp]= stft(x,NFFT,Fe,fen,pas);

% Affichage
figure(2)
% Transform ?ee de Fourier `a court terme
subplot('position',[0.07 0.1 0.67 0.58])
imagesc(tp,fp,abs(X)); axis xy
title('Spectrogramme')
xlabel('Temps (sec)')
ylabel('fréquence (Hz)')
Xlim = axis; Xlim=Xlim(1:2);

% P ?eriodogramme
subplot('position',[0.78 0.1 0.2 0.58])
N2 = 8*N;
freq = (0:N2-1)/N2*Fe; TFx = fft(x,N2);
plot(abs(TFx(1:end/2)).^2/N,freq(1:end/2))
xlabel('Puissance')
title('Spectre')

% Repr ?esentation temporelle
subplot('position',[0.07 0.75 0.67 0.2])
plot(1/Fe*(1:length(x)),x)
xlabel('Temps (sec)')
ylabel('Amplitude')
title('Représentation temporelle du signal')
ax = axis; ax(1:2) = Xlim; axis(ax)



% P ?eriodogramme
subplot('position',[0.78 0.1 0.2 0.58])
N2 = 8*N;
freq = (0:N2-1)/N2*Fe; TFx = fft(x,N2);
plot(abs(TFx(1:end/2)).^2/N,freq(1:end/2))
xlabel('Puissance')
title('Spectre')

% Repr ?esentation temporelle
subplot('position',[0.07 0.75 0.67 0.2])
plot(1/Fe*(1:length(x)),x)
xlabel('Temps (sec)')
ylabel('Amplitude')
title('Représentation temporelle du signal')
ax = axis; ax(1:2) = Xlim; axis(ax)



N2 = 8*N;
freq = (0:N2-1)/N2*Fe; TFx = fft(x,N2);
plot(abs(TFx(1:end/2)).^2/N,freq(1:end/2))
xlabel('Puissance')
title('Spectre')

% Repr ?esentation temporelle
subplot('position',[0.07 0.75 0.67 0.2])
plot(1/Fe*(1:length(x)),x)
xlabel('Temps (sec)')
ylabel('Amplitude')
title('Représentation temporelle du signal')
ax = axis; ax(1:2) = Xlim; axis(ax)


%4
Fe=8000;
N=1024;
dt=1/Fe;
t=(0:N-1)/Fe;
f=(-(N-1)/2:(N-1)/2)*Fe/N ;
a=100;
b=1010;
x= cos(2*pi*(a + b*t).*t);
t= 1/Fe*(1:length(x));


% Parametres de la STFT
NFFT = N; % taille des fft
pas = 1; % pas d2 ?echantillonnage en temps de la STFT
% Choix de la fen^etre
L = 100;% longueur de la fen^etre
l=(0:L-1);

%rect
fen1 = ones(1,L);

%Hanning:
fen2 = 0.5*(1- cos((2*pi*l)/(L-1)));

%Blackman:
fen3 =(27/64)- 0.5*cos((2*pi*l)/(L-1)) + (5/64)*cos((4*pi*l)/(L-1));

% Calcul de la Transform ?ee de Fourier à court terme
[X,tp,fp]= stft(x,NFFT,Fe,fen1,pas);

% Affichage
figure(2)
% Transform ?ee de Fourier `a court terme
subplot('position',[0.07 0.1 0.67 0.58])
imagesc(tp,fp,abs(X)); axis xy
title('Spectrogramme')
xlabel('Temps (sec)')
ylabel('fréquence (Hz)')
Xlim = axis; Xlim=Xlim(1:2);

% Periodogramme
subplot('position',[0.78 0.1 0.2 0.58])
N2 = 8*N;
freq = (0:N2-1)/N2*Fe; TFx = fft(x,N2);
plot(abs(TFx(1:end/2)).^2/N,freq(1:end/2))
xlabel('Puissance')
title('Spectre')

% Representation temporelle
subplot('position',[0.07 0.75 0.67 0.2])
plot(1/Fe*(1:length(x)),x)
xlabel('Temps (sec)')
ylabel('Amplitude')
title('Représentation temporelle du signal')
ax = axis; ax(1:2) = Xlim; axis(ax)


%% 4.2 Analyse de signaux réels



Fe=8000;
N=1024;
dt=1/Fe;
t=(0:N-1)/Fe;
f=(-(N-1)/2:(N-1)/2)*Fe/N ;

%[x fs]= audioread('haendel.wav');
%[x fs]= audioread('gong.wav');
[x fs]= audioread('laughter.wav');
t= 1/Fe*(1:length(x));
figure;
plot(t,x);
title('Signal laughter')
Xfft= fft(x);
figure;
N1=length(x);
f2=(-(N1-1)/2:(N1-1)/2)*Fe/N1 ;
plot(f,fftshift(abs(Xfft)));
title('FFT laughter')

% Parametres de la STFT
NFFT = N; % taille des fft
pas = 1; % pas d'echantillonnage en temps de la STFT
% Choix de la fenêtre
L = 110;% longueur de la fenêtre
l=(0:L-1);

%rect
fen1 = ones(1,L);

%Hanning:
fen2 = 0.5*(1- cos((2*pi*l)/(L-1)));

%Blackman:
fen3 =(27/64)- 0.5*cos((2*pi*l)/(L-1)) + (5/64)*cos((4*pi*l)/(L-1));

% Calcul de la Transform ?ee de Fourier à court terme
[X,tp,fp]= stft(x,NFFT,Fe,fen1,pas);

% Affichage
figure ;
% Transformee de Fourier à court terme
subplot('position',[0.07 0.1 0.67 0.58])
imagesc(tp,fp,abs(X)); axis xy
title('Spectrogramme')
xlabel('Temps (sec)')
ylabel('fréquence (Hz)')
Xlim = axis; Xlim=Xlim(1:2);

% Periodogramme
subplot('position',[0.78 0.1 0.2 0.58])
N2 = 8*N;
freq = (0:N2-1)/N2*Fe; TFx = fft(x,N2);
plot(abs(TFx(1:end/2)).^2/N,freq(1:end/2))
xlabel('Puissance')
title('Spectre')

% Representation temporelle
subplot('position',[0.07 0.75 0.67 0.2])
plot(1/Fe*(1:length(x)),x)
xlabel('Temps (sec)')
ylabel('Amplitude')
title('Représentation temporelle du signal')
ax = axis; ax(1:2) = Xlim; axis(ax)







