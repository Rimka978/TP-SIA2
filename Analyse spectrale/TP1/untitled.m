clear all ;
close all;
clc;

%%
Fe = 1000;
Te = 1/Fe;
f0 = 100;
N= 100;

t= (0:N-1)/Fe; 
f=(-(N-1)/2: (N-1)/2)*Fe/N ;
% Tf d'une sinusoide 

x1 = cos(2*pi*f0*t);
figure;
plot(t,x1); title("Signal Sinusoidal ");

TF_x1 = fft(x1);
figure;
plot(f,fftshift(abs(TF_x1))); title("Spectre Signal Sinusoidal ")

%normalisation 
x1_norm = TF_x1 * 2/N;
figure;
plot(f,fftshift(abs(x1_norm)));
title("Spectre Signal Sinusoidal normalisé")

%  angle
%avec la phase 

x1_phase = angle(TF_x1);


figure; plot(f,fftshift(angle(x1_phase)));
hold all
plot(f,abs(fftshift((x1_phase)))); title("Signal Sinusoidal f=100Hz ");
legend("phase signal", "spectre");

%% 3 pour f0 = 95

f0 = 95;
N= 100;

t= (0:N-1)/Fe; 
f=(-(N-1)/2: (N-1)/2)*Fe/N ;
% Tf d'une sinusoide 

x2 = cos(2*pi*f0*t);
figure;
plot(t,x2); title("Signal Sinusoidal ");

TF_x2 = fft(x2);
figure;
plot(f,fftshift(abs(TF_x2))); title("Spectre Signal Sinusoidal ")

%normalisation 
x1_norm = TF_x2 * 2/N;
figure;
plot(f,fftshift(abs(x1_norm)));
title("Spectre Signal Sinusoidal normalisé")

%avec la phase 

x2_phase = angle(TF_x2);


figure; plot(f,fftshift(angle(x2_phase)));
hold all
plot(f,abs(fftshift((x2_phase)))); title("Signal Sinusoidal f=95Hz");
legend("phase signal", "spectre");

%calcul de l'erreur
erreur_relatif_x1 = abs((1 - max(TF_x1)*2/N)/1);
erreur_relatif_x2 = abs((1 - max(TF_x2)*2/N)/1);

%% 5

N1=1000;
t1= (0:N1-1)/Fe; 
f1=(-(N1-1)/2: (N1-1)/2)*Fe/N1 ;

%Signal f=100Hz
x1_1=cos(2*pi*100*t1);
TF_x1_1=fft(x1_1);

figure;
plot(t1,x1_1);title("Signal Sinusoidal");

figure;
plot(f1,fftshift(abs(TF_x1_1))); title("Spectre Signal Sinusoidal ")

X1_1_norm1=fft(x1_1)*(2/N1);
figure;
plot(f1,fftshift(abs(X1_1_norm1))); title("Spectre Signal Sinusoidal normalisé")

%Signal f=99.5Hz
x2_2=cos(2*pi*99.5*t1);
TF_x2_2=fft(x2_2);

figure;
plot(t1,x2_2); title("Signal Sinusoidal");

figure;
plot(f1,fftshift(abs(TF_x2_2))); title("Spectre Signal Sinusoidal")

x2_2_norm=fft(x2_2)*(2/N1);
figure;
plot(f1,fftshift(abs(x2_2_norm))); title("Spectre Signal Sinusoidal normalisé")

%calcul de l'erreur
erreur_relatif2_x1_1 = abs((1 - max(TF_x1_1)*2/N1)/1);
erreur_relatif2_x2_2 = abs((1 - max(TF_x2_2)*2/N1)/1);


%% 9

N=100;
t= (0:N-1)/Fe; 
f=(-(N-1)/2: (N-1)/2)*Fe/N ;

%Signal avec zero padding f=100Hz
x1_100= zeros(1,N*10);
x1_100(1:N)= cos(2*pi*100*t);

X1_100=fft(x1_100);
X9_100_norm=fft(x1_100)*2/N ;

%pour N=1000
t3= (0:N*10-1)/Fe; 
f3=(-((N*10)-1)/2: ((N*10)-1)/2)*Fe/(N*10) ;


figure;
plot(t3,x1_100);title("Signal Sinusoidal avec zero padding f=100Hz");
figure;
plot(f3,fftshift(abs(X1_100))); title("Spectre Signal Sinusoidal avec zero padding f=100Hz")
figure;
plot(f3,fftshift(abs(X9_100_norm))); title("Spectre Signal Sinusoidal normalisé avec zero padding f=100Hz")


%Signal avec zero padding f=95Hz
x2_95= zeros(1,N*10);
x2_95(1:N)=cos(2*pi*95*t);

X2_95=fft(x2_95);
X9_95_norm=fft(x2_95)*2/N ;

figure;
plot(t3,x2_95); title("Signal Sinusoidal avec zero padding f=95Hz");
figure;
plot(f3,fftshift(abs(X2_95))); title("Spectre Signal Sinusoidal avec zero padding f=95Hz")
figure;
plot(f3,fftshift(abs(X9_95_norm))); title("Spectre Signal Sinusoidal normalisé avec zero padding f=95Hz")

%calcul de l'erreur
erreur_relatif3_100= abs((1 - max(X9_100_norm))/1)
erreur_relatif3_95= abs((1 - max(X9_95_norm))/1)

%% 10

N=100; 
t= (0:(N/2)-1)/Fe; 
f=(-((N/2)-1)/2: ((N/2)-1)/2)*Fe/(N/2) ;

t10= (0:N-1)/Fe; 
f10=(-((N)-1)/2: ((N)-1)/2)*Fe/N ;

%Signal avec zero padding N/2 f=100Hz
x10_100= zeros(1,N);
x10_100(1:N/2)= cos(2*pi*100*t);

X10_100=fft(x10_100);
X10_100_norm=X10_100*4/N

figure;
plot(t10,x10_100);title("Signal Sinusoidal avec zero padding f=100Hz");
figure;
plot(f10,fftshift(abs(X10_100))); title("Spectre Signal Sinusoidal avec zero padding N/2 f=100Hz")
figure;
plot(f10,fftshift(abs(X10_100_norm))); title("Spectre Signal Sinusoidal normalisé avec zero padding N/2 f=100Hz")


%Signal avec zero padding N/2 f=95Hz
x10_95= zeros(1,N);
x10_95(1:N/2)=cos(2*pi*95*t);

X10_95=fft(x10_95);
X10_95_norm=fft(x10_95)*4/N ;

figure;
plot(t10,x10_95); title("Signal Sinusoidal avec zero padding N/2 f=95Hz");
figure;
plot(f10,fftshift(abs(X10_95))); title("Spectre Signal Sinusoidal avec zero padding N/2 f=95Hz")
figure;
plot(f10,fftshift(abs(X10_95_norm))); title("Spectre Signal Sinusoidal normalisé avec zero padding N/2 f=95Hz")

%calcul de l'erreur
erreur_relatif4_100= abs((1 - max(X10_100_norm))/1)
erreur_relatif4_95= abs((1 - max(X10_95_norm))/1)
