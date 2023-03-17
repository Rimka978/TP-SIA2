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




%%  II-2 fenetrage de la TFD

%II.3 3)
N=100;
t= (0:N-1)/Fe; 
n = (0:N-1);
f=(-(N-1)/2: (N-1)/2)*Fe/N ;

%fenetre rectangulaire
wr= ones(1,N);
Wr= fft(wr);
%fenetre Hanning
wh= 0.5*(1- cos(2*pi*n/(N-1)));
Wh= fft(wh);
%fenetre Blackman
wb=(27/64)- 0.5*cos((2*pi*n)/(N-1)) + (5/64)*cos((4*pi*n)/(N-1));
Wb= fft(wb);


figure;
subplot(3,1,1); plot(t,wr); title("fenetre rectangulaire");
subplot(3,1,2); plot(t,wh); title("fenetre Hanning");
subplot(3,1,3); plot(t,wb); title("fenetre Blackman");

figure;
subplot(3,1,1); plot(f,fftshift(abs(Wr))); title("spectre fenetre rectangulaire");
subplot(3,1,2); plot(f,fftshift(abs(Wh))); title("spectre fenetre Hanning");
subplot(3,1,3); plot(f,fftshift(abs(Wb))); title("spectre fenetre Blackman");

%avec zero padding

N=100;
t= (0:N-1)/Fe; 
f=(-(N-1)/2: (N-1)/2)*Fe/N ;
n = (0:N-1);

n10 = (0:N*10-1);
f10=(-(N*10-1)/2: (N*10-1)/2)*Fe/(N*10) ;

w=zeros(1,N*10);

%fenetre rectangulaire
wr=w;
wr(1:N)= ones(1,N);
Wr= fft(wr);
%fenetre Hanning
wh=w;
wh(1:N)= 0.5*(1- cos(2*pi*n/(N-1)));
Wh= fft(wh);
%fenetre Blackman
wb=w;
wb(1:N)=(27/64)- 0.5*cos((2*pi*n)/(N-1)) + (5/64)*cos((4*pi*n)/(N-1));
Wb= fft(wb);

figure;
subplot(3,1,1); plot(n10,wr); title("fenêtre rectangulaire avec zero padding");
subplot(3,1,2); plot(n10,wh); title("fenêtre Hanning avec zero padding");
subplot(3,1,3); plot(n10,wb); title("fenêtre Blackman avec zero padding");

figure;
subplot(3,1,1); plot(f10,fftshift(abs(Wr))); title("spectre fenêtre rectangulaire avec zero padding");
subplot(3,1,2); plot(f10,fftshift(abs(Wh))); title("spectre fenêtre Hanning avec zero padding");
subplot(3,1,3); plot(f10,fftshift(abs(Wb))); title("spectre fenêtre Blackman avec zero padding");


%normalisation des fentres

%Normalisation des fenetre
w=zeros(1,N*10);
%fenetre rectangulaire
wr=w;
wr(1:N)= ones(1,N);
Wr= fft(wr);
Wr_norm= Wr/max(Wr);

%fenetre Hanning
wh=w;
wh(1:N)= 0.5*(1- cos(2*pi*n/(N-1)));
Wh= fft(wh);
Wh_norm= Wh/max(Wh);

%fenetre Blackman
wb=w;
wb(1:N)=(27/64)- 0.5*cos((2*pi*n)/(N-1)) + (5/64)*cos((4*pi*n)/(N-1));
Wb= fft(wb);
Wb_norm= Wb/max(Wb);

figure;
subplot(3,1,1); plot(f10,fftshift(abs(Wr_norm))); title("spectre fenêtre rectangulaire normalisé avec zero padding");
subplot(3,1,2); plot(f10,fftshift(abs(Wh_norm))); title("spectre fenêtre Hanning normalisé avec zero padding");
subplot(3,1,3); plot(f10,fftshift(abs(Wb_norm))); title("spectre fenêtre Blackman normalisé avec zero padding");

%%
%II.3 5)

N=100;
t= (0:N-1)/Fe; 
f=(-(N-1)/2: (N-1)/2)*Fe/N ;
n = (0:N-1);

%fenetre rectangulaire
wr= ones(1,N);
Wr= fft(wr);
Wr_norm= Wr/max(Wr);

%fenetre Hanning
wh= 0.5*(1- cos(2*pi*n/(N-1)));
Wh= fft(wh);
Wh_norm= Wh/max(Wh);

%fenetre Blackman
wb=(27/64)- 0.5*cos((2*pi*n)/(N-1)) + (5/64)*cos((4*pi*n)/(N-1));
Wb= fft(wb);
Wb_norm= Wb/max(Wb);

%Signal 100Hz
x_100=cos(2*pi*100*t);
X_100=fft(x_100);

%Multiplication par la fenetre
vr_100=x_100.*wr;
vh_100=x_100.*wh;
vb_100=x_100.*wb;

Vr_100=fft(vr_100);
Vh_100=fft(vh_100);
Vb_100=fft(vb_100);


figure;
subplot(3,1,1); plot(t,real(vr_100)); title("Signal f=100Hz fenêtré rectangulaire");
subplot(3,1,2); plot(t,real(vh_100)); title("Signal f=100Hz fenêtré Hanning");
subplot(3,1,3); plot(t,real(vb_100)); title("Signal f=100Hz fenêtré Blackman");

figure;
subplot(3,1,1); plot(f,fftshift(abs(Vr_100))); title("spectre Signal f=100Hz fenêtré rectangulaire");
subplot(3,1,2); plot(f,fftshift(abs(Vh_100))); title(" spectre Signal f=100Hz fenêtré Hanning");
subplot(3,1,3); plot(f,fftshift(abs(Vb_100))); title("spectre Signal  f=100Hz fenêtré Blackman");

%Signal 95Hz
x_95=cos(2*pi*95*t);
X_95=fft(x_95);

%Multiplication par la fenetre
vr_95=x_95.*wr;
vh_95=x_95.*wh;
vb_95=x_95.*wb;

Vr_95=fft(vr_95);
Vh_95=fft(vh_95);
Vb_95=fft(vb_95);


figure;
subplot(3,1,1); plot(t,real(vr_95)); title("Signal f=95Hz fenêtré rectangulaire");
subplot(3,1,2); plot(t,real(vh_95)); title("Signal f=95Hz fenêtré Hanning");
subplot(3,1,3); plot(t,real(vb_95)); title("Signal f=95Hz fenêtré Blackman");

figure;
subplot(3,1,1); plot(f,fftshift(abs(Vr_95))); title("spectre Signal f=95Hz fenêtré rectangulaire");
subplot(3,1,2); plot(f,fftshift(abs(Vh_95))); title(" spectre Signal f=95Hz fenêtré Hanning");
subplot(3,1,3); plot(f,fftshift(abs(Vb_95))); title("spectre Signal  f=95Hz fenêtré Blackman");


%% 6 
Vr_100_norm=Vr_100*2/N;
Vh_100_norm=Vh_100*(2/max(Wh));
Vb_100_norm=Vb_100*(2/max(Wb));

figure;
subplot(3,1,1); plot(f,fftshift(abs(Vr_100_norm))); title("spectre Signal fenêtré normalisé f=100Hz, rectangulaire");
subplot(3,1,2); plot(f,fftshift(abs(Vh_100_norm))); title(" spectre Signal fenêtré normalisé f=100Hz, Hanning");
subplot(3,1,3); plot(f,fftshift(abs(Vb_100_norm))); title("spectre Signal fenêtré normalisé f=100Hz, Blackman");

Vr_95_norm=Vr_95*2/N;
Vh_95_norm=Vh_95*(2/max(Wh));
Vb_95_norm=Vb_95*(2/max(Wb));

figure;
subplot(3,1,1); plot(f,fftshift(abs(Vr_95_norm))); title("spectre Signal fenêtré normalisé f=95Hz, rectangulaire");
subplot(3,1,2); plot(f,fftshift(abs(Vh_95_norm))); title(" spectre Signal fenêtré normalisé f=95Hz, Hanning");
subplot(3,1,3); plot(f,fftshift(abs(Vb_95_norm))); title("spectre Signal fenêtré normalisé f=95Hz, Blackman");

%TFD sans zero padding
Vr_100_1=Vr_100_norm;
Vh_100_1=Vh_100_norm;
Vb_100_1=Vb_100_norm;

Vr_95_1=Vr_95_norm;
Vh_95_1=Vh_95_norm;
Vb_95_1=Vb_95_norm;

%TFD avec zero padding
Vr_100_2=Vr_100_norm;
Vh_100_2=Vh_100_norm;
Vb_100_2=Vb_100_norm;

Vr_95_2=Vr_95_norm;
Vh_95_2=Vh_95_norm;
Vb_95_2=Vb_95_norm;


%%
%
%Superposition

figure;
subplot(3,1,1); plot(f,fftshift(abs(Vr_100_1))); hold on ; plot(f10,fftshift(abs(Vr_100_2))) ;
legend("sans zero padding","avec zero padding");
title("spectre Signal fenêtré normalisé f=100Hz, rectangulaire");
subplot(3,1,2); plot(f,fftshift(abs(Vh_100_1)));  hold on ; plot(f10,fftshift(abs(Vh_100_2))) ;
legend("sans zero padding","avec zero padding");
title(" spectre Signal fenêtré normalisé f=100Hz, Hanning");
subplot(3,1,3); plot(f,fftshift(abs(Vb_100_1)));  hold on ; plot(f10,fftshift(abs(Vb_100_2))); 
legend("sans zero padding","avec zero padding");
title("spectre Signal fenêtré normalisé f=100Hz, Blackman");

figure;ss
subplot(3,1,1); plot(f,fftshift(abs(Vr_95_1))); hold on ; plot(f10,fftshift(abs(Vr_95_2))) ;
legend("sans zero padding","avec zero padding");
title("spectre Signal fenêtré normalisé f=95Hz, rectangulaire");
subplot(3,1,2); plot(f,fftshift(abs(Vh_95_1)));  hold on ; plot(f10,fftshift(abs(Vh_95_2))) ;
legend("sans zero padding","avec zero padding");
title(" spectre Signal fenêtré normalisé f=95Hz, Hanning");
subplot(3,1,3); plot(f,fftshift(abs(Vb_95_1)));  hold on ; plot(f10,fftshift(abs(Vb_95_2))); 
legend("sans zero padding","avec zero padding");
title("spectre Signal fenêtré normalisé f=95Hz, Blackman");

%% III
N=100;
Fe=1000;
f=(-(N-1)/2: (N-1)/2)*Fe/N ;
n = (0:N-1);

A1=1;
f1=95;
phi1=0;


%A2=0.1;
f2=140;
phi2=0;

y= A1*sin(2*pi*n*f1/Fe + phi1) + A2*sin(2*pi*n*f2/Fe + phi2);
Y=fft(y); 

w=zeros(1,N);
%fenetre rectangulaire
wr=w;
wr(1:N)= ones(1,N);
Wr= fft(wr);
Wr_norm= Wr/max(Wr);

%fenetre Hanning
wh=w;
wh(1:N)= 0.5*(1- cos(2*pi*n/(N-1)));
Wh= fft(wh);
Wh_norm= Wh/max(Wh);

%fenetre Blackman
wb=w;
wb(1:N)=(27/64)- 0.5*cos((2*pi*n)/(N-1)) + (5/64)*cos((4*pi*n)/(N-1));
Wb= fft(wb);
Wb_norm= Wb/max(Wb);


yfr= y.*wr ;
yfh= y.*wh ;
yfb= y.*wb ;

Yfr= fft(yfr)*(2/max(Wr));
Yfh= fft(yfh)*(2/max(Wh));
Yfb= fft(yfb)*(2/max(Wb));

figure;
subplot(3,1,1); plot(n,real(yfr)); title(sprintf("Signal A1=1 A2= %d fenetré rectangulaire", A2));
subplot(3,1,2); plot(n,real(yfh)); title(sprintf("Signal A1=1 A2= %d fenetré Hanning",A2));
subplot(3,1,3); plot(n,real(yfb)); title(sprintf("Signal A1=1 A2= %d fenetré Blackman", A2));

figure;
subplot(3,1,1); plot(f,fftshift(abs(Yfr))); title(sprintf("spectre Signal A1=1 A2= %d fenêtré rectangulaire", A2));
subplot(3,1,2); plot(f,fftshift(abs(Yfh))); title(sprintf(" spectre Signal A1=1 A2= %d fenêtré Hanning",A2));
subplot(3,1,3); plot(f,fftshift(abs(Yfb))); title(sprintf("spectre Signal  A1=1 A2= %d fenêtré Blackman", A2));


%%
%Q3:


N=100;
Fe=1000;
f=(-(N-1)/2: (N-1)/2)*Fe/N ;
n = (0:N-1);

A1=1;
f1=95;
phi1=0;


A2=0.3;
f2=140;
phi2=0;
G=0.05;
bruit=  sqrt(G)*randn(1,N);
yb= A1*sin(2*pi*n*f1/Fe + phi1) + A2*sin(2*pi*n*f2/Fe + phi2) + bruit;
Yb=fft(yb); 

w=zeros(1,N);
%fenetre rectangulaire
wr=w;
wr(1:N)= ones(1,N);
Wr= fft(wr);
Wr_norm= Wr/max(Wr);

%fenetre Hanning
wh=w;
wh(1:N)= 0.5*(1- cos(2*pi*n/(N-1)));
Wh= fft(wh);
Wh_norm= Wh/max(Wh);

%fenetre Blackman
wb=w;
wb(1:N)=(27/64)- 0.5*cos((2*pi*n)/(N-1)) + (5/64)*cos((4*pi*n)/(N-1));
Wb= fft(wb);
Wb_norm= Wb/max(Wb);


yfr= yb.*wr ;
yfh= yb.*wh ;
yfb= yb.*wb ;

Yfr= fft(yfr)*(2/max(Wr));
Yfh= fft(yfh)*(2/max(Wr));
Yfb= fft(yfb)*(2/max(Wr));

figure;
subplot(3,1,1); plot(n,real(yfr)); title(sprintf("Signal bruité de sigma= %d  fenetré rectangulaire", G));
subplot(3,1,2); plot(n,real(yfh)); title(sprintf("Signal bruité de sigma= %d  fenetré Hanning",G));
subplot(3,1,3); plot(n,real(yfb)); title(sprintf("Signal bruité de sigma= %d  fenetré Blackman", G));

figure;
subplot(3,1,1); plot(f,fftshift(abs(Yfr))); title(sprintf("spectre Signal bruité de sigma= %d  fenêtré rectangulaire", G));
subplot(3,1,2); plot(f,fftshift(abs(Yfh))); title(sprintf(" spectre Signal bruité de sigma= %d  Hz fenêtré Hanning",G));
subplot(3,1,3); plot(f,fftshift(abs(Yfb))); title(sprintf("spectre Signal  bruité de sigma= %d  fenêtré Blackman", G));

%variance trouvé = 0.05


%%
N=100;
Fe=1000;
f=(-(N-1)/2: (N-1)/2)*Fe/N ;
n = (0:N-1);

A1=1;
%f1=100;
f1=95;
phi1=0;


A2=A1;
f2=180;
phi2=0;

y = A1*sin(2*pi*n*f1/Fe + phi1) %+ A2*sin(2*pi*n*f2/Fe + phi2);
y(33:66)=0;
Y=fft(y); 



w=zeros(1,N);
%fenetre rectangulaire
wr=w;
wr(1:N)= ones(1,N);
Wr= fft(wr);
Wr_norm= Wr/max(Wr);

%fenetre Hanning
wh=w;
wh(1:N)= 0.5*(1- cos(2*pi*n/(N-1)));
Wh= fft(wh);
Wh_norm= Wh/max(Wh);

%fenetre Blackman
wb=w;
wb(1:N)=(27/64)- 0.5*cos((2*pi*n)/(N-1)) + (5/64)*cos((4*pi*n)/(N-1));
Wb= fft(wb);
Wb_norm= Wb/max(Wb);


yfr= y.*wr ;
yfh= y.*wh ;
yfb= y.*wb ;

Yfr= fft(yfr);
Yfh= fft(yfh);
Yfb= fft(yfb);

figure;
subplot(3,1,1); plot(n,real(yfr)); title(sprintf("Signal f0=%d Hz fenetré rectangulaire", f1));
subplot(3,1,2); plot(n,real(yfh)); title(sprintf("Signal f0=%d Hz fenetré Hanning",f1));
subplot(3,1,3); plot(n,real(yfb)); title(sprintf("Signal f0=%d Hz fenetré Blackman", f1));

figure;
subplot(3,1,1); plot(f,fftshift(abs(Yfr))); title(sprintf("spectre Signal f0=%d Hz fenêtré rectangulaire", f1));
subplot(3,1,2); plot(f,fftshift(abs(Yfh))); title(sprintf(" spectre Signal f0=%d Hz fenêtré Hanning",f1));
subplot(3,1,3); plot(f,fftshift(abs(Yfb))); title(sprintf("spectre Signal  f0=%d Hz fenêtré Blackman", f1));

%Q4: enlever la seconde sinusoide et chercher le spectre
%Q6: ajouter la deuxième sinusoïde avec une même amplitude 


%%% un chirp

%%
%Chirp 

N=100;
Fe=1000;
t= (0:N-1)/Fe; 
f=(-(N-1)/2: (N-1)/2)*Fe/N ;
n = (0:N-1);

a= 100;
b= 1010;


x= sin(2*pi*( a + b*t).*t);

X=fft(x); 


figure;
plot(t,x); title("Signal chirp");
figure;
plot(f,fftshift(abs(X))); title("Spectre Signal chirp"); 

%fenetre rectangulaire
wr= ones(1,N);
Wr= fft(wr);
Wr_norm= Wr/max(Wr);

%fenetre Hanning
wh= 0.5*(1- cos(2*pi*n/(N-1)));
Wh= fft(wh);
Wh_norm= Wh/max(Wh);

%fenetre Blackman
wb=(27/64)- 0.5*cos((2*pi*n)/(N-1)) + (5/64)*cos((4*pi*n)/(N-1));
Wb= fft(wb);
Wb_norm= Wb/max(Wb);

%Multiplication par la fenetre
vr=x.*wr;
vh=x.*wh;
vb=x.*wb;

Vr=fft(vr);
Vh=fft(vh);
Vb=fft(vb);


figure;
subplot(3,1,1); plot(t,real(vr)); title("Signal chirp fenêtré rectangulaire");
subplot(3,1,2); plot(t,real(vh)); title("Signal chirp fenêtré Hanning");
subplot(3,1,3); plot(t,real(vb)); title("Signal chirp fenêtré Blackman");

figure;
subplot(3,1,1); plot(f,fftshift(abs(Vr))); title("spectre Signal chirp fenêtré rectangulaire");
subplot(3,1,2); plot(f,fftshift(abs(Vh))); title(" spectre Signal chirp fenêtré Hanning");
subplot(3,1,3); plot(f,fftshift(abs(Vb))); title("spectre Signal  chirp fenêtré Blackman");

% QUESTION 2 
% APRES AVOIR TRACE LES DIFFERENTES TRANSFORMEES POUR LE SIGNAL ON VOIT
% QU'AVEC LA FENETRE RECTANGULAIRE ON NE VOIT RIEN TANDISQU 'AVEC LES DEUX
% AUTRES FENETREES ON ARRIVE A VOIR QUELQUE CHOSE AVEC CES DEUX FREQUENCES 



