close all 
clear all

%3)

N=50;
%%Estimation de l'autocorrelation
for i=1:100
x1(:,i)=randn(1,N);
x2(:,i)=x1(:,i)+1;
R1_b(i,:)=xcorr(x1(:,i),'biased');
R2_b(i,:)=xcorr(x2(:,i),'biased');
R1_nb(i,:)=xcorr(x1(:,i),'unbiased');
R2_nb(i,:)=xcorr(x2(:,i),'unbiased');
end 

%moyenne
M1_b= mean(R1_b,1);
M2_b= mean(R2_b,1);
M1_nb= mean(R1_nb,1);
M2_nb= mean(R2_nb,1);

%Variance
V1_b= var(R1_b,1);
V2_b= var(R2_b,1);
V1_nb= var(R1_nb,1);
V2_nb= var(R2_nb,1);   

%4
%non biasé
% moyenne + ecart type
M1_ba=M1_b+std(R1_b);
M2_ba=M2_b+std(R2_b);
M1_nba=M1_nb+std(R1_nb);
M2_nba=M2_nb+std(R2_nb);
% moyenne - ecart type
M1_br=M1_b-std(R1_b,1);
M2_br=M2_b-std(R2_b,1);
M1_nbr=M1_nb-std(R1_nb,1);
M2_nbr=M2_nb-std(R2_nb,1);

figure;
subplot(2,1,1);
plot(M1_b); 
hold all
plot(M1_ba);
plot(M1_br);
legend('moyenne','écart ajouté','écart retranché');
title('estimateur 1 biaisé');

subplot(2,1,2);
plot(M2_b); 
hold all
plot(M2_ba);
plot(M2_br);
legend('moyenne','écart ajouté','écart retranché');
title('estimateur 2 biaisé');

figure;
subplot(2,1,1);
plot(M1_nb); 
hold all
plot(M1_nba);
plot(M1_nbr);
legend('moyenne','écart ajouté','écart retranché');
title('estimateur 1 non-biaisé');

subplot(2,1,2);
plot(M2_nb); 
hold all
plot(M2_nba);
plot(M2_nbr);
legend('moyenne','écart ajouté','écart retranché');
title('estimateur 2 non-biaisé');

%%
%III)1)
N=1024;
x1=randn(1,N);
%14)
N=1024;
n=0:N-1;

x1= randn(1,N);

f0= 0.1;
f1=0.15;
phi0= 2*pi*rand ;
phi1= 2*pi*rand ;
x3= cos(2*pi*f0*n + phi0) + cos(2*pi*f1*n + phi1);

figure;
subplot(2,1,1); plot(x1); title("Signal x1");
subplot(2,1,2); plot(x3); title("Signal x3");

L=[1 16 32 128];

    Vp3=zeros(4,1024);
    Mp3=zeros(4,1024);

    Vp1=zeros(4,1024);
    Mp1=zeros(4,1024);

for i=1:4
    M=N/L(i);
    for k=1:L(i)
        F3= fft(x3((k-1)*M+1:k*M),N);
        P3=(abs(F3).^2)/M;
        H3(i,:)=(1/L(i))*P3;

        F1= fft(x1((k-1)*M+1:k*M),N);
        P1=(abs(F1).^2)/M;
        H1(i,:)=(1/L(i))*P1;
    end
    Vp3(i,:)=var(H3(i,:),0,1);
    Mp3(i,:)=mean(H3(i,:),1);

    Vp1(i,:)=var(H1(i,:),0,1);
    Mp1(i,:)=mean(H1(i,:),1);
end
  
figure;
subplot(2,2,1); plot(1:1024,Mp1(1,:)); title(" périodogramme moyenné x1 L=1");
subplot(2,2,2) ; plot(1:1024,Mp1(2,:)); title("périodogramme moyenné x1 L=16");
subplot(2,2,3); plot(1:1024,Mp1(3,:)); title("périodogramme moyenné x1 L=32");
subplot(2,2,4); plot(1:1024,Mp1(4,:)); title("périodogramme moyenné x1 L=128");

figure;
subplot(2,2,1); plot(1:1024,Mp3(1,:)); title(" périodogramme moyenné x3 L=1");
subplot(2,2,2) ; plot(1:1024,Mp3(2,:)); title("périodogramme moyenné x3 L=16");
subplot(2,2,3); plot(1:1024,Mp3(3,:)); title("périodogramme moyenné x3 L=32");
subplot(2,2,4); plot(1:1024,Mp3(4,:)); title("périodogramme moyenné x3 L=128");




figure; plot(x1); title("signal X1");


R1_b=xcorr(x1,'biased');
R1_nb=xcorr(x1,'unbiased');

S1_b= fft(R1_b([N:end, 1:(N-1)]));
S1_nb= fft(R1_nb([N:end, 1:(N-1)]));




f=-(N-1):1:N-1 ;

figure;
subplot(2,1,1); plot(R1_b); title('Autocorrélation baisée');
subplot(2,1,2); plot(R1_nb); title('Autocorrélation non-baisée');

figure;
subplot(2,1,1); plot(abs(S1_b)); title('DSP baisée');
subplot(2,1,2); plot(abs(S1_nb)); title('DSP non-baisée');

%%
%III)2)
%11
N=1024;
P= 2*N -1 ;
Xfft=fft(x1, P);

Sx_p=(1/N) *(abs(Xfft).^2);
f= -(N-1):1:N-1;

figure; 
plot(f,fftshift(abs(Sx_p)),'*'); 
hold on
plot(f,fftshift(abs(S1_b))); 
legend('périodogramme Sp','corrélogramme Sb');
axis([-1010 1010  0  8]);


%% 13
N=128; 
P= 2*N -1 ;
f= -(N-1):1:N-1;

St= ones(1,255);

for i=1:100
x(:,i)=randn(1,N);
Xfft=fft(x(:,i), P);
Sx_pr(i,:)= fftshift(abs((1/N) *(abs(Xfft).^2)));
end 

MS_p= mean(Sx_pr,1);

figure; 
plot(f,MS_p); 
hold on
plot(f,St); 
legend('périodogramme Sp','DSP theo');
title("périodogramme 100 réalisations")

for i=1:1000
x(:,i)=randn(1,N);
Xfft=fft(x(:,i), P);
Sx_pr(i,:)= fftshift(abs((1/N) *(abs(Xfft).^2)));
end 

MS_p= mean(Sx_pr,1);

figure; 
plot(f,MS_p); 
hold on
plot(f,St); 
legend('périodogramme Sp','DSP theo');
title("périodogramme 1000 réalisations") 


%% periodogramme moyenné


%14)
N=1024;
n=0:N-1;

x1= randn(1,N);

f0= 0.1;
f1=0.15;
phi0= 2*pi*rand ;
phi1= 2*pi*rand ;
x3= cos(2*pi*f0*n + phi0) + cos(2*pi*f1*n + phi1);

figure;
subplot(2,1,1); plot(x1); title("Signal x1");
subplot(2,1,2); plot(x3); title("Signal x3");

L=[1 16 32 128];

    Vp3=zeros(4,1024);
    Mp3=zeros(4,1024);

    Vp1=zeros(4,1024);
    Mp1=zeros(4,1024);

for i=1:4
    M=N/L(i);
    for k=1:L(i)
        F3= fft(x3((k-1)*M+1:k*M),N);
        P3=(abs(F3).^2)/M;
        H3(i,:)=(1/L(i))*P3;

        F1= fft(x1((k-1)*M+1:k*M),N);
        P1=(abs(F1).^2)/M;
        H1(i,:)=(1/L(i))*P1;
    end
    Vp3(i,:)=var(H3(i,:),0,1);
    Mp3(i,:)=mean(H3(i,:),1);

    Vp1(i,:)=var(H1(i,:),0,1);
    Mp1(i,:)=mean(H1(i,:),1);
end
  
figure;
subplot(2,2,1); plot(1:1024,Mp1(1,:)); title(" périodogramme moyenné x1 L=1");
subplot(2,2,2) ; plot(1:1024,Mp1(2,:)); title("périodogramme moyenné x1 L=16");
subplot(2,2,3); plot(1:1024,Mp1(3,:)); title("périodogramme moyenné x1 L=32");
subplot(2,2,4); plot(1:1024,Mp1(4,:)); title("périodogramme moyenné x1 L=128");

figure;
subplot(2,2,1); plot(1:1024,Mp3(1,:)); title(" périodogramme moyenné x3 L=1");
subplot(2,2,2) ; plot(1:1024,Mp3(2,:)); title("périodogramme moyenné x3 L=16");
subplot(2,2,3); plot(1:1024,Mp3(3,:)); title("périodogramme moyenné x3 L=32");
subplot(2,2,4); plot(1:1024,Mp3(4,:)); title("périodogramme moyenné x3 L=128");


%%

%16

N = 1024;
n = (0:N-1);

f0 = 0.1;
phi0 = 2*pi*rand;
x4 = cos(2*pi*f0*n+phi0) + 4*randn(1,N);

figure;
plot(n,x4); title('Signal x4');

figure;
periodogram(x4); title('Périodogramme');
figure; 
pwelch(x4); title('Périodogramme de Welch');


%%

%21

N = 10000;
a = -0.7;
W = randn(2*N,1);
X=filter(1,[1 a],W);
X=X(N+1:end);

figure;
plot(X); title("Signal X");

Y=X + 4*randn(N,1);

figure;
plot(Y); title("Signal Y");


RY = xcorr(Y,'unbiased');
RX = xcorr(X,'unbiased');

RX0 = RX(N);
RX1 = RX(N+1);
RX2 = RX(N+2);

%estimateur non-robuste
a_enr= -RX1/RX0
sigma_enr = RX0 - (RX1^2)/RX0
%Estimateur robuste
a_er = -RX2/RX1
sigma_er = (RX1^2-RX2^2)/(RX2)




