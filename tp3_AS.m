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

