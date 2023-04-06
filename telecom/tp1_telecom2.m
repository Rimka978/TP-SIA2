clear all ;
clc ;


load('SignalModule.mat')

N = length(signalModule);
t = (0 : N-1)/fe;
fc = 1000 ;
f0 = 1000000;


f = (-(N-1)/2 : (N-1)/2)*fe/N;

figure;
plot(t(1:40000),signalModule(1:40000))


%transformée de fourier du signal 

Tf = fft(signalModule);
figure;
plot(f,fftshift(abs(Tf)))


% extraction des ak  et bk par demodulation 
%ak
a = 2*(signalModule.*cos(2*pi*f0*t));
filt = fir1(100,fc/(fe/2));
ak = filter(filt,1,a);
figure ;
plot(t(1:100000),ak(1:100000))

%bk
b = 2*(signalModule.*sin(2*pi*f0*t));
bk = filter(filt,1,b);
figure;
plot(t(1:100000),bk(1:100000))
%figure;
% tracé de ak en fonction de bk
plot(ak(1:100000),bk(1:100000))

%3
%prendre un seul point pour chaque symbole Ak ,Bk
p = 1600; % p = fe/10000
figure;
scatter(ak(p/2:p:N),bk(p/2:p:N),'filled')

%% avec le bruit 

signalModule2 = signalModule + 5*rand(1,N);

% extraction des ak  et bk par demodulation 
%ak
aa = 2*(signalModule2.*cos(2*pi*f0*t));
filt = fir1(100,fc/(fe/2));
akk = filter(filt,1,aa);
figure ;
plot(t(1:100000),ak(1:100000))

%bk
bb = 2*(signalModule2.*sin(2*pi*f0*t));
bkk = filter(filt,1,bb);
figure;
plot(t(1:100000),bk(1:100000))
figure;
% tracé de ak en fonction de bk
plot(akk(1:100000),bkk(1:100000))

%prendre un seul point pour chaque symbole Ak ,Bk
p = 1600;
figure;
scatter(akk(p/2:p:N),bkk(p/2:p:N),'filled')

%% calul de la distance 
A = [1 ,2*sqrt(3/4),max(signalModule)];
S = zeros(8,2);

S(8,:) = [-A(3),A(2)];
S(7,:) = [0,A(2)];
S(6,:) = [A(3),A(2)];
S(5,:) = [-A(1),0];
S(4,:) = [A(1),0];
S(3,:) = [-A(3),A(2)];
S(2,:) = [0,A(2)];
S(1,:) = [A(3),A(2)];

symbole = zeros(1,N/1600);
h = 1 ;

for k = 1 : 1600 : N
    dmin = 1;
    for m = 1 : 1 :8
        d = sqrt((ak(k)- S(m,1))^2 +((bk(k)- S(m,2))^2 ));

        if d <dmin 
            dmin = d ;
            symbole(h) = m;
        end
        Dist(h) = dmin ;
    end
    h= h+1 ;

end

%% decodage 
C = zeros(1,2292/3);
for k = 1: 2292-3
    b = [dec2bin(symbole(1,k),3)  dec2bin(symbole(1,k+1),3)  dec2bin(symbole(1,k+2),3)];
    C(1,k) = bin2dec(b);
    code = char(C);

end