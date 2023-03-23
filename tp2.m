clear all;
close all ;
clc ;


load('Données.mat');

%% 3 estimation par moindres carres


M = length(h);
n=0:M-1;
%x1 = [1 zeros(1,length(h) - 1)];
y = syst(x,h,0.3);

% Estimation par moindres carrés 
% 3 
R1 = toeplitz([0;y(1:end-1)],[0  0]);
R2 = toeplitz(x,[x(1) 0]);
R_ = [R1 R2];
teta = (inv(R_'*R_))*R_'*y;  
a1 = teta(1);
a2 = teta(2);
b0 = teta(3);
b1 = teta(4);

% 4
h1= zeros(M,1);       %initialisation de h[n]
h1(1) = b0;            %premier element de h[n];
h1(2) = b1 + a1*b0;    %deuxiement element de h[n];

for i = 3 : M
    h1(i) = a1*h1(i-1) + a2*h1(i-2);
end

% 5
figure(1)
plot(h);
title('Représentation de la réponse impulsionnelle estimée par les moindres carrés');
hold on;
plot(h1,'r');
legend('réponse impulsionnelle originale','réponse impulsionnelle estimée');

d = norm(h - h1)/norm(h); % on trouve une distance d = 0.6591


%6 cette estimation n'est pas du tout satisfaisante

%% methode de steiglitz et Mc Bride 

% initialisation
k = 0;
%y = syst(x,h);
y_tilde = y;
x_tilde = x;
 h2= zeros(M,1);
dis = 5
M = 32;

while dis > 0.3
    R1_ = toeplitz([0;y_tilde(1:end-1)],[0  0]);
   
    R2_ = toeplitz(x_tilde,[x_tilde(1) 0]);
    R_ = [R1_ R2_];
    tetha = (inv(R_'*R_))*R_'*y_tilde;

    a1_ = tetha(1);
    a2_ = tetha(2);
    b0_ = tetha(3);
    b1_ = tetha(4);
    
    G = [1 -a1_ -a2_];

          %initialisation de h[n]
    h2(1) = b0 ;            %premier element de h[n];
    h2(2) = (b1 + a1*b0) ;    %deuxiement element de h[n];


    y_tilde = filter(1,G,y);
    x_tilde = filter(1,G,x);
    
   
    for i = 3 : M
        h2(i) = a1_*h2(i-1) + a2_*h2(i-2);
    end

    dis = norm(h - h2)/norm(h);

    

figure(2);
plot(n,h);
title(sprintf('Réponse Impulsionnelle estimé par Steiglitz et Mc Bride itération %d', k));xlabel('n');ylabel('h(n)'); hold on
plot(n,h2);legend('h','h estimée');
hold off; pause(1)
   k = k+1 ;
end

%%  Biais et variance sur les paramètres estimés
% estimation de biais et de la variance par la méthode de Steiglitz et Mc Bride 

Nr= 400 ;
y_tilde=y;
x_tilde =x;
%declaration des filtres

h3=zeros(length(h),Nr);
 for k=1:Nr
    y_tilde=syst(x,h,0.3);
    x_tilde=x;
     for l=1:Nr
        
            R1 = toeplitz([0;y_tilde(1:end-1)],[0 0]);
            R2 = toeplitz(x_tilde,[x_tilde(1) 0]);
            R = [R1 R2];
            
           
            %calculer les parametres estimes pour Nr
            teeta(:,k)=inv(R'*R)*R'*y_tilde ;
            a1 = teeta(1,k); 
            a2 = teeta(2,k); 
            b0 = teeta(3,k); 
            b1 = teeta(4,k);
    
            h3(1)=b0;
            h3(2)=(b1+a1*b0);

            a=[1 -a1 -a2];
            % Filtrage
            y_tilde=filter(1,a,y);
            x_tilde=filter(1,a,x);
            
            for i=3 :length(h) 
                h3(i)= a1*h3(isa-1)+a2*h3(i-2);
            end
        
            
        end
    

 end

 moyteta= mean(teeta,2);
varteta= var(teeta,0,2);

figure;
subplot(1,2,1); plot(moyteta); title("moyenne paramètres");
subplot(1,2,2); plot(varteta); title("variance paramètres"); 



moyh3= mean(h3,2);
varh3= var(h3,0,2);

figure;
subplot(2,1,1); plot(moyh3); hold on; legend("h estimé", "h theorique"); plot(h); title("moyenne h estimé");
subplot(2,1,2); plot(varh3); title("variance h estimé"); 

biais= moyh3 - h ;
figure; plot(biais);


%% boostrap

%V.2)Estimation du biais et de la variance d'estimation par la methode du Boostrap

%1
yg = y;
xg = x;

for k=1:7
    
    R=[toeplitz([0;yg(1:end-1)],[0 0]) toeplitz(xg,[xg(1) 0])];
 
    teta = (inv(R'*R))*R'*yg;
    a1 = teta(1); 
    a2 = teta(2); 
    b0 = teta(3); 
    b1 = teta(4); 
 
    G = [1 -a1 -a2];
    yg = filter(1,G,y);
    xg = filter(1,G,x);
end

ym = R*teta;
epsilon = y - ym;

figure;
plot(epsilon),title('erreur du modèle');
 

%%
%2 3

he1 = zeros(M,Nr);
teta = zeros(4,1000);
for p=1:1000
    eps_k = epsilon(randi(300,1,300));
    yg = y + eps_k;
    xg = x;
    
    for k=1:7
        R=[toeplitz([0;yg(1:end-1)],[0 0]) toeplitz(xg,[xg(1) 0])];
 
        teta(:,p)=(inv(R'*R))*R'*yg;
        a1 = teta(1,p);
        a2 = teta(2,p);
        b0 = teta(3,p);
        b1 = teta(4,p); 
 
        G = [1 -a1 -a2];

        yg = filter(1,G,yg);
        xg = filter(1,G,xg);

        he1(1,p) = b0;
        he1(2,p) = (b1 + a1*b0);
        G= [1 -a1 -a2];

        for m=3:M
        he1(m,p) = a1*he1(m-1,p)+a2*he1(m-2,p);
        end
    end
end



moyteta1= mean(teta,2);
varteta1= var(teta,0,2);

moyhe1= mean(he1,2);
varhe1= var(he1,0,2);

figure;
subplot(2,1,1); plot(moyhe1); title("moyenne h estimé");
subplot(2,1,2); plot(varhe1); title("variance h estimé"); 




    
    






