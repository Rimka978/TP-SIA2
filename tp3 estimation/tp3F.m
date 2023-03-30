[N,T,Z,u,F,G,H,mX0,PX0,Qw,Rv,X] = simulationDonnees;


X_pred = cell(1,N);  
P_pred = cell(1,N);  
X_est  = cell(1,N);  
P_est  = cell(1,N); 
K  = cell(1,N); 
Z_est = cell(1,N); 
Gamma = cell(1,N);  
gamma = zeros(2,N); 


% initialisation :
X_est{1} = mX0;
P_est{1} = PX0;
Gamma{:,1} = [0,0;0,0];

% les equation de KALMAN
for i=2:N
     
    X_pred{i} = F*X_est{i-1} + G*u(:,i-1); % Prediction   %u(2*4);
    P_pred{i} = F*P_est{i-1}*F' + Qw;
    
    K{i} = P_pred{i}*H' * inv(H*P_pred{i}*H' + Rv);
    X_est{i} = X_pred{i} + K{i}*(Z(:,i) - H*X_pred{:,i});
    P_est{i} = P_pred{i} - K{i}*H*P_pred{i};
    
    Gamma{i} = Rv + H*P_pred{i}*H';
    Z_est{i} = H*X_pred{i};
    gamma(:,i) = Z(:,i) - Z_est{i};

end


% Affichage :

%pour x

x = cell2mat(X_est) ;

pX = zeros(1,N);
MX = zeros(1,N);
for k=1:N
    pX(k) = X_est{k}(1)+3*sqrt(P_est{k}(1,1)) ;
    MX(k) = X_est{k}(1)-3*sqrt(P_est{k}(1,1)) ;
end

figure;
plot(X(1,:),'k');hold on
plot(pX,'r'); hold on;
plot(MX,'r');
plot(x(1,:),'b');
title('trajet éstimé sur l''axe x');
legend('X réel','limite superieure', 'limite inferieure','X estimé ');


% pour y
PY = zeros(1,N);
MY = zeros(1,N);
for k=1:N
    PY(k) = X_est{k}(2)+3*sqrt(P_est{k}(2,2)) ;
    MY(k) = X_est{k}(2)-3*sqrt(P_est{k}(2,2)) ;
end

figure;
plot(X(2,:),'g');hold on;
plot(PY,'r'); hold on;
plot(MY,'r');
plot(x(2,:),'b');
title('le trajet éstimé sur l''axe Y')
legend('y réel','limite superieure', 'limite inferieure','Y estimé ')

% Ellipse :

figure,
for i = 1:N

    mx =  x(1:2,i);
        Px = P_est{i}(1:2,1:2);
        ellipse(mx,Px,'b'); hold on
        plot(X(1,i), X(2,i),'r+'); 
        plot(x(1,i), x(2,i),'g+'); 
        legend('Ellipse','Parametre reel','parametre estime')
        title(' ellipse de confiance ')
        %plot(px,PY,'r'), Plot(MX,MY,'r')

% Innovation :

G = Gamma{:,i};

gammax_p(i) = 3*sqrt(G(1,1));
gammax_m(i) = - 3*sqrt(G(1,1));

gammay_p(i) = 3*sqrt(G(2,2));
gammay_m(i) = - 3*sqrt(G(2,2));
    

end
   % plot(pX,PY,'r'), plot(MX,MY,'r')



    
    
%figure,
% plot(gamma(1,:)); hold on
% plot(gammax_p,'m'); hold on
% plot(gammax_m,'m'); 
% 
% figure,
% plot(gamma(2,:)); hold on
% plot(gammay_p,'m'); hold on
% plot(gammay_m,'m'); 
% 
% 
% zest = cell2mat(Z_est);
% figure,
% plot(Z(1,:),'r');hold on ;
% plot(zest(1,:),'k');
% legend('Z reel','Z estimé')
% title('Realisation du Z')


%% 6 innovation
Z_est = cell2mat(Z_est);
figure; 
plot(Z(1,:));
hold on 
plot(Z_est(1,:));
legend("Z","Z est");
title("1");


figure; 
plot(Z(2,:));
hold on 
plot(Z_est(2,:));
legend("Z","Z est");
title("2");


%%
Mz=zeros(1,120);
Pz=zeros(1,120);

for k=1:119

    Pz(1,k)= Z_est(1,k) + 3*sqrt(G(1,1));
    Mz(1,k)= Z_est(1,k) - 3*sqrt(G(1,1));
end


figure;
plot(Z(1,:),'k');hold on
plot(Z_est(1,:),'r'); hold on;
plot(Pz,'c'); hold on;
plot(Mz,'c');
title('Z1 et Z1 estimé');
legend('Z réel','Z estimé ','couloir superieure', 'couloir inferieure');

%%

Mz=zeros(1,120);
Pz=zeros(1,120);

for k=1:119

    Pz(1,k)= Z_est(2,k) + 3*sqrt(G(2,2));
    Mz(1,k)= Z_est(2,k) - 3*sqrt(G(2,2));
end

figure;
plot(Z(2,:),'k');hold on
plot(Z_est(2,:),'r'); hold on;
plot(Pz,'c'); hold on;
plot(Mz,'c');
title('Z2 et Z2 estimé');
legend('Z réel','Z estimé ','couloir limite superieure', 'couloir limite inferieure');



      