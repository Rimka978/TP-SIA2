close all ;
clc;

%image formes
Img0=imread('formes.pgm');
im1=double(Img0) ;
% 
% %figure;
% subplot(1,2,1);imagesc(im1);
% subplot(1,2,2);imhist(uint8(im1));
% colormap gray 
% axis([0 300 0 20000])

nc=5;

%choix aléatoire des classes


iclasse=zeros(1,nc);
for k=1:nc
    iclasse(1,k)= randi([round((k-1)*(255/nc)),round((255/nc)*k)],1,1);
end

[c,l]=size(im1);
nclasse= zeros(c,l);
nbP= zeros(1,nc);
somme= zeros(1,nc);
Nclasse=zeros(1,nc);
t=0;

nb=0;


while(t==0)
    nclasse= zeros(c,l);
    nbP= zeros(1,nc);
    somme= zeros(1,nc);subplot(1,2,2);imagesc(im1); title('image initiale');
    for k=1:c
        for p=1:l
            d= abs(im1(k,p)*ones(1,nc) - iclasse);
            [val,I] = min(d);
            nclasse(k,p)=I;

            nbP(1,I)= nbP(1,I)+1 ;
            somme(1,I)= somme(1,I)+ im1(k,p);

        end
    end
    nbP= zeros(1,nc);
    somme= zeros(1,nc);
    for k=1:c
        for p=1:l
            d= abs(im1(k,p)*ones(1,nc) - iclasse);
            [val,I] = min(d);
            nclasse(k,p)=I;

            nbP(1,I)= nbP(1,I)+1 ;
            somme(1,I)= somme(1,I)+ im1(k,p);

        end
    end
    Nclasse=round(somme./nbP) ;
    T = isnan(Nclasse);
    Nclasse(isnan(Nclasse))= iclasse(1,find(T==1));
  
    t=isequal(Nclasse,iclasse);
    iclasse=Nclasse;
    nb=nb+1;
end

imG=nclasse;
for k=1:nc
     imG(imG(:)==k)= iclasse(1,k);
end

figure;
subplot(1,2,1);imagesc(imG); title('image reconstruite');
subplot(1,2,2);imagesc(im1); title('image initiale');
colormap gray
figure;
subplot(1,2,1);imhist(uint8(imG)); title('hist image reconstruite'); 
subplot(1,2,2);imhist(uint8(im1));  title('hist image initiale');



%%
%%
%Image brain

 IM1=imread('brain.pgm');
 im1=double(IM1);
 [c1,l1]=size(im1);
 im1=double(IM1) + 3*randn(c1,l1);

figure;
subplot(1,2,1);imagesc(im1);
subplot(1,2,2);imhist(IM1);
colormap gray 

nc=5;


%choix aléatoire mais avec valeurs bien espacées
iclasse=zeros(1,nc);
F=round(230/nc);
for k=1:nc
    iclasse(1,k)= randi([(k-1)*F,F*k],1,1);
end

nclasse= zeros(c1,l1);
nbP= zeros(1,nc);
somme= zeros(1,nc);
nwclasse=zeros(1,nc);
tf=0;

nb=0;

while(tf==0)
    nclasse= zeros(c1,l1);
    nbP= zeros(1,nc);
    somme= zeros(1,nc);
    for k=1:c1
        for p=1:l1
            d= abs(im1(k,p)*ones(1,nc) - iclasse);
            [val,I] = min(d);
            nclasse(k,p)=I;

            nbP(1,I)= nbP(1,I)+1 ;
            somme(1,I)= somme(1,I)+ im1(k,p);

        end
    end

    nwclasse=round(somme./nbP) ;
    TF = isnan(nwclasse);
    nwclasse(isnan(nwclasse))= iclasse(1,find(TF==1));
  
    tf=isequal(nwclasse,iclasse);
    iclasse=nwclasse;
    nb=nb+1;
end

imG=nclasse;
for k=1:nc
     imG(imG(:)==k)= iclasse(1,k);
end

figure;
imagesc(imG); 

colormap gray
figure;
imhist(uint8(imG)); 


%%

IM1=imread('carre_bruit.png');
im1=double(IM1); %+ rand(256);
 
figure;
subplot(1,2,1);imagesc(IM1);
subplot(1,2,2);imhist(IM1);


red=IM1(:,:,1);
green=IM1(:,:,2);
blue=IM1(:,:,3);

figure;
subplot(1,3,1);imhist(red); title('Red');
subplot(1,3,2);imhist(green); title('Green');
subplot(1,3,3);imhist(blue); title('Blue');

nc=3;

iclasse= randi([0,255],3,nc); %choix aléatoire de classes  RGB

[c1,l1,j]=size(im1);
nclasse= zeros(c1,l1,3);
nbP= zeros(3,nc);
somme= zeros(3,nc);
nwclasse=zeros(3,nc);
tf=0;

nb=0;
d=zeros(3,nc);
while(tf==0)
    nclasse= zeros(c1,l1,3);
    nbP= zeros(3,nc);
    somme= zeros(3,nc);
    for m=1:3
        for k=1:c1
                d(m,:)= abs(im1(k,p,m)*ones(1,nc) - iclasse(m,:));
            for p=1:l1
                d(m,:)= abs(im1(k,p,m)*ones(1,nc) - iclasse(m,:));
               
                [val,I] = min(d');
                nclasse(k,p,m)=I(1,m);
    
                nbP(m,I(1,m))= nbP(m,I(1,m))+1 ;
    
                somme(m,I(1,m))= somme(m,I(1,m))+ im1(k,p,m);
    
            end
        end
    end
    nwclasse=round(somme./nbP) ;
    TF = isnan(nwclasse);
    nwclasse(isnan(nwclasse))= iclasse(find(TF==1));
  
    tf=isequal(nwclasse,iclasse);
    iclasse=nwclasse;
    nb=nb+1;
                d(m,:)= abs(im1(k,p,m)*ones(1,nc) - iclasse(m,:));
end

imG=im1;
for m=1:3
    imC=nclasse(:,:,m);
    for k=1:nc
         imC(imC(:,:)==k)= iclasse(m,k);
    end
    imG(:,:,m)=imC(:,:);
end

figure;
imagesc(uint8(imG));

colormap gray
figure;
imhist(uint8(imG)); 

