clear all; close all; clc;

load e1.mat, load e2.mat, load f1.mat, load f2.mat;

%figure, imshow(uint8(e1));
%figure, imshow(uint8(e2));
%figure, imshow(uint8(F1));
%figure, imshow(uint8(F2));

img1 = e1; % img1 = F1
img2 = e2; % img2 = F2

%% Calcul mouvement
b = 16;
f = floor(15 / 2);
[M, N] = size(img1);
mouvement = zeros(ceil(M/b), ceil(N/b), 2);

for i = 1:b:M
    for j = 1:b:N
        block = img1(i:i+b-1, j:j+b-1);
        min_eqm = inf;
        mv = [0 0];
        for u = -f:f
            for v = -f:f
                fi = i+u;
                fj = j+v;
                if (fi >= 1 && fi+b-1 <= M && fj >= 1 && fj+b-1 <= N)
                    block_rech = img2(fi:fi+b-1, fj:fj+b-1);
                    eqm = sum((block(:) - block_rech(:)).^2 / b^2);
                    if eqm < min_eqm
                        min_eqm = eqm;
                        mv = [u, v];
                    end
                end
            end
        end
        mouvement((i-1)/b+1, (j-1)/b+1, :) = mv;
    end
end

%% Reconstruction
img3 = zeros(M, N);
for i = 1:b:M
    for j = 1:b:N
        mv = mouvement((i-1)/b+1, (j-1)/b+1, :);
        x = mv(1);
        y = mv(2);
        img3(i:i+b-1, j:j+b-1) = img1(i+x:i+x+b-1, j+y:j+y+b-1);
    end
end
imshow(uint8(img3));