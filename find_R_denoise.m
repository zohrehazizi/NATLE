function [R,hueDen,satDen] = find_R_denoise(Rhat,Hue,Sat,S,epsG,sigma,lambdaG,Dx,Dy,sp_eye,beta)

%[h,w] = size(Rhat);
Rhat_hsv = cat(3,Hue,Sat,Rhat);
Rhat_rgb = hsv2rgb(Rhat_hsv);

Rden = Rhat_rgb(:,:,1);
Gden = Rhat_rgb(:,:,2);
Bden = Rhat_rgb(:,:,3);


RdenB4 = medfilt2(Rden);                          % Median Filter Denoising 
GdenB4 = medfilt2(Gden);
BdenB4 = medfilt2(Bden);

disp('Median Filter done');

N = 1;
rho = 5;
Rden = FABF(RdenB4,rho,N);
Gden = FABF(GdenB4,rho,N);
Bden = FABF(BdenB4,rho,N);


rgbnew = cat(3,Rden,Gden,Bden);
hsvnew = rgb2hsv(rgbnew);
hueDen = hsvnew(:,:,1);
satDen = hsvnew(:,:,2);
RhatDen = hsvnew(:,:,3);

[h,w] = size(RhatDen);                                              % vectorizing
hw = h*w;
RhatDenVec = reshape(RhatDen,[hw,1]);
disp('Rhat Denoised vector');

[Gx,Gy] = find_G(S,epsG,sigma,lambdaG);                             % G
%G = abs(find_grad(S));
[h,w] = size(Gx);                                              % vectorizing
hw = h*w;
Gx_vec = reshape(Gx,[hw,1]);
Gy_vec = reshape(Gy,[hw,1]);
disp('G vect');

A = sp_eye+beta*(Dx'*Dx+Dy'*Dy);
b = RhatDenVec+beta*(Dx'*Gx_vec+Dy'*Gy_vec);
%b = RhatDenVec+beta*(Dx'+Dy')*G_vec;
R_vec = bicgstabl(A,b,[],100);

R = reshape(R_vec,[h,w]);
disp('R computed');

     









