function [L,S,Hue,Sat,Rhat,Dx,Dy,sp_eye] = find_L(Img,alpha,epsA,epsR)

HSV = rgb2hsv(Img);                                %convert to hsv 
disp('converted to hsv');
S = HSV(:,:,3);                                    %v component
%S = imresize(S,0.5);
Hue = HSV(:,:,1);
%Hue = imresize(Hue,0.5);
Sat = HSV(:,:,2);
%Sat = imresize(Sat,0.5);
Lhat = (0.2989*Img(:,:,1)+0.5870*Img(:,:,2)+0.1140*Img(:,:,3));%/3;       %L estimation

[grad_Lhat_x, grad_Lhat_y] = find_grad(Lhat);        %L hat gradients
disp('L^ grad found');

A_x = (alpha)./(abs(grad_Lhat_x)+epsA);    %Ad(x) 
A_y = (alpha)./(abs(grad_Lhat_y)+epsA);

disp('A computed');

[h,w] = size(Lhat);         %vectorizing
hw = h*w;
Lhat_vec = reshape(Lhat,[hw,1]);

A_xvec = reshape(A_x,[hw,1])';       %diag(ad)
A_yvec = reshape(A_y,[hw,1])';
i_ind = [1:hw];
diag_ax = sparse(i_ind,i_ind,A_xvec,hw,hw);
diag_ay = sparse(i_ind,i_ind,A_yvec,hw,hw);


D_val = repmat([-1,1],1,hw);        %Dy
Dy_indx = [];
for i = 1:w
    Dy_indx = [Dy_indx,(i-1)*h+1,(i-1)*h+2];
    for j = 1:h-2
        Dy_indx = [Dy_indx,(i-1)*h+j,(i-1)*h+j+2];
    end
    Dy_indx = [Dy_indx,i*h-1,i*h];
end
D_indy = [];
for i = 1:hw
    D_indy = [D_indy,i,i];
end
Dy = sparse(D_indy,Dy_indx,D_val,hw,hw,2*hw);


disp('Dy computed');

Dx_indx = [];
for j = 1:h
    Dx_indx = [Dx_indx,j,h+j];
end
for j = 1:(w-2)*h
    Dx_indx = [Dx_indx,j,2*h+j];
end
for j = ((w-2)*h+1):(hw-h)
    Dx_indx = [Dx_indx,j,h+j];
end
Dx = sparse(D_indy,Dx_indx,D_val,hw,hw,2*hw);



disp('Dx computed');

vec3 = ones(hw,1)';
ind_i = [1:hw];
sp_eye = sparse(ind_i,ind_i,vec3,hw,hw);

A = sp_eye+Dx'*diag_ax*Dx+Dy'*diag_ay*Dy;
b = Lhat_vec;
L_vec = bicgstabl(A,b,[],100);

L = reshape(L_vec,[h,w]);
disp('L computed');

L(L==0) = epsR;
Rhat = S./(L);                                              % R^
disp('R^ computed');
