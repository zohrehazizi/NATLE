t1 = now;

alpha = 0.015;
epsA = 0.001;
epsR = 0.001;
epsG = 0.05;
sigma = 10;
lambdaG = 1.5;
beta = 3;
gamma = 2.2;

addpath('input');
addpath('output');

files = dir('./input/*.bmp');

for i = 1:length(files)
    
    img = im2double(imread(['./input/' files(i).name]));      
    [L,S,Hue,Sat,Rhat,Dx,Dy,sp_eye] = find_L(img,alpha,epsA,epsR);   
    [R,hueDen,satDen] = find_R_denoise(Rhat,Hue,Sat,S,epsG,sigma,lambdaG,Dx,Dy,sp_eye,beta);
    
    L_new = L.^(1/gamma);    
    
    S_cmplx = R.*L_new; 
    S_real = real(S_cmplx);
    max_s = max(max(S_real));
    S_new = S_real/max_s; 
   
    HSV_n = cat(3,hueDen,satDen,S_real);
    RGB = hsv2rgb(double(HSV_n));
    figure, imshow(RGB);
    imwrite(RGB, './output/result.png');
end
    
t2 = now;
dt = t2-t1;
d = datetime(dt,'convertFrom','datenum');
disp(d);



