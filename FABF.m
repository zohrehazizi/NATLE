function Xden = FABF(X,rho,N)


[zeta,sigma_r] = logClassifier(X,rho,[23,33]);
Xden = fastABF(X,rho,sigma_r,X+zeta,N);

Rtf = isnan(Xden);
    if Rtf == 1
        Xden = X;
%        Rden = padarray(RdenB4,[2 2]);
%         Rden = Bilateral_Filter(Rden,h,w,25,300,5);
    end




