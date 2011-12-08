function I_pm=perona_malik(I,iter)
% This function implement for Perona,Malik,Nonlinear Diffusion
% 
lambda=1/4;
I_pm=I;
[xpm,ypm]=size(I_pm);
for i=1:iter
    idxx=2:xpm-1;
    idxy=2:ypm-1;
    I_n=I_pm(idxx-1,idxy)-I_pm(idxx,idxy);
    I_s=I_pm(idxx+1,idxy)-I_pm(idxx,idxy);
    I_e=I_pm(idxx,idxy-1)-I_pm(idxx,idxy);
    I_w=I_pm(idxx,idxy+1)-I_pm(idxx,idxy);
    max_n=max(I_n(:));
    max_s=max(I_s(:));
    max_e=max(I_e(:));
    max_w=max(I_w(:));
    k=max([max_n,max_s,max_e,max_w]);
    k=0.05*k;
    C_n=1./(1.+(abs(I_n))./k);
    C_s=1./(1.+(abs(I_s))./k);
    C_e=1./(1.+(abs(I_e))./k);
    C_w=1./(1.+(abs(I_w))./k);
    I_pm(idxx,idxy)=I_pm(idxx,idxy)+lambda.*(C_n.*I_n+...
                    C_s.*I_s+C_e.*I_e+C_w.*I_w);
end