function [Xhat_LpS,Time_Otazo]=lpluss_otazo_func(kdata,b1c,samp,Xtrue)
tic;
[nx,ny,nt,nc]=size(kdata);
param.E=getE(b1c,nt,'samp',kdata(:,:,:,1)~=0);
param.d=kdata;
param.T=TempFFT(3);
param.lambda_L=0.01;
param.lambda_S=0.01;
param.nite=50;
param.tol=0.0025;
M=param.E'*param.d;
M=reshape(M,[nx*ny,nt]);
Lpre=M;
S=zeros(nx*ny,nt);
ite=0;
while(1)
    ite=ite+1;

    M0=M;
    [Ut,St,Vt]=svd(M-S,0);
    St=diag(SoftThresh(diag(St),St(1)*param.lambda_L));
    L=Ut*St*Vt';
    S=reshape(param.T'*(SoftThresh(param.T*reshape(M-Lpre,[nx,ny,nt]),param.lambda_S)),[nx*ny,nt]);
    resk=param.E*reshape(L+S,[nx,ny,nt])-param.d;
    M=L+S-reshape(param.E'*resk,[nx*ny,nt]);
    Lpre=L;
    tmp2=param.T*reshape(S,[nx,ny,nt]);
     if (ite > param.nite) || (norm(M(:)-M0(:))<param.tol*norm(M0(:))), break;end
    Xhat_LpS1=L+S;
    Xhat_LpS=reshape(Xhat_LpS1,[nx,ny,nt]);
%     Time1_Otazo(ite)= toc;
%     Error_Otazo(ite)=RMSE_modi(Xhat_LpS,Xtrue);
end
Xhat_LpS1=L+S;
Xhat_LpS=reshape(Xhat_LpS1,[nx,ny,nt]);
Time_Otazo= toc;
end