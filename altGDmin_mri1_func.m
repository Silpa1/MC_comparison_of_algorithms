function [X_MRI,Time_MRI]=altGDmin_mri1_func(kdata,b1c,samp,Xtrue)
global r mk
tic;
y=kdata;
[nx,ny,nt,nc]=size(kdata);
param.E=getE(b1c,nt,'samp',kdata(:,:,:,1)~=0);
n=nx*ny;
[zbar_hat,flag,resNE,iter] = cgls_mean(@E_forw_for_mean,@E_back_for_mean, y,0,1e-3,10);
ybar_hat=param.E*repmat(zbar_hat,[1,1,nt]);
 yinter=y-ybar_hat;
sum_mk=nnz(y);
for k=1:1:nt
    mk(k)=nnz(y(:,:,k,:));
end
m=max(mk);
C_tilda=36;
alpha=C_tilda*norm(yinter(:))^2/sum_mk;
Y_trunc=yinter;
Y_trunc(abs(yinter)>sqrt(alpha))=0;
X0_temp=param.E'*Y_trunc;
DiagMat=diag(ones(1,nt)./sqrt(mk*m));
X0_image=X0_temp;
X0=(reshape(X0_image,[nx*ny,nt]))*(DiagMat);
r_big=floor(min([n/10,nt/10,m/10]));
[Utemp, Stemp,~]=svds(X0,r_big);
SS=diag(Stemp);
% tic;
E=sum(SS.^2);
Esum=0;
for i=1:1:r_big
    Esum=Esum+((SS(i))^2);
    if Esum >(E*0.85)
        break
    end
end
r=i+1;
r=min(r,r_big);
U0=Utemp(:,1:r);
Uhat=U0;
T=70;
  y_temp=reshape(yinter,[nx*ny,nt,nc]);

for t = 1 : T
  
    Uhatm=reshape(Uhat,[nx,ny,r]);
    B = E_forw_for_AU_new(Uhatm,y_temp);
    X=reshape(Uhat*B,[nx,ny,nt]);
    Z=param.E'*((param.E*X)-yinter);
    Z_mat=reshape(Z,[nx*ny,nt]);
    Grad_U=Z_mat*B';
    if t==1
        eta=1/(7*norm(Grad_U));
    end
    Uhat_t0=Uhat;
    Uhat=Uhat-eta*Grad_U;
    [Qu,~]  =  qr(Uhat,0);
    Uhat  =  Qu(:, 1:r);
    Uhat_t1=Uhat;
    Subspace_d= ( norm((Uhat_t0 - Uhat_t1*(Uhat_t1'*Uhat_t0)), 'fro')/sqrt(r));
    if  (Subspace_d <=.01)
        break;
    end
 
end
X_GD=X+zbar_hat;
yk=y-(param.E*X_GD);
global k

for k=1:nt
    Ehat(:,:,k)=cgls_modi(@E_forw_for_Ak,@E_back_for_Ak,squeeze(yk(:,:,k,:)),0,1e-36,3);
end
Time_MRI=toc;
X_MRI=X+zbar_hat+Ehat;
end