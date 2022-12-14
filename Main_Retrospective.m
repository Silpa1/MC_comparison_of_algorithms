
clc;clear all;close all;
global nx ny nt nc  b1c samp mk r
[fid,msg] = fopen('Comparison_error.txt','wt');

fprintf(fid, '%s(%s) & %s & %s & %s & %s & %s & %s   \n','Dataset','Radial','ktslr','L+S-Otazo','L+S-Lin','altGDmin_mean','altGDmin-MRI','altGDmin-MRI2');
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kdata=[];
filenames={ 'cardiac_cine_R6.mat'}
load('Xinf.mat')
Xtrue = double(Xinf.cine);
load(filenames{1})
[~,name,~] = fileparts(filenames{1});
[nx,ny,nt,nc]=size(kdata);
n=nx*ny;



tmp = sqrt(sum(abs((b1)).^2,3));
b1c = b1./tmp;


samp = kdata(:,:,:,1)~=0;
param.E=getE(b1c,nt,'samp',samp);

radial=0;
[Xhat_ktslr,Time_ktslr ]=ktslr_func(kdata,b1c,samp);
NMSE_ktslr=RMSE_modi(Xhat_ktslr,Xtrue);
[Xhat_LpS_otazo,Time_otazo ]=lpluss_otazo_func(kdata,b1c,samp);
NMSE_otazo=RMSE_modi(Xhat_LpS_otazo,Xtrue);
[Xhat_LpS_lin,Time_lin]=lpluss_lin_func(kdata,b1c,samp);
NMSE_lin=RMSE_modi(Xhat_LpS_lin,Xtrue);
[Xhat_mri,Time_mri]=altGDmin_mri1_func(kdata,b1c,samp);
NMSE_mri=RMSE_modi(Xhat_mri,Xtrue);
[Xhat_mri2,Time_mri2]=altGDmin_mri2_func(kdata,b1c,samp);
NMSE_mri2=RMSE_modi(Xhat_mri2,Xtrue);
fprintf(fid, '%s(%d) & %8.4f (%5.2f)&%8.4f (%5.2f)& %8.4f (%5.2f)&%8.4f (%5.2f)& %8.4f (%5.2f) \n', name, radial,NMSE_ktslr,Time_ktslr,NMSE_otazo,Time_otazo,NMSE_lin,Time_lin,NMSE_mri,Time_mri,NMSE_mri2,Time_mri2);

% % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

kdata=[];
filenames={ 'cardiac_perf_R8.mat'}
load('Xinf.mat')
Xtrue = double(Xinf.perf);

load(filenames{1})
[~,name,~] = fileparts(filenames{1});
[nx,ny,nt,nc]=size(kdata);
n=nx*ny;

tmp = sqrt(sum(abs((b1)).^2,3));
b1c = div0(b1,tmp);


samp = kdata(:,:,:,1)~=0;
param.E=getE(b1c,nt,'samp',samp);

radial=0;
[Xhat_ktslr,Time_ktslr ]=ktslr_func(kdata,b1c,samp);
NMSE_ktslr=RMSE_modi(Xhat_ktslr,Xtrue);
[Xhat_LpS_otazo,Time_otazo ]=lpluss_otazo_func(kdata,b1c,samp);
NMSE_otazo=RMSE_modi(Xhat_LpS_otazo,Xtrue);
[Xhat_LpS_lin,Time_lin]=lpluss_lin_func(kdata,b1c,samp);
NMSE_lin=RMSE_modi(Xhat_LpS_lin,Xtrue);
[Xhat_mri,Time_mri]=altGDmin_mri1_func(kdata,b1c,samp);
NMSE_mri=RMSE_modi(Xhat_mri,Xtrue);
[Xhat_mri2,Time_mri2]=altGDmin_mri2_func(kdata,b1c,samp);
NMSE_mri2=RMSE_modi(Xhat_mri2,Xtrue);
fprintf(fid, '%s(%d) & %8.4f (%5.2f)&%8.4f (%5.2f)& %8.4f (%5.2f)&%8.4f (%5.2f)& %8.4f (%5.2f) \n', name, radial,NMSE_ktslr,Time_ktslr,NMSE_otazo,Time_otazo,NMSE_lin,Time_lin,NMSE_mri,Time_mri,NMSE_mri2,Time_mri2);


% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kdata=[];
filenames={ 'brain_T2_T1rho.mat'}
load(filenames{1})
[~,name,~] = fileparts(filenames{1});
[nx,ny,nc,nt]=size(img1);
n=nx*ny;
k= (1/sqrt(n))*fft2(img1);
no_comp=8;
[k]= coil_compress_withpca(k,no_comp);
[nx,ny,nc,nt]=size(k);
im = sqrt(n)*ifft2(k);
xm = mean(im,4);
b1 = ismrm_estimate_csm_walsh_modified(xm);
Xtrue = squeeze(sum(im.*repmat(conj(b1),[1 1 1 nt]),3));
for i =1:1:nc
    kdata(:,:,:,i)=k(:,:,i,:);
end
tmp = sqrt(sum(abs((b1)).^2,3));
b1c = div0(b1,tmp);
radial=[4,8,16];

for ii=1:1:length(radial)
    samp = goldencart(nx,ny,nt,radial(ii));
    mask=repmat(samp,[1,1,1,nc]);
    param.E=getE(b1c,nt,'samp',mask(:,:,:,1)~=0);
    kdata=param.E*Xtrue;
    [Xhat_ktslr,Time_ktslr ]=ktslr_func(kdata,b1c,samp);
    NMSE_ktslr=RMSE_modi(Xhat_ktslr,Xtrue);
    [Xhat_LpS_otazo,Time_otazo ]=lpluss_otazo_func(kdata,b1c,samp);
    NMSE_otazo=RMSE_modi(Xhat_LpS_otazo,Xtrue);
    [Xhat_LpS_lin,Time_lin]=lpluss_lin_func(kdata,b1c,samp);
    NMSE_lin=RMSE_modi(Xhat_LpS_lin,Xtrue);
    [Xhat_mri,Time_mri]=altGDmin_mri1_func(kdata,b1c,samp);
    NMSE_mri=RMSE_modi(Xhat_mri,Xtrue);
    [Xhat_mri2,Time_mri2]=altGDmin_mri2_func(kdata,b1c,samp);
    NMSE_mri2=RMSE_modi(Xhat_mri2,Xtrue);
    fprintf(fid, '%s(%d) & %8.4f (%5.2f)&%8.4f (%5.2f)& %8.4f (%5.2f)&%8.4f (%5.2f)& %8.4f (%5.2f) \n', name, radial(ii),NMSE_ktslr,Time_ktslr,NMSE_otazo,Time_otazo,NMSE_lin,Time_lin,NMSE_mri,Time_mri,NMSE_mri2,Time_mri2);
end
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kdata=[];
filenames={ 'multi_coil_lowres_speech.mat'}
load(filenames{1})
[~,name,~] = fileparts(filenames{1});
no_comp=8;
[k]= coil_compress_withpca(k,no_comp);
[nx,ny,nc,nt]=size(k);
n=nx*ny;
im = sqrt(n)*ifft2(k);
xm = mean(im,4);
b1 = ismrm_estimate_csm_walsh_modified(xm);
Xtrue = squeeze(sum(im.*repmat(conj(b1),[1 1 1 nt]),3));
tmp = sqrt(sum(abs((b1)).^2,3));
b1c = div0(b1,tmp);

radial=[4,8,16];
for ii=1:1:length(radial)
    samp = goldencart(nx,ny,nt,radial(ii));
    mask=repmat(samp,[1,1,1,nc]);
    param.E=getE(b1c,nt,'samp',mask(:,:,:,1)~=0);
    kdata=param.E*Xtrue;
    [Xhat_ktslr,Time_ktslr ]=ktslr_func(kdata,b1c,samp);
    NMSE_ktslr=RMSE_modi(Xhat_ktslr,Xtrue);
    [Xhat_LpS_otazo,Time_otazo ]=lpluss_otazo_func(kdata,b1c,samp);
    NMSE_otazo=RMSE_modi(Xhat_LpS_otazo,Xtrue);
    [Xhat_LpS_lin,Time_lin]=lpluss_lin_func(kdata,b1c,samp);
    NMSE_lin=RMSE_modi(Xhat_LpS_lin,Xtrue);
    [Xhat_mri,Time_mri]=altGDmin_mri1_func(kdata,b1c,samp);
    NMSE_mri=RMSE_modi(Xhat_mri,Xtrue);
    [Xhat_mri2,Time_mri2]=altGDmin_mri2_func(kdata,b1c,samp);
    NMSE_mri2=RMSE_modi(Xhat_mri2,Xtrue);
    fprintf(fid, '%s(%d) & %8.4f (%5.2f)&%8.4f (%5.2f)& %8.4f (%5.2f)&%8.4f (%5.2f)& %8.4f (%5.2f) \n', name, radial(ii),NMSE_ktslr,Time_ktslr,NMSE_otazo,Time_otazo,NMSE_lin,Time_lin,NMSE_mri,Time_mri,NMSE_mri2,Time_mri2);
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kdata=[];
filenames={ 'multicoil_ungated_cmr_perf.mat'}
load(filenames{1})
[~,name,~] = fileparts(filenames{1});
no_comp=8;
[kSpace]= coil_compress_withpca(kSpace,no_comp);
[nx,ny,nc,nt]=size(kSpace);
n=nx*ny;
im = sqrt(n)*ifft2(kSpace);
xm = mean(im,4);
b1= ismrm_estimate_csm_walsh_modified(xm);
Xtrue = squeeze(sum(im.*repmat(conj(b1),[1 1 1 nt]),3));
tmp = sqrt(sum(abs((b1)).^2,3));
b1c = div0(b1,tmp);
radial=[4];
for ii=1:1:length(radial)
    samp = goldencart(nx,ny,nt,radial(ii));
    mask=repmat(samp,[1,1,1,nc]);
    param.E=getE(b1c,nt,'samp',mask(:,:,:,1)~=0);
    kdata=param.E*Xtrue;
    [Xhat_ktslr,Time_ktslr ]=ktslr_func(kdata,b1c,samp);
    NMSE_ktslr=RMSE_modi(Xhat_ktslr,Xtrue);
    [Xhat_LpS_otazo,Time_otazo ]=lpluss_otazo_func(kdata,b1c,samp);
    NMSE_otazo=RMSE_modi(Xhat_LpS_otazo,Xtrue);
    [Xhat_LpS_lin,Time_lin]=lpluss_lin_func(kdata,b1c,samp);
    NMSE_lin=RMSE_modi(Xhat_LpS_lin,Xtrue);
    [Xhat_mri,Time_mri]=altGDmin_mri1_func(kdata,b1c,samp);
    NMSE_mri=RMSE_modi(Xhat_mri,Xtrue);
    [Xhat_mri2,Time_mri2]=altGDmin_mri2_func(kdata,b1c,samp);
    NMSE_mri2=RMSE_modi(Xhat_mri2,Xtrue);
    fprintf(fid, '%s(%d) & %8.4f (%5.2f)&%8.4f (%5.2f)& %8.4f (%5.2f)&%8.4f (%5.2f)& %8.4f (%5.2f) \n', name, radial(ii),NMSE_ktslr,Time_ktslr,NMSE_otazo,Time_otazo,NMSE_lin,Time_lin,NMSE_mri,Time_mri,NMSE_mri2,Time_mri2);
    
end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kdata=[];
filenames={ 'fs_0016_3T_ksp_csm.mat'}
load(filenames{1})
[~,name,~] = fileparts(filenames{1});
k = permute(kspace,[3 4 2 1]);
[nx,ny,nc,nt]=size(k);
for c = 1:nc,
    for t = 1:nt,
        k(:,:,c,t) = fft2(ifftshift(ifft2(fftshift(k(:,:,c,t)))));
    end
end
no_comp=8;
[k]= coil_compress_withpca(k,no_comp);
n=nx*ny;
im = sqrt(n)*ifft2(k);
xm = mean(im,4);
b1 = ismrm_estimate_csm_walsh_modified(xm);
Xtrue= squeeze(sum(im.*repmat(conj(b1),[1 1 1 nt]),3));
tmp = sqrt(sum(abs((b1)).^2,3));
b1c = div0(b1,tmp);
radial=[4,8,16];
for ii=1:1:length(radial)
    samp = goldencart(nx,ny,nt,radial(ii));
    mask=repmat(samp,[1,1,1,nc]);
    param.E=getE(b1c,nt,'samp',mask(:,:,:,1)~=0);
    kdata=param.E*Xtrue;
    [Xhat_ktslr,Time_ktslr ]=ktslr_func(kdata,b1c,samp);
    NMSE_ktslr=RMSE_modi(Xhat_ktslr,Xtrue);
    [Xhat_LpS_otazo,Time_otazo ]=lpluss_otazo_func(kdata,b1c,samp);
    NMSE_otazo=RMSE_modi(Xhat_LpS_otazo,Xtrue);
    [Xhat_LpS_lin,Time_lin]=lpluss_lin_func(kdata,b1c,samp);
    NMSE_lin=RMSE_modi(Xhat_LpS_lin,Xtrue);
    [Xhat_mri,Time_mri]=altGDmin_mri1_func(kdata,b1c,samp);
    NMSE_mri=RMSE_modi(Xhat_mri,Xtrue);
    [Xhat_mri2,Time_mri2]=altGDmin_mri2_func(kdata,b1c,samp);
    NMSE_mri2=RMSE_modi(Xhat_mri2,Xtrue);
    fprintf(fid, '%s(%d) & %8.4f (%5.2f)&%8.4f (%5.2f)& %8.4f (%5.2f)&%8.4f (%5.2f)& %8.4f (%5.2f) \n', name, radial(ii),NMSE_ktslr,Time_ktslr,NMSE_otazo,Time_otazo,NMSE_lin,Time_lin,NMSE_mri,Time_mri,NMSE_mri2,Time_mri2);
    
end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kdata=[];
filenames={ 'fs_0019_3T_ksp_csm.mat'}
load(filenames{1})
[~,name,~] = fileparts(filenames{1});
k = permute(kspace,[3 4 2 1]);
[nx,ny,nc,nt]=size(k);
for c = 1:nc,
    for t = 1:nt,
        k(:,:,c,t) = fft2(ifftshift(ifft2(fftshift(k(:,:,c,t)))));
    end
end
no_comp=8;
[k]= coil_compress_withpca(k,no_comp);
[nx,ny,nc,nt]=size(k);
n=nx*ny;
im = sqrt(n)*ifft2(k);
xm = mean(im,4);
b1 = ismrm_estimate_csm_walsh_modified(xm);
Xtrue= squeeze(sum(im.*repmat(conj(b1),[1 1 1 nt]),3));
tmp = sqrt(sum(abs((b1)).^2,3));
b1c = div0(b1,tmp);

radial=[4,8,16];
for ii=1:1:length(radial)
    samp = goldencart(nx,ny,nt,radial(ii));
    mask=repmat(samp,[1,1,1,nc]);
    param.E=getE(b1c,nt,'samp',mask(:,:,:,1)~=0);
    kdata=param.E*Xtrue;
    [Xhat_ktslr,Time_ktslr ]=ktslr_func(kdata,b1c,samp);
    NMSE_ktslr=RMSE_modi(Xhat_ktslr,Xtrue);
    [Xhat_LpS_otazo,Time_otazo ]=lpluss_otazo_func(kdata,b1c,samp);
    NMSE_otazo=RMSE_modi(Xhat_LpS_otazo,Xtrue);
    [Xhat_LpS_lin,Time_lin]=lpluss_lin_func(kdata,b1c,samp);
    NMSE_lin=RMSE_modi(Xhat_LpS_lin,Xtrue);
    [Xhat_mri,Time_mri]=altGDmin_mri1_func(kdata,b1c,samp);
    NMSE_mri=RMSE_modi(Xhat_mri,Xtrue);
    [Xhat_mri2,Time_mri2]=altGDmin_mri2_func(kdata,b1c,samp);
    NMSE_mri2=RMSE_modi(Xhat_mri2,Xtrue);
    fprintf(fid, '%s(%d) & %8.4f (%5.2f)&%8.4f (%5.2f)& %8.4f (%5.2f)&%8.4f (%5.2f)& %8.4f (%5.2f) \n', name, radial(ii),NMSE_ktslr,Time_ktslr,NMSE_otazo,Time_otazo,NMSE_lin,Time_lin,NMSE_mri,Time_mri,NMSE_mri2,Time_mri2);
    
end

fclose(fid);