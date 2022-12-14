function B = E_forw_for_AU_new(U_im,y_temp)
global nx ny nt nc r b1c samp mk
B=zeros(r,nt);
n=nx*ny;
smaps = b1c;
mask1=samp;
% U_im = reshape(Uhat,[nx,ny,r]);
U_im2 = reshape(U_im,[ nx, ny, 1,r]);  % nx ny 1 r  
smaps = reshape(smaps, [nx,ny,nc,1]); % nx ny nc 1
AUtmp = bsxfun(@times, U_im2, smaps);  % nx ny nc r
AUtmp1= fft2c_mri(AUtmp); %nx ny nc r
AUtmp2=reshape(AUtmp1,[nx*ny,nc,r]);
    for k=1:1:nt
        ObsFreq=find(mask1(:,:,k));
        m1 = mk(k);
        AUk_tmp = AUtmp2(ObsFreq,:,: ) ;
        yk_tmp=y_temp(ObsFreq,k,: ) ;
        AUk = reshape(AUk_tmp, [m1, r]);
        yk = reshape(yk_tmp, [m1,1]);
        B(:,k)=AUk\yk;
    end
end