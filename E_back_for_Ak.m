
function x = E_back_for_Ak(zk)
global nx ny nt nc r b1c samp k
zk_2 = squeeze(zk); %zk: nx ny nc
x = zeros(nx,ny); 
% smaps_k = squeeze(arg.smaps(:,:,k,:)); %nx ny nc 
samp_k = samp(:,:,k) ;  %nx ny 

 ztmp = bsxfun(@times,zk_2,samp_k); % nx ny nc 
s = ifft2c_mri(ztmp); %nx ny nc 
x = sum(bsxfun(@times,s,conj(b1c)),3);
end