function Ak = E_forw_for_Ak(xk)
global nx ny nt nc r b1c samp k
s = zeros(nx,ny,nc); 

% smaps_k = squeeze(arg.smaps(:,:,k,:)); %nx ny nc 
samp_k = squeeze(samp(:,:,k)) ;  %nx ny 
s = bsxfun(@times,xk,b1c); % nx ny  nc

S=fft2c_mri(s); % nx ny 1 nc  
Ak = bsxfun(@times,S,samp_k);
end