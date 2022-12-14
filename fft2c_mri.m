function X=fft2c_mri(x)
% Claire Lin: Otazo's fft
X=fftshift(ifft(ifftshift(x,1),[],1),1)*sqrt(size(x,1));
X=fftshift(ifft(ifftshift(X,2),[],2),2)*sqrt(size(x,2));
% X1=ifft(x,[],1)*sqrt(size(x,1));
% X1=ifft(ifftshift(X1,2),[],2),2)*sqrt(size(x,2));
% rrr=2;