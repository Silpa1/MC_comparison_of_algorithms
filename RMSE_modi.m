function[Error,U]= RMSE_modi(U,I)
U = (U); I = (I);
E=0;E1=0;
for ii=1:size(U,3)
alpha = sum(dot(U(:,:,ii),I(:,:,ii)))/(sum(dot(U(:,:,ii),U(:,:,ii))));

U(:,:,ii)=(alpha)*U(:,:,ii);

E=E+sum(sum(abs((I(:,:,ii)-U(:,:,ii))).^2));         
E1=E1+(sum(sum(abs((I(:,:,ii))).^2)));

end
E2=E*1/E1;

Error=(E2);

