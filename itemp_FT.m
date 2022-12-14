function[out]= itemp_FT(in)
 out = double(ifft(double(in),[],3));
end