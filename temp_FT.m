function[out]= temp_FT(in)
 out = double(fft(double(in),[],3));
end