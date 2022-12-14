function [Zhat_ktslr,Time_ktslr ]=ktslr_func(kdata,b1c,samp,Xtrue)
tic;
[nx,ny,nt,nc]=size(kdata);
param.E=getE(b1c,nt,'samp',kdata(:,:,:,1)~=0);
b=kdata;

x_init = param.E'*kdata;
mu =20e-6;
opts.mu = mu;
opts.betarate = 35;
opts.outer_iter =12;
opts.inner_iter = 50;
opts.beta=1e3; %

U=x_init; [m n d] = size(U);
Lam= zeros(m,n,d);

o=0;
earray=zeros(1,opts.outer_iter*opts.inner_iter); cost=zeros(1,opts.outer_iter*opts.inner_iter);

fac=length((find(b~=0)));
U = double(U); b=double(b); fac=double(fac);
for out = single(1:opts.outer_iter),
    o=o+1;

    for in = single(1:opts.inner_iter)
        
        z = itemp_FT(U);
        
        
        z1 = (abs(z)-1/opts.beta); z1 = z1.*(z1>0);
        z2 = abs(z) + (abs(z)<1e-12);
        z = z1.*z./z2;
        z = temp_FT(z);
        
        
        
        [U,earray1] = xupdateCG(b,param,z,opts,U,Lam, 1e-7,10);
        
        e = param.E*U - b;
        
        e_ft = itemp_FT(U);
        
        
        cost = [cost, sum(abs(e(:)).^2)  +  sum(abs(e_ft(:)))*opts.mu];
        
        if in>1
            if abs(cost(end) - cost(end-1))/abs(cost(end)) < 1e-3
                break;
            end
        end
        
        
        
    end

    opts.beta=opts.beta*opts.betarate;
end
Zhat_ktslr=U;
Time_ktslr=toc;
end