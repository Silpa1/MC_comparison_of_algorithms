% CG solution to region limited problem
% [X,Potential,ErrorArray,ErrorIndex,efinal] = xupdateHOTV(A,b,baseline,mask,kappa,lambda,mu,Niter,Potential)
% Solves {X*} = arg min_{X} ||Af-b||^2 + mu ||Rf||_{l_1}
%--------------------------------------------------------------------------

function [X,earray1] = xupdateCG(b,param,z,C,X, Lam,THRESHOLD,Niter)


oldcost = 0;
earray1 = [];
lam1 = 0.5*C.mu*C.beta;
el1=double(0);
for i=1:Niter,
    
    resY = (param.E*X - b);
    eY = sum(abs(resY(:)).^2);
    
    resl1 =(X)-(z); 
    el1 = lam1*sum(abs(resl1(:)).^2);
     
    
    cost1 = eY + el1;
    
    earray1 = [earray1,cost1];
    
    if(abs(cost1-oldcost)/abs(cost1) < THRESHOLD)
       % i
        break;
    end
    oldcost = cost1;
    
  %  conjugate gradient direction
   % ------------------------------
    
    % gradient: gn
    
    gn = 2*param.E'*resY + 2*lam1*(resl1) ;%+ C.mu*Lam ;
   
    
    % search direction: sn  
    if(i==1)
        sn = gn;                                          
        oldgn = gn;
    else
        gamma = abs(sum(gn(:)'*gn(:))/sum(oldgn(:)'*oldgn(:)));
        sn = gn + gamma*sn; 
        oldgn = gn;
    end
    
    % line search
    %-------------
    Asn =param.E*sn;  
   
    
    %numer = Asn(:)'*resY(:) + lam1*sn(:)'*resl1(:) ;%+0.5* C.mu*sn(:)'*Lam(:) ;
    numer = sum(conj(Asn(:)).*resY(:)) + lam1*(sum(conj(sn(:)).*resl1(:)));% + 0.5* C.mu1*sum(conj(sn(:)).*Lam4(:));
  
    denom = sum(conj(Asn(:)).*Asn(:)) + lam1*sum(conj(sn(:)).*sn(:)); 
   
    
   % denom = Asn(:)'*Asn(:) + lam1*sn(:)'*sn(:); 
    
    if(denom < 1e-18)
       % break;
    end
    alpha = -real(numer)/real(denom);
   
    % updating
    %-------------
    
    X = (X + alpha*sn);
end

    
