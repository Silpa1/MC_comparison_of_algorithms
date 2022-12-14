function [x,flag,resNE,iter] = cgls_new(A,At,b,shift,tol,maxit,prnt,x0)

if (nargin < 4 || isempty(shift)), shift = 0;    end
if (nargin < 5 || isempty(tol))  , tol   = 1e-6; end
if (nargin < 6)                  , maxit = [];   end
if (nargin < 7 || isempty(prnt)) , prnt  = 0;    end
if (nargin < 8)                  , x0 = [];      end


if isa(A, 'numeric')
    explicitA = true;
elseif isa(A, 'function_handle')
    explicitA = false;
else
    error('A must be numeric or a function handle.');
end

if explicitA
    [m,n] = size(A);

    if ~isempty(x0)
        x = x0;
    else
        x = zeros(n,1);
    end
    
    r = b - A*x;
    s = A'*r-shift*x;
else
    m = size(b,1);
    
    if ~isempty(x0)
        x = x0;
        r = b - A(x,1);
        s = A(r,2) - shift*x;
        n = size(s,1);
    else
        r = b;
        %s = A(b,2);
        s=sum(At(b),3);
%         s=sum(,3);
        [n1,n2] = size(s);
        x = zeros(n1,n2);
    end
    
end


if isempty(maxit)
    maxit = min([m,n,20]);
end

p      = s;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
norms0 = norm(s,'fro');
gamma  = norms0^2;
normx  = norm(x,'fro');
xmax   = normx;
k      = 0;
flag   = 0;

if prnt
    head = '    k       x(1)             x(n)           normx        resNE';
    form = '%5.0f %16.10g %16.10g %9.2g %12.5g\n';
    disp('  ');   disp(head);
    fprintf(form, k, x(1), x(n), normx, 1);
end

indefinite = 0;

%--------------------------------------------------------------------------
% Main loop
%--------------------------------------------------------------------------
while (k < maxit) && (flag == 0)
    
    k = k+1;
    
    if explicitA
        q = A*p;
    else
        %%%%%%%%%%%%%%%%%%%
        q = A(p);
    end
        
    delta = norm(q(:))^2  +  shift*norm(p(:))^2;
    if delta <= 0, indefinite = 1;   end
    if delta == 0, delta      = eps; end
    alpha = gamma / delta;
    
    x     = x + alpha*p;
    r     = r - alpha*q;
       
    if explicitA
        s = A'*r - shift*x;
    else
        %%%%%%%%%%%%%%%%%%
        s = sum(At(r),3) - shift*x;
    end
   
    norms  = norm(s(:));
    gamma1 = gamma;
    gamma  = norms^2;
    beta   = gamma / gamma1;
    p      = s + beta*p;
    
    % Convergence
    normx = norm(x(:));
    xmax  = max(xmax, normx);
    flag  = (norms <= norms0 * tol) || (normx * tol >= 1);
    
    % Output
    resNE = norms / norms0;
    if prnt, fprintf(form, k, x(1), x(n), normx, resNE); end
end % while

iter = k;

shrink = normx/xmax;
if k == maxit,          flag = 2; end
if indefinite,          flag = 3; end
if shrink <= sqrt(tol), flag = 4; end