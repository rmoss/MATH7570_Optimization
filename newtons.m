function [x,k] = newtons(fgrad,jac,x,tol,nmax)
%
% function ]x,k] = newtons(f,x,tol,nmax)
%
% This function returns in x a column vector x_k such that 
%   || x_k - x_{k-1} || < tol (1 + ||x_k||)
% and in k the number of iterations (Jacobian evaluations) required.
% On entry, x contains an innitial guess.
% If k equals nmax then no convergence has been reached.
%
% The iterates ||f(x_k)|| are recorded.  This option 
% can be easily turned off.

    %Initialize
    x = x(:); %ensure x is a column vector
    fprintf ('k      ||f(x_k)||         x1               x2                 error \n')
    format long g
    e = [];
    c = [];
    %Newton
    %syms x1 x2;
    %[f_grad,jac] = gradjac(f);
    
    for k=1:nmax
        [fx,Jx] = evaljacfunc(fgrad,jac,x);
        p = -Jx \ fx;
        x_next = x + p;
        fprintf ('%d    %e      %e      %e      %e      \n',k-1,norm(fx),x(1),x(2),norm(x_next-x))
        e = [e norm(x_next-x)];
        x = x_next;
        if norm(p) < tol*(1+norm(x))
            [fx,Jx] = evaljacfunc(fgrad,jac,x);
            fprintf ('%d    %e      \n',k,norm(fx) )
            for i=1:k-2
                c(i) = (log(e(i+2)/e(i+1)))/(log(e(i+1)/e(i)));
            end
    
%             c(k-1) = nan;
%             c(k) = nan;
            T = table(c.');
            T
            return
        end
    end
    
    for i=1:k-2
        c(i) = (log(e(i+2)/e(i+1)))/(log(e(i+1)/e(i)));
    end
    
%     c(k-1) = nan;
%     c(k) = nan;
    T = table(c.');
    T
    k = nmax;
end