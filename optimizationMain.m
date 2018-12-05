%Example 9.1/9.2
tol = 1.0e-8;
nmax = 50;

syms x1 x2
fprintf('Example 9.1/9.2 \n')
f1 = x1^2-2*x1-x2+1;
f2 = x1^2+x2^2-1;
g = [f1;f2]
jac = jacobian([f1,f2], [x1,x2]);
x = [1;1]
[xtst1,ktst1] = newtons(g,jac,x,tol,nmax)

fprintf('((x1 - 3)^2 + (x2 - 2)^2)^4 \n')
fn = ((x1 - 3)^2 + (x2 - 2)^2)^4
x = [1;1]
[xtst2,ktst2] = newtonsM(fn,x,tol,nmax)
