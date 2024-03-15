%Suppose the torsional angle φ(x)=a+bx
function K = element_torsion(G,Ir,L)
syms f1 f2 a b x
M = [1,0;1,L];
eq = M*[a;b]==[f1;f2];
[a b] = solve(eq,a,b);
fx = a+b*x;
N = jacobian(fx,[f1,f2]);
K(1,:) = -G*Ir*subs(diff(N,x),x,0);
K(2,:) = G*Ir*subs(diff(N,x),x,L);
end