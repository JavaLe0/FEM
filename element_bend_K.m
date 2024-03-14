% clc;clear;
function K = element_bend(e,I,l)
%%
%The shear displacement v(x)=a*x^3+b*x^2+c*x+d
syms x a b c d L v1 v2 theta1 theta2;
A = [0 0 0 1;
     L^3 L^2 L 1;
     0 0 1 0;
     3*L^2 2*L 1 0];
B = [v1;v2;theta1;theta2];
Z = [a;b;c;d];
eq = A*Z==B;
[a b c d] = solve(eq,Z);
vx = [a b c d]*[x^3;x^2;x;1];
N = jacobian(vx,B);%N is shape function matrix
%%
%creat stiffness matrix
syms E Iz;% E is elastic modulus,Iz is inertia moment of z axis
K = sym('k',[4 4]);
K(1,:) = E*Iz*diff(N,x,2);
K(2,:) = K(1,:);
K(3,:) = E*Iz*diff(N,x,3);
K(4,:) = K(3,:);
for i=1:4
    if i==1 | i==3
        K(i,:)=subs(K(i,:),x,0);
    else
        K(i,:)=subs(K(i,:),x,l);
    end
end
K = subs(K,[E Iz L],[e I l]);
end



