% clc;clear;
function K = element_bend(E,Iz,L)% E is elastic modulus 
                                 % Iz is inertia moment of z axis
                                 % L is the length of the beam
%%
%Suppose the shear displacement v(x)=a*x^3+b*x^2+c*x+d
syms x a b c d v1 v2 theta1 theta2;
A = [0 0 0 1;
     0 0 1 0;
     L^3 L^2 L 1;
     3*L^2 2*L 1 0];
B = [v1;theta1;v2;theta2];
Z = [a;b;c;d];
eq = A*Z==B;
[a b c d] = solve(eq,Z);
vx = [a b c d]*[x^3;x^2;x;1];
N = jacobian(vx,B);%N is shape function matrix
%%
%creat stiffness matrix
K = sym('k',[4 4]);
K(1,:) = E*Iz*diff(N,x,3);
K(3,:) = -K(1,:);
K(4,:) = E*Iz*diff(N,x,2);
K(2,:) = -K(4,:);
for i=1:4
    if i==1 | i==2
        K(i,:)=subs(K(i,:),x,0);
    else
        K(i,:)=subs(K(i,:),x,L);
    end
end
end



