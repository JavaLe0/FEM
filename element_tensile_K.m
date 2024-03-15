function K = element_tensile(E,A,L)
M = [1,-1;-1,1];
K = E*A/L*M;
end