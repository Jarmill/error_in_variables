function gamma2 = compute_H2_K(sysd, K)
% compute closed loop H2 given K
E = sysd.E;
A_K = sysd.A+sysd.B*K;
C_K = sysd.C+sysd.D1*K;
P = dlyap(A_K',C_K'*C_K);
gamma2 = trace(E'*P*E);   %
gamma2 = sqrt(gamma2);
end