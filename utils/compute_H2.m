function gamma2 = compute_H2(sysd, tol, opts)
% compute H2 given sysd
A = sysd.A;
B = sysd.B;
E = sysd.E;
C = sysd.C;
D = sysd.D1;

n = size(A,2);
m = size(B,2);
r = size(C,1);

Y = sdpvar(n);
S = sdpvar(m,n,'full');
Z = sdpvar(r);
gamma2 = sdpvar(1);

C1 = [Y-E*E' A*Y + B*S; Y*A' + S'*B' Y];
C2 = [Z C*Y+D*S; Y*C'+S'*D' Y];
C3 = gamma2-trace(Z);

cons = [Y >= tol; Z >= tol;...
    C1 >= tol; C2 >= tol; C3 >= tol];
optimize(cons, gamma2, opts);
gamma2 = sqrt(value(gamma2));
end