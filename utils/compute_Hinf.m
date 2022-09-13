function gamma = compute_Hinf(sysd,tol,opts)
% compute Hinf given sysd
A = sysd.A;
B = sysd.B;
E = sysd.E;
C = sysd.C;
D1 = sysd.D1;
D2 = sysd.D2;

n = size(A,2);
m = size(B,2);
r = size(C,1);
e = size(E,2);

Y = sdpvar(n);
S = sdpvar(m,n,'full');
gamma = sdpvar(1);

F = [Y A*Y-B*S E zeros(n,r);
    (A*Y-B*S)' Y zeros(n,e) (C*Y-D1*S)';
    E' zeros(e,n) gamma*eye(e) D2';
    zeros(r,n) C*Y-D1*S D2 gamma*eye(r)];
cons = [Y >= tol; F >= tol];
optimize(cons, gamma, opts);
gamma = value(gamma);
end