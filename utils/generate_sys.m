function sysd = generate_sys(n,m,A,B)
r = n+m;
e = n;
if nargin < 3
    sysd.A = rand(n);
    sysd.B = rand(n,m);
else
    sysd.A = A;
    sysd.B = B;
end
sysd.E = eye(n,e);
sysd.C = eye(r,n);
sysd.D1 = [zeros(n,m); eye(m)];
sysd.D2 = eye(r,e);
sysd.F = 0;
sysd.r = r;
sysd.e = e;
end