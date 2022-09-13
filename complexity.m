%% used to check size of psatz, follows from table 1,2
clc
n = 3;
m = 2;
T = 6;
r = n+m;
e = n;
%% SS   1 for full method, 2 for alternative method
p1 = n*(n+m+T);     % number of variables
p2 = n*(n+m);

d1 = 2;             % degree of full method
d2 = 1;             % degree of alternative method
% full size                             % P = sigma0 + \sum sigmai*g + \sum muj*h
q1 = nchoosek(p1 + 2*d1, 2*d1)          % size of P
s01 = nchoosek(p1 + d1, d1)             % size of sigma0
s11 = nchoosek(p1 + d1-1,d1-1)          % size of sigmai
mu1 = nchoosek(p1 + 2*d1-2, 2*d1-2)     % size of muj
% alternative size
q2 = nchoosek(p2 + 2*d2, 2*d2)
s02 = nchoosek(p2 + d2, d2)
s12 = nchoosek(p2 + d2, d2)
mu2 = nchoosek(p2 + 2*d2-1, 2*d2-1)
% total number of variables
total1 = (q1 + s01*(s01+1)/2 + s11*(s11+1)*n*T + mu1*2*n*(T-1))*(n^2+n)
total2 = (q2 + s02*(s02+1)/2 + s12*(s12+1)*n*T + mu2*2*n*(T-1))*(n^2+n)

%% QS
l = 2*n;            % size of psd P
nu = (l+1)*l/2;     % number of variables in P

q1_QS = nu*nchoosek(p1 + 2*d1, 2*d1)      
s01_QS = l*nchoosek(p1 + d1, d1)
s11_QS = l*nchoosek(p1 + d1-1,d1-1)
mu1_QS = nu*nchoosek(p1 + 2*d1-2, 2*d1-2)

q2_QS = nu*nchoosek(p2 + 2*d2, 2*d2)
s02_QS = l*nchoosek(p2 + d2, d2)
s12_QS = l*nchoosek(p2 + d2, d2)
mu2_QS = nu*nchoosek(p2 + 2*d2-1, 2*d2-1)

total1_QS = nu*q1_QS + s01_QS*(s01_QS+1)/2 + s11_QS*(s11_QS+1)*n*T + mu1_QS*2*n*(T-1)
total2_QS = nu*q2_QS + s02_QS*(s02_QS+1)/2 + s12_QS*(s12_QS+1)*n*T + mu2_QS*2*n*(T-1)
%% Hinf
l = 2*n+e+r;
nu = (l+1)*l/2;

q1_inf = nu*nchoosek(p1 + 2*d1, 2*d1)      
s01_inf = l*nchoosek(p1 + d1, d1)
s11_inf = l*nchoosek(p1 + d1-1,d1-1)
mu1_inf = nu*nchoosek(p1 + 2*d1-2, 2*d1-2)

q2_inf = nu*nchoosek(p2 + 2*d2, 2*d2)
s02_inf = l*nchoosek(p2 + d2, d2)
s12_inf = l*nchoosek(p2 + d2, d2)
mu2_inf = nu*nchoosek(p2 + 2*d2-1, 2*d2-1)

total1_inf = nu*q1_inf + s01_inf*(s01_inf+1)/2 + s11_inf*(s11_inf+1)*n*T + mu1_inf*2*n*(T-1)
total2_inf = nu*q2_inf + s02_inf*(s02_inf+1)/2 + s12_inf*(s12_inf+1)*n*T + mu2_inf*2*n*(T-1)