function out = Dual_SS_manual(sim, d, T, obj)
% this is an alternative and efficienet way implementing the algorithm 3
% instead of defining the polynomial and multiply them,
% we directly define the multiplication of the coefficient
% current version does not support prior info

% extract data and size
sim.X_noise = sim.X_noise(:,1:T);   % state data
sim.U = sim.U(:,1:T-1);             % input data
tol = sim.tolerance;                % delta in paper
n = size(sim.X_noise,1);            % dim of state
m = size(sim.U,1);                  % dim of input

sM = nchoosek(n^2+n*m+2*d,2*d);     % size of coefficient vector of each entry in M(A,B)
cM = sdpvar(sM,n^2,'full');         % coefficient matrices of M(A,B)
K = sdpvar(m,n,'full');

Acl = zeros(sM,n^2,'like',sdpvar);         % coefficient of A+BK with variable A,B
Acl(2+n^2:1+n^2+n*m,:) = kron(K,eye(n));   % BK
Acl(2:n^2+1,:) = eye(n^2);                 % A

switch obj
    case '[]'
        lambda = 1;
        Cons = [];
    case 'lambda'
        lambda = sdpvar(1);
        Cons = lambda <= 1-tol;
    otherwise
        disp('You need to define obj')
end

% define coefficients of non-negative polynomials
cp1 = cM-Acl;   % coefficient of M-A-BK
cp2 = cM+Acl;   % coefficient of M+A+BK
for i = 1:n   
    cp3(:,i) = lambda*eye(sM,1) - sum(cM(:,i:n:end),2) - tol*eye(sM,1);   % 1-delta-sum(M,2) >= 0
end
cP = [cp1 cp2 cp3];      % cP contains coefficients of M-A-BK >= 0, M+A+BK >= 0, 1-delta-sum(M,2) >= 0

% get psatz constraints
for i = 1:2*n^2+n       
    cons = Dual_SS_psatz_manual(sim, cP(:,i), d);
    Cons = [Cons; cons];
end
out.cons = Cons;
out.K = K;
out.cM = cM;
out.obj = lambda;
end






