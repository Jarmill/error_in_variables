function out = process_bound_EIV(sim, d, T, type)
%estimate the minimum process-bound radius needed to contain the EIV
%feasible set
%
%  sim:     sampled trajectories
%    d:     degree of the psatz
%    T:     # of samples for design
% type:     'no_prior', system unknown
%           'prior',    system partically known

% extract data and size
X = sim.X_noise(:,1:T);     % state data
U = sim.U_noise(:,1:T-1);   % input data
eps = sim.epsilon;          % noise bound epsilon
tol = sim.tolerance;        % delta in paper
n = size(X,1);              % dim of state
m = size(U,1);              % dim of input

switch type
    case 'no_prior'
        A = sdpvar(n,n,'full');
        B = sdpvar(n,m,'full');
%         K = sdpvar(m,n,'full');
        vars = [A(:);B(:)];
%         Acl = A+B*K;
    case 'prior'     % second column of A known
        A_var = sdpvar(1,n,'full');
        B_var = sdpvar(n,m,'full');
%         K = sdpvar(m,n,'full');
        vars = [A_var(:);B_var(:)];
        A = [A_var(1) 0.3968; A_var(2) 1.0388];
        B = B_var;
%         Acl = A+B*K;
end
eps_wrap = sdpvar(1, 1);


% define constraints
g = [];
if eps(1) ~= 0
    g = [g; eps(1)*ones(2*n*T,1)];
end
if eps(2) ~= 0
    g = [g; eps(2)*ones(2*m*(T-1),1)];
end
if eps(3) ~= 0
    g = [g; eps(3)*ones(2*n*(T-1),1)];
end
for i = 1:T-1
    h(1+n*(i-1):n*i,:) = X(:,i+1) - A*X(:,i) - B*U(:,i);
end
cons_data = struct('ineq', g, 'eq', h);  

% define M
M = zeros(n, n, 'like', sdpvar);
for i = 1:n
    for j = 1:n
        M(i, j) = polynomial(vars, 2*d);
    end
end

% define non-negative polynomial

% P = [p1(:);p2(:);p3(:)];    % P contains M-A-BK >= 0, M+A+BK >= 0, 1-delta-sum(M,2) >= 0
P = [eps_wrap - h; eps_wrap+h];
% get constraints
Cons = [];
Gram = cell(size(P));
Coef = cell(size(P));
for i = 1:length(P)
    [p_psatz, cons_psatz, Gram{i,1}, Coef{i,1}] = Dual_SS_psatz_all_noise(P(i), cons_data, d, vars, A, B, eps, n, m ,T);
    Cons = [Cons; cons_psatz];
end

out.poly = p_psatz;
out.cons = Cons;
out.Gram = Gram;
out.Mu = Coef;
out.eps_wrap = eps_wrap;
out.obj = eps_wrap;
out.vars = vars;
end