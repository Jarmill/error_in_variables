function out = Dual_SS_Extended_all_noise(sim, d, T, type)
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

v = sdpvar(n, 1);
Y = sdpvar(m,n,'full');
Acl = A*diag(v)+B*Y;


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
p1 = M-Acl;
p2 = M+Acl;
p3 = v-sum(M,2)-tol;
P = [p1(:);p2(:);p3(:)];    % P contains M-A-BK >= 0, M+A+BK >= 0, 1-delta-sum(M,2) >= 0

% get constraints
Cons = [v>= tol];
Gram = cell(2*n^2+n,1);
Coef = cell(2*n^2+n,1);
for i = 1:2*n^2+n
    [p_psatz, cons_psatz, Gram{i,1}, Coef{i,1}] = Dual_SS_psatz_all_noise(P(i), cons_data, d, vars, A, B, eps, n, m ,T);
    Cons = [Cons; cons_psatz];
end

out.poly = p_psatz;
out.cons = Cons;
out.Gram = Gram;
out.Mu = Coef;
out.v = v;
out.Y = Y;
out.obj = 0;
out.vars = vars;
end