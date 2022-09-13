function out = Full_SS(sim, d, T, type, obj)
%  sim:     sampled trajectories
%    d:     degree of the psatz
%    T:     # of samples for design
% type:     'no_prior', system unknown
%           'prior',    system partically known
%  obj:     '[]',       stabilization
%           'lambda',   minimize lambda

% extract data and size
X = sim.X_noise(:,1:T);     % state data
U = sim.U(:,1:T-1);         % input data
eps = sim.epsilon;       % noise bound epsilon
tol = sim.tolerance;        % delta in paper
n = size(X,1);              % dim of state
m = size(U,1);              % dim of input
Cons = [];  

switch type
    case 'no_prior'   % A,B,Dx,K are all variables
        A = sdpvar(n,n,'full');
        B = sdpvar(n,m,'full');
        K = sdpvar(m,n,'full');
        Dx = sdpvar(n,T);           % delta x in paper
        vars = [A(:);B(:);Dx(:)];  
        Acl = A+B*K;
    case 'prior'
        % some entries are known
        A_var = sdpvar(n,1,'full');
        K = sdpvar(m,n,'full');
        Dx = sdpvar(n,T);
        vars = [A_var(:);Dx(:)];
        A = [A_var(1) 0.3968; A_var(2) 1.0388];
        B = [0.4192; 0.6852];
        Acl = A+B*K;
end

switch obj
    case '[]'
        lambda = 1;
        obj = [];
    case 'lambda'
        lambda = sdpvar(1);
        obj = lambda;
end

% constraints
g = [eps-Dx(:);eps+Dx(:)];        % eq.(9.2) inequality constraint
for i = 1:T-1                     % eq.(9.1) equality constraint
    ht(1+n*(i-1):n*i,:) = X(:,i+1) - A*X(:,i) - B*U(:,i);              
    h(1+n*(i-1):n*i,:) = -Dx(:,i+1) + A*Dx(:,i) + ht(1+n*(i-1):n*i,:);  
end
cons_data = struct('ineq', g, 'eq', h);    % constraints from data

% define M
M = zeros(n, n, 'like', sdpvar);
for i = 1:n
    for j = 1:n
        M(i, j) = polynomial(vars, 2*d);
    end
end
                  
p1 = M-Acl;
p2 = M+Acl;
p3 = lambda-sum(M,2)-tol;
P = [p1(:);p2(:);p3(:)];    % P contains M-A-BK >= 0, M+A+BK >= 0, 1-delta-sum(M,2) >= 0

% get constraints
Gram = cell(2*n^2+n,1);
Coef = cell(2*n^2+n,1);
for i = 1:2*n^2+n
    [p_psatz, cons_psatz, Gram{i,1}, Coef{i,1}] = Full_SS_psatz(P(i), cons_data, d, vars);
    Cons = [Cons; cons_psatz];      % constraints from psatz
end
out.poly = p_psatz;
out.cons = Cons;
out.Gram = Gram;
out.Coef = Coef;
out.K = K;
out.obj = obj;
out.vars = vars;
end