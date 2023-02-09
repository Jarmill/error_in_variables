function out = Full_Pos(sim, d, T, type, obj)
%DUAL_POS positive systems stabilization
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
eps = sim.epsilon;          % noise bound epsilon
tol = sim.tolerance;        % delta in paper
n = size(X,1);              % dim of state
m = size(U,1);              % dim of input


Dx = sdpvar(n,T);           % delta x in paper
switch type
    case 'no_prior'
        A = sdpvar(n,n,'full');
        B = sdpvar(n,m,'full');        
        vars = [A(:);B(:);Dx(:)];  

    case 'prior'     % second column of A known
        A_var = sdpvar(1,n,'full');
        B_var = sdpvar(n,m,'full');
        vars = [A_var(:);Dx(:)];
        A = [A_var(1) 0.3968; A_var(2) 1.0388];
        B = B_var;        
end


%design parameter variable definitions
v = sdpvar(n, 1);    
Y = sdpvar(m,n,'full');
V = diag(v);
Acl = A*V+B*Y;

switch obj
    %ignore objective for now. add it later
    case '[]'
        lambda = 1;
        Cons = [];
    case 'lambda'
        lambda = sdpvar(1);
        Cons = (lambda <= 1-tol);
    otherwise
        disp('You need to define objective')
end
Cons = [Cons; sum(v)==1; v>=tol];

% constraints
% constraints
g = [eps-Dx(:);eps+Dx(:)];        % eq.(9.2) inequality constraint
for i = 1:T-1                     % eq.(9.1) equality constraint
    ht(1+n*(i-1):n*i,:) = X(:,i+1) - A*X(:,i) - B*U(:,i);              
    h(1+n*(i-1):n*i,:) = -Dx(:,i+1) + A*Dx(:,i) + ht(1+n*(i-1):n*i,:);  
end
cons_data = struct('ineq', g, 'eq', h);    % constraints from data



% define non-negative polynomials
p1 = Acl;
p2 = v -sum(Acl,2)  -tol;
P = [p1(:);p2(:)];    % P contains M-A-BK >= 0, M+A+BK >= 0, 1-delta-sum(M,2) >= 0

% get psatz constraints
Gram = cell(n^2+n,1);
Coef = cell(n^2+n,1);
for i = 1:n^2+n
    [p_psatz, cons_psatz, Gram{i,1}, Coef{i,1}] = Full_SS_psatz(P(i), cons_data, d, vars);
    Cons = [Cons; cons_psatz];      % constraints from psatz
end
out.poly = p_psatz;
out.cons = Cons;
out.Gram = Gram;
out.Coef = Coef;
out.v = v;
out.Y = Y;
out.obj = lambda;
out.vars = vars;
end