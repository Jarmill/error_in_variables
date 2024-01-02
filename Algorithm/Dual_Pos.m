function out = Dual_Pos(sim, d, T, sysd, type, obj)
%DUAL_POS positive systems stabilization
%  sim:     sampled trajectories
%    d:     degree of the psatz
%    T:     # of samples for design
% sysd:     system description (performance)
% type:     'no_prior', system unknown
%           'prior',    system partically known
%  obj:     '[]',       stabilization
%           'p2p',      minimize peak-to-peak gain

% extract data and size
X = sim.X_noise(:,1:T);     % state data

if isfield(sim, 'U_noise')
    U = sim.U_noise(:, 1:T-1);
else
    U = sim.U(:,1:T-1);         % input data
end
eps = sim.epsilon;          % noise bound epsilon

if length(sim.epsilon)==2
    epsu = sim.epsilon(2);
else
    
    epsu = 0;   
end
eps = eps(1);

tol = sim.tolerance;        % delta in paper
n = size(X,1);              % dim of state
m = size(U,1);              % dim of input



switch type
    case 'no_prior'
        A = sdpvar(n,n,'full');
        B = sdpvar(n,m,'full');        
        vars = [A(:);B(:)];

    case 'prior'     % second column of A known
        A_var = sdpvar(1,n,'full');
        B_var = sdpvar(n,m,'full');
        vars = [A_var(:);B_var(:)];
        A = [A_var(1) 0.3968; A_var(2) 1.0388];
        B = B_var;        
end

%design parameter variable definitions
v = sdpvar(n, 1);    
Y = sdpvar(m,n,'full');
V = diag(v);
Acl = A*V+B*Y;

Cons = [v >= tol];
switch obj
    %ignore objective for now. add it later
    case '[]'
        lambda = 1;
        Cons = [Cons; sum(v)==1];
        E = 0;        
    case 'p2p'
        lambda = sdpvar(1);
        Ccl = sysd.C*V + sysd.D1*Y;
        E = sysd.E;
        Cons = [Cons; Ccl(:) >= 0; ...
            lambda - sum(Ccl, 2) - sum(sysd.F, 2) >= tol];
        
%     case 'lambda'
%         lambda = sdpvar(1);
%         Cons = (lambda <= 1-tol);
    otherwise
        disp('You need to define objective')
end


% constraints
g = eps*ones(2*n*T,1);            % eq.(16) inequality constraint
if epsu
    gu = epsu*ones(2*m*(T-1),1);            % eq.(16) inequality constraint
else
    gu = [];
end

for i = 1:T-1                     % eq.(16) equality constraint
    h(1+n*(i-1):n*i,:) = X(:,i+1) - A*X(:,i) - B*U(:,i);
end
cons_data = struct('ineq', g, 'inequ', gu, 'eq', h);  



% define non-negative polynomials
p1 = Acl;
p2 = v -sum(Acl,2) - sum(E, 2)  -tol;
P = [p1(:);p2(:)];    % P contains M-A-BK >= 0, M+A+BK >= 0, 1-delta-sum(M,2) >= 0

% get psatz constraints
Gram = cell(n^2+n,1);
Coef = cell(n^2+n,1);
for i = 1:n^2+n
    [p_psatz, cons_psatz, Gram{i,1}, Coef{i,1}] = Dual_SS_psatz(P(i), cons_data, d, vars, A, B);
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