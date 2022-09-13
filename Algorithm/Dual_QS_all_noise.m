function out = Dual_QS_all_noise(sim, d, T, sysd, type, obj)
%  sim:     sampled trajectories
%    d:     degree of the psatz
%    T:     # of samples for design
% sysd:     contains C,D,E
% type:     'no_prior', system unknown
%           'prior',    system partically known
%  obj:     '[]',       stabilization
%           'H2',       minimize H2
%           'Hinf',     minimize Hinf

% extract data and size
X = sim.X_noise(:,1:T);     % state data
U = sim.U_noise(:,1:T-1);         % input data
eps = sim.epsilon;       % noise bound epsilon
tol = sim.tolerance;        % delta in paper
n = size(X,1);              % dim of state
m = size(U,1);              % dim of input

% extract performance matrix C,D,E
r = sysd.r;
e = sysd.e;
E = sysd.E;         
C = sysd.C;
D1 = sysd.D1;       % D1 used in H2
D2 = sysd.D2;       % D2 used in Hinf

switch type
    case 'no_prior'
        A = sdpvar(n,n,'full');
        B = sdpvar(n,m,'full');
        vars = [A(:);B(:)];
    case 'prior'        % first row of A known
        A_var = sdpvar(1,n,'full');
        B_var = sdpvar(n,m,'full');
        vars = [A_var(:);B_var(:)];
        A = [0.6863 0.3968; A_var(1) A_var(2)];
        B = B_var;
end

% constraints
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

% define sdp F
Y = sdpvar(n);
S = sdpvar(m,n,'full');
switch obj
    case '[]'           % Algotithm 5
        F = [Y A*Y+B*S; Y*A'+S'*B' Y]-tol*eye(2*n);
        Cons = (Y >= tol);
        gamma = 0;
    case 'H2'           % Algorithm 6
        Z = sdpvar(r);
        eta = sdpvar(1);
        gamma = sdpvar(1);
        C1 = [Y-E*E' A*Y+B*S; Y*A'+S'*B' Y];
        C2 = [Z C*Y+D1*S; Y*C'+S'*D1' Y];
        C3 = gamma-trace(Z);
        F = C1-eta*eye(size(C1));
        Cons = [Y >= tol; Z >= tol; eta >= tol; C2 >= tol; C3 >= tol];
    case 'Hinf'         % Algorithm 7   does not work, bug to be fixed
        eta = sdpvar(1);
        gamma = sdpvar(1);
        F = [Y A*Y-B*S E zeros(n,r);
            (A*Y-B*S)' Y zeros(n,e) (C*Y-D1*S)';
            E' zeros(e,n) gamma*eye(e) D2';
            zeros(r,n) C*Y-D1*S D2 gamma*eye(r)];
        Cons = [Y >= tol; eta >= tol];
        F = F - eta*eye(size(F));
end
[p_psatz, cons_psatz, Gram, Coef] = Dual_QS_psatz_all_noise(F, cons_data, d, vars, A, B, eps, n, m, T);
Cons = [Cons; cons_psatz];

out.poly = p_psatz;
out.cons = Cons;
out.Gram = Gram;
out.Coef = Coef;
out.Y = Y;
out.S = S;
out.obj = gamma;
out.vars = vars;
end