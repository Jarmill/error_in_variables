% Implementation of
% "Data-Driven Superstabilizing Control of Error-in-Variables Discrete-Time Linear Systems"

clear
yalmip('clear')
rng(1, 'twister')
tic
%% parameter
umax = 1;                   % input bound
eps = 0.01;                 % noise bound    dx
T = 5;                      % # of samples
d = 2;                      % degree of psatz
tol = 1e-6;                 % delta in paper
opts = sdpsettings('solver','mosek','verbose', 0);

%% generate system
A = 1.1;
B = 0.5;
n = size(A,1);              % dim of state
m = size(B,2);              % dim of input
sysd = generate_sys(n,m,A,B);

%% generate trajectory
U = (2*rand(m, T-1)-1)*umax;
X = zeros(n,T);
X(:, 1) = rand(n,1);
noise_x = (2*rand(size(X))-1)*eps;
for t = 1:(T-1)
    X(:,t+1) = A*X(:,t) + B*U(:,t);
end
X_noise = X + noise_x;
sim = struct('U',U,'X_noise',X_noise,'epsilon',eps,'tolerance',tol);

%% solve SS
type = 'no_prior';
obj = '[]';
out = Full_Pos(sim, d, T, sysd, type, obj);
sol = optimize(out.cons, out.obj, opts)

%% extract solution
v = value(out.v);
Y = value(out.Y);
K = Y * diag(1./v)
t = toc


