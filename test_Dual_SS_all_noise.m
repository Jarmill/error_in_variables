% Implementation of
% "Data-Driven Stabilizing and Robust Control of Discrete-Time Linear Systems with Error in Variables"

clear
yalmip('clear')
rng(1, 'twister')
tic
%% parameter
umax = 1;                   % input bound
% eps = [0.01; 0.02; 0.01];   % noise bound    dx,du,w
eps = [0.1; 0.2; 0.1];   % noise bound    dx,du,w
T = 5;                      % # of samples
d = 1;                      % degree of psatz
tol = 1e-6;                 % delta in paper
opts = sdpsettings('solver','mosek','verbose', 0);

%% generate system
% % 1st order system, used to compare with full method
% A = 1.1;
% B = 0.5;
% % 2nd order system, works for dual method
A = [0.6863    0.3968
    0.3456    1.0388];
B = [0.4170    0.0001
    0.7203    0.3023];
n = size(A,1);              % dim of state
m = size(B,2);              % dim of input
sysd = generate_sys(n,m,A,B);

%% generate trajectory
U = (2*rand(m, T-1)-1)*umax;
X = zeros(n,T);
X(:, 1) = rand(n,1);
noise_x = (2*rand(size(X))-1)*eps(1);
noise_u = (2*rand(size(U))-1)*eps(2);
noise_w = (2*rand(size(X))-1)*eps(3);
for t = 1:(T-1)
    X(:,t+1) = A*X(:,t) + B*U(:,t) + noise_w(:,t);
end
X_noise = X + noise_x;
U_noise = U + noise_u;
sim = struct('X_noise',X_noise,'U_noise',U_noise,'epsilon',eps,'tolerance',tol);

%% solve SS
type = 'no_prior';
obj = 'lambda';
out = Dual_SS_all_noise(sim, d, T, type, obj);
% out = Dual_SS_all_noise_manual(sim, d, T, obj);  % does not support prior
sol = optimize(out.cons, out.obj, opts)

%% extract solution
K = value(out.K);
lambda = value(out.obj)
t = toc


