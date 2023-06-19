load('all_sim_test.mat')


clear
yalmip('clear')
rng(1, 'twister')
tic
%% parameter
umax = 1;                   % input bound
% eps = [0.01; 0.02; 0.01];   % noise bound    dx,du,w
eps = [0.1; 0.2; 0.1];   % noise bound    dx,du,w
T = 10;                      % # of samples
opts = sdpsettings('solver','mosek','verbose', 0);
tol = 1e-6;

%% generate system
% % 2nd order system, works for dual method
A_true = [0.6863    0.3968
    0.3456    1.0388];
B_true = [0.4170    0.0001
    0.7203    0.3023];
n = size(A_true,1);              % dim of state
m = size(B_true,2);              % dim of input
sysd = generate_sys(n,m,A_true,B_true);

%% generate trajectory
U = (2*rand(m, T-1)-1)*umax;
X = zeros(n,T);
X(:, 1) = rand(n,1);
noise_x = (2*rand(size(X))-1)*eps(1);
noise_u = (2*rand(size(U))-1)*eps(2);
noise_w = (2*rand(size(X))-1)*eps(3);
for t = 1:(T-1)
    X(:,t+1) = A_true*X(:,t) + B_true*U(:,t) + noise_w(:,t);
end
X_noise = X + noise_x;
U_noise = U + noise_u;
sim = struct('X_noise',X_noise,'U_noise',U_noise,'epsilon',eps,'tolerance',tol);


%% find a solution (assuming process noise alone)
A= sdpvar(n, n, 'full');
B= sdpvar(n, m, 'full');

cons = 0;

%fails
% eps_w = sim.epsilon(3); %fails
eps_w = 2.20*sim.epsilon(3); %fails
eps_w = 2.21*sim.epsilon(3); %succeeds

T = size(sim.U_noise, 2);
R = zeros(n, T, 'like', sdpvar)
for i = 1:T
    r = sim.X_noise(:, i+1) - A*sim.X_noise(:, i) - B*sim.U_noise(:, i);
    R(:, i) = r;
    cons = [cons; norm(r, 'inf') <= eps_w];
end

R_true = replace(R, [A, B], [A_true, B_true])

