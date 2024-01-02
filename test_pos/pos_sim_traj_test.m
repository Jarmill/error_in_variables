clear
yalmip('clear')
rng(1, 'twister')
tic
%% parameter
umax = 1;                   % input bound
% eps = 0.1;                 % noise bound    dx
% T = 6;                      % # of samples
% eps = 0.14;                 % noise bound    dx

eps = [0.08, 0.04];
T = 8;                      % # of samples

% epslist = [0.05; 0.08; 0.11; 0.14];
epslist = [0.05; 0.08; 0.11; 0.14]*[1, 1];
Tlist = [8; 10; 12; 14];
% T = 8;
d = 1;                      % degree of psatz
tol = 1e-5;                 % delta in paper
opts = sdpsettings('solver','mosek','verbose', 1);

%% generate system
% % 1st order system, used to compare with full method
% A = 1.1;
% B = 0.5;
% % 2nd order system, works for dual method
%
% stabilizable but not p2p
% A = [0.6863    0.3968
%     0.3456    1.0388];
% B = [0.4170    0.0001
%     0.7203    0.3023];

%works with p2p
A = [0.302566731143229,1.013365160613880;0.642592615777238,0.341971254054088];
B = [-0.230751776840409,0.0677675585788196;-0.459150458710386,0.871935186168177];
n = size(A,1);              % dim of state
m = size(B,2);              % dim of input
sysd = generate_sys(n,m,A,B);

%% generate trajectories
Ntraj = 50;
sim = pos_sim_eiv(Ntraj, A, B, umax, eps, T, tol);

sim_time = cell(length(Tlist), 1);
sim_eps = cell(length(epslist), 1);

for i = 1:length(Tlist)
    sim_time{i} = pos_sim_eiv(Ntraj, A, B, umax, eps, Tlist(i), tol);
end
for i = 1:length(epslist)
    sim_eps{i} = pos_sim_eiv(Ntraj, A, B, umax, epslist(i, :), T, tol);
end

save('pos_traj_std_u.mat', 'sim', 'sysd', 'sim_time', 'sim_eps');

function sim_out = pos_sim_eiv(Ntraj, A, B, umax, eps, T, tol)
%generate samples of the trajectory starting from A, B
n = size(A,1);              % dim of state
m = size(B,2);              % dim of input

if length(eps)<2
    epsx = eps;
    epsu = 0;
else
    epsx = eps(1);
    epsu = eps(2);
end

sim_out = cell(Ntraj, 1);
for i = 1:Ntraj
    U = (2*rand(m, T-1)-1)*umax;
    X = zeros(n,T);
    X(:, 1) = rand(n,1);
    noise_x = (2*rand(size(X))-1)*epsx;
    noise_u = (2*rand(size(U))-1)*epsu;
    for t = 1:(T-1)
        X(:,t+1) = A*X(:,t) + B*U(:,t);
    end
    X_noise = X + noise_x;
    U_noise = U + noise_u;
    sim_out{i} = struct('U',U,'U_noise', U_noise, 'X_noise',X_noise,'epsilon',eps,'tolerance',tol, 'T', T);
end
end