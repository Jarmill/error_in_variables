% Implementation of
% "Data-Driven Superstabilizing Control of Error-in-Variables Discrete-Time Linear Systems"

clear
yalmip('clear')
rng(1, 'twister')
tic

load('data_eps0.05-0.14.mat')
sim = data{4, 1};
d = 1;
% T = size(sim.X, 2);
A = [0.6863    0.3968
    0.3456    1.0388];
B = [0.4170    0.0001
    0.7203    0.3023];

T = 14;
simc = sim;
simc.X = sim.X(:, 1:T);
simc.U = sim.U(:, 1:T);
simc.X_noise = sim.X_noise(:, 1:T);
%% solve SS
type = 'no_prior';
obj = '[]';
out = Dual_SS_Extended(sim, d, T, type, obj);
% out = Dual_SS_manual(sim, d, T, obj);   % does not support prior
opts = sdpsettings('solver','mosek','verbose', 1);
sol = optimize(out.cons, out.obj, opts)

%% extract solution
if sol.problem==0
v = value(out.v);
Y = value(out.Y);
K = Y * diag(1./v)
% lambda = value(out.obj)
t = toc
else
    disp('infeasible');
end


