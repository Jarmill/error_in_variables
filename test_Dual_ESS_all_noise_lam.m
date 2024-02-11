% Implementation of
% "Data-Driven Stabilizing and Robust Control of Discrete-Time Linear Systems with Error in Variables"

clear
yalmip('clear')
rng(1, 'twister')
tic
%% parameter
umax = 1;                   % input bound
% eps = [0.01; 0.02; 0];   % noise bound    dx,du,w
% eps = [0.01; 0.0; 0];   % noise bound    dx,du,w
eps = [0.04; 0.04; 0];   % noise bound    dx,du,w
T = 12;                      % # of samples
% T = 6;
d = 1;                      % degree of psatz
tol = 1e-6;                 % delta in paper
opts = sdpsettings('solver','mosek','verbose', 1);

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
% lam = 0.5;
% 
lam_use = 0.9;
% out = Dual_SS_all_noise_manual(sim, d, T, obj);  % does not support prior


% USE_OPTIMIZER = true;
% if USE_OPTIMIZER
    lam = sdpvar(1, 1);
    out = Dual_SS_Extended_all_noise_lam(sim, d, T, lam, type);

    P = optimizer(out.cons, out.obj, opts, lam, {out.v, out.Y, out.beta});

disp('Compiled model\n')
%% iterate over the free parameter
Nlam = 100;
% Nlam = 12;
% Nlam = 5;
lam_list = linspace(1, 0, Nlam+1);
lam_list = lam_list(2:end);


v_list = cell(Nlam, 1);
K_list = cell(Nlam, 1);
p2p_list = zeros(Nlam,1);
abs_eig_list = zeros(Nlam, n);
err_list = zeros(Nlam, 1);

for i = 1:Nlam
    [solc, err_list(i)] = P(lam_list(i));
    v_list{i} = value(solc{1});
    Y = value(solc{2});
    K_list{i} = Y * diag(1./v_list{i});
    beta = solc{3};
    p2p_list(i) = value(beta)/(1-lam_list(i));
    Acl = A + B*K_list{i};
    abs_eig_list(i, :) = abs(eig(Acl));
    disp(i)
end

save('p2p_ess_result.mat', 'eps', 'T', 'sim', 'v_list', 'K_list', 'p2p_list', 'abs_eig_list', 'd', 'err_list')
    


%% test on superstability alone

out_ss = Dual_SS_all_noise(sim, d, T, type, obj)
% out = Dual_SS_manual(sim, d, T, obj);   % does not support prior
sol_ss = optimize(out_ss.cons, out_ss.obj, opts)

%% extract solution
K_ss = value(out_ss.K);
lambda_ss = value(out_ss.obj)

[mm, ii] = min(p2p_list);
%% plot the figure
figure(2)
clf
hold on
scatter(lam_list(err_list==0), p2p_list(err_list==0), 'filled')
scatter(lam_list(ii), mm, 200, 'k')
ylim([4, max(p2p_list(err_list==0))])
set(gca, 'Yscale', 'log')
xlabel('$\lambda$', 'interpreter', 'latex', 'FontSize', 14)
ylabel('$\ell_\infty$ upper-bound', 'interpreter', 'latex', 'FontSize', 14)
title('Worst-Case Peak-to-Peak Regulation', 'Interpreter', 'latex', 'FontSize', 16)
legend({'d=1 bounds', sprintf('best bound %0.4f', mm)}, 'location', 'southeast')
    % solc = P(lam_use);
    % obj_rec = value(out.obj);
    % v = value(solc{1});
    % Y = value(solc{2});
    % beta = solc{3};
    % K = Y * diag(1./v);
    % p2p = value(beta)/(1-lam_use);

% else
%     out = Dual_SS_Extended_all_noise_lam(sim, d, T, lam_use, type);
% 
%     sol = optimize(out.cons, out.obj, opts);
%     obj_rec = value(out.obj);
%     v = value(out.v);
%     Y = value(out.Y);
%     K = Y * diag(1./v);
%     p2p = value(out.obj)/(1-lam_use);
% end
% 





%% extract solution

t = toc


