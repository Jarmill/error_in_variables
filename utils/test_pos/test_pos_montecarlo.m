load('pos_traj_std.mat')

% success_time = zeros(length

%% solve positive-stabilization 
type = 'no_prior';
obj = '[]';
% obj = 'p2p';
d = 1;
    
opts = sdpsettings('solver','mosek','verbose', 0);

%% experiments for changing epsilon while T is fixed

success_eps = zeros(length(sim_eps));
for i = 1:length(sim_eps)
    for k = 1:length(sim_eps{i})
        yalmip('clear')
        sys_curr = sim_eps{i}{k};
        T = sys_curr.T;
        out = Dual_Pos(sys_curr, d, sys_curr.T, sysd, type, obj);
        sol = optimize(out.cons, out.obj, opts);

        if sol.problem==0
            success_eps(i) = success_eps(i) + 1;
        end
        fprintf('\t\t eps %d out of %d: %d\n', k, length(sim_eps{i}), sol.problem)
    end    
    save('pos_experiment_stab_std.mat', 'success_eps');
    fprintf('eps %d out of %d\n', i, length(sim_eps))
end

success_T = zeros(length(sim_eps));
for i = 1:length(sim_time)
    for k = 1:length(sim_time{i})
        yalmip('clear')
        sys_curr = sim_time{i}{k};
        T = sys_curr.T;
        out = Dual_Pos(sys_curr, d, sys_curr.T, sysd, type, obj);
        sol = optimize(out.cons, out.obj, opts);

        if sol.problem==0
            success_T(i) = success_T(i) + 1;
        end
        fprintf('\t\t eps %d out of %d: %d\n', k, length(sim_time{i}), sol.problem)
    end
    save('pos_experiment_stab_std.mat', 'success_eps', 'success_T');    
    fprintf('eps %d out of %d\n', i, length(sim_time))
end