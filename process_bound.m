function bnd = process_bound(sim, opts)
%% find a solution (assuming process noise alone)
n = size(sim.X_noise, 1);
m = size(sim.U_noise, 1);
T = size(sim.U_noise, 2);

eps_w = sdpvar(1);
A = sdpvar(n, n, 'full');
B = sdpvar(n, m, 'full');
R = zeros(n, T, 'like', sdpvar);
cons = 0;

for i = 1:T
    r = sim.X_noise(:, i+1) - A*sim.X_noise(:, i) - B*sim.U_noise(:, i);
    R(:, i) = r;
    cons = [cons; norm(r, 'inf') <= eps_w];
end

optimize(cons, eps_w, opts)
bnd = value(eps_w);

end
