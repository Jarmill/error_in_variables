function [p_psatz, cons_psatz, Gram, Coef] = Dual_SS_psatz_all_noise(p, C, d, vars, A, B, eps, n, m, T)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this code covers Dual_SS_psatz.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n_g = length(C.ineq);       % # of inequality
n_h = length(C.eq);         % # of equality
n_var = length(vars);       % # of variables
d_g = degree(C.ineq);       % degree of inequality
d_h = degree(C.eq);         % degree of equality
Gram = cell(1+n_g, 1);      % gram coefficients of sigma0, sigmai
Coef = cell(n_h, 1);        % coefficients of phi
p_psatz = 0;
s_x = n*T;                  % size of constraint in dx,du,dw
s_u = m*(T-1);
s_w = n*(T-1);
Gram_px = cell(s_x,1);      % Gram matrices corresponding to dx,du,dw
Gram_nx = cell(s_x,1);
Gram_pu = cell(s_u,1);
Gram_nu = cell(s_u,1);
Gram_pw = cell(s_w,1);
Gram_nw = cell(s_w,1);
zeta_px = zeros(s_x,1,'like',sdpvar);       % zeta corresponding to dx,du,dw
zeta_nx = zeros(s_x,1,'like',sdpvar);
zeta_pu = zeros(s_u,1,'like',sdpvar);
zeta_nu = zeros(s_u,1,'like',sdpvar);
zeta_pw = zeros(s_w,1,'like',sdpvar);
zeta_nw = zeros(s_w,1,'like',sdpvar);
mu = zeros(n_h,1,'like',sdpvar);

% generate -Q
[pow0, ~] = momentPowers(0, n_var, d);  
vect_0 = recovermonoms(pow0, vars);         % basis v(x)
l_0 = length(vect_0);
Gram{1} = sdpvar(l_0);
Q = -vect_0'*Gram{1}*vect_0;
p_psatz = p_psatz - Q;
cons_psatz = (Gram{1} >= 0):'Gram_0';

% generate zeta
if mod(d_g,2) == 0
    [pow_zeta, ~] = momentPowers(0, n_var, floor((2*d-d_g)/2));
else
    [pow_zeta, ~] = momentPowers(0, n_var, floor((2*d-d_g-1)/2));
end
vect_zeta = recovermonoms(pow_zeta, vars);
l_zeta = length(vect_zeta);

if eps(1) ~= 0
    kx = 1;             % indicator of eps(1)~=0
    for i = 1:s_x
        Gram_px{i} = sdpvar(l_zeta);
        Gram_nx{i} = sdpvar(l_zeta);
        zeta_px(i) = vect_zeta'*Gram_px{i}*vect_zeta;
        zeta_nx(i) = vect_zeta'*Gram_nx{i}*vect_zeta;
        p_psatz = p_psatz + zeta_px(i)*C.ineq(i) + zeta_nx(i)*C.ineq(i);
        cons_psatz = [cons_psatz; (Gram_px{i} >= 0):'Gram_zeta_px'];
        cons_psatz = [cons_psatz; (Gram_nx{i} >= 0):'Gram_zeta_nx'];
    end
    zeta_px = reshape(zeta_px, n, T);
    zeta_nx = reshape(zeta_nx, n, T);
else
    kx = 0;
    Gram_px = [];
    Gram_nx = [];
end

if eps(2) ~= 0
    ku = 1;             % indicator of eps(2)~=0
    for i = 1:s_u
        Gram_pu{i} = sdpvar(l_zeta);
        Gram_nu{i} = sdpvar(l_zeta);
        zeta_pu(i) = vect_zeta'*Gram_pu{i}*vect_zeta;
        zeta_nu(i) = vect_zeta'*Gram_nu{i}*vect_zeta;
        p_psatz = p_psatz + zeta_pu(i)*C.ineq(i+kx*2*s_x) + zeta_nu(i)*C.ineq(i+kx*2*s_x);
        cons_psatz = [cons_psatz; (Gram_pu{i} >= 0):'Gram_zeta_pu'];
        cons_psatz = [cons_psatz; (Gram_nu{i} >= 0):'Gram_zeta_nu'];
    end
    zeta_pu = reshape(zeta_pu, m, T-1);
    zeta_nu = reshape(zeta_nu, m, T-1);
else
    ku = 0;
    Gram_pu = [];
    Gram_nu = [];
end

if eps(3) ~=0
    for i = 1:s_w
        Gram_pw{i} = sdpvar(l_zeta);
        Gram_nw{i} = sdpvar(l_zeta);
        zeta_pw(i) = vect_zeta'*Gram_pw{i}*vect_zeta;
        zeta_nw(i) = vect_zeta'*Gram_nw{i}*vect_zeta;
        p_psatz = p_psatz + zeta_pw(i)*C.ineq(i+kx*2*s_x+ku*2*s_u) + zeta_nw(i)*C.ineq(i+kx*2*s_x+ku*2*s_u);
        cons_psatz = [cons_psatz; (Gram_pw{i} >= 0):'Gram_zeta_pw'];
        cons_psatz = [cons_psatz; (Gram_nw{i} >= 0):'Gram_zeta_nw'];
    end
    zeta_pw = reshape(zeta_pw, n, T-1);
    zeta_nw = reshape(zeta_nw, n, T-1);
else
    Gram_pw = [];
    Gram_nw = [];
end
Gram = {Gram;Gram_px;Gram_nx;Gram_pu;Gram_nu;Gram_pw;Gram_nw};

% generate mu
[pow_mu, ~] = momentPowers(0, n_var, 2*d-d_h);
vect_mu = recovermonoms(pow_mu, vars);
l_mu = length(vect_mu);
for i = 1:n_h  
    Coef{i} = sdpvar(l_mu,1);
    mu(i) = vect_mu'*Coef{i};
    p_psatz = p_psatz + mu(i)*C.eq(i);
end
mu = reshape(mu, n, T-1);

% p = p_psatz and eq.(47b-e), not following exactly but leads to same result
cons_psatz = [cons_psatz; (coefficients(p - p_psatz, vars) == 0):'eq'];
if eps(1) ~= 0
    cons_psatz = [cons_psatz; (coefficients(zeta_px(:,1) - zeta_nx(:,1) - A'*mu(:,1), vars) == 0):'eq_A_k=1'];
    for i = 2:T-1
        cons_psatz = [cons_psatz; (coefficients(zeta_px(:,i) - zeta_nx(:,i) - (A'*mu(:,i) - mu(:,i-1)), vars) == 0):'eq_A_k=2~T-1'];
    end
    cons_psatz = [cons_psatz; (coefficients(zeta_px(:,end) - zeta_nx(:,end) + mu(:,end), vars) == 0):'eq_A_k=T'];
end
if eps(2) ~= 0
    for i = 1:T-1
        cons_psatz = [cons_psatz; (coefficients(zeta_pu(:,i) - zeta_nu(:,i) - B'*mu(:,i), vars) == 0):'eq_B'];
    end
end
if eps(3) ~= 0
    for i = 1:T-1
        cons_psatz = [cons_psatz; (coefficients(zeta_pw(:,i) - zeta_nw(:,i) + mu(:,i), vars) == 0):'eq_w'];
    end
end
end
