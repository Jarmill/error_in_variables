function [p_psatz, cons_psatz, Gram, Coef] = Dual_QS_psatz_all_noise(p, C, d, vars, A, B, eps, n, m, T)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this code covers Dual_QS_psatz.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n_g = length(C.ineq);       % # of inequality
n_h = length(C.eq);         % # of equality
n_var = length(vars);       % # of variables
d_g = degree(C.ineq);       % degree of inequality
d_h = degree(C.eq);         % degree of equality
Gram = cell(1+n_g, 1);      % gram coefficients of sigma0, sigmai
Coef = cell(n_h, 1);        % coefficients of phi
p_psatz = 0;
s_x = n*T;
s_u = m*(T-1);
s_w = n*(T-1);
s_p = length(p);            % size of polynomial matrix
Gram_px = cell(s_x,1);      % Gram matrices corresponding to dx,du,dw
Gram_nx = cell(s_x,1);
Gram_pu = cell(s_u,1);
Gram_nu = cell(s_u,1);
Gram_pw = cell(s_w,1);
Gram_nw = cell(s_w,1);
zeta_px = cell(s_x,1);      % zeta corresponding to dx,du,dw
zeta_nx = cell(s_x,1);
zeta_pu = cell(s_u,1);
zeta_nu = cell(s_u,1);
zeta_pw = cell(s_w,1);
zeta_nw = cell(s_w,1);
mat_mu = cell(n_h,1);

% generate -Q
[pow0, ~] = momentPowers(0, n_var, d);   
vect_0 = recovermonoms(pow0, vars);
l_0 = length(vect_0);
temp = kron(vect_0, eye(s_p));
Gram{1} = sdpvar(l_0*s_p);
Q = -temp'*Gram{1}*temp; 
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
temp = kron(vect_zeta, eye(s_p));

if eps(1) ~= 0
    kx = 1;
    for i = 1:s_x
        Gram_px{i} = sdpvar(l_zeta*s_p);
        Gram_nx{i} = sdpvar(l_zeta*s_p);
        zeta_px{i} = temp'*Gram_px{i}*temp;
        zeta_nx{i} = temp'*Gram_nx{i}*temp;
        p_psatz = p_psatz + zeta_px{i}*C.ineq(i) + zeta_nx{i}*C.ineq(i);
        cons_psatz = [cons_psatz; (Gram_px{i} >= 0):'Gram_zeta_px'];
        cons_psatz = [cons_psatz; (Gram_nx{i} >= 0):'Gram_zeta_mx'];
    end
    zeta_px = reshape(zeta_px, n, T);
    zeta_nx = reshape(zeta_nx, n, T);
else
    kx = 0;
    Gram_px = [];
    Gram_nx = [];
end

if eps(2) ~= 0
    ku = 1;
    for i = 1:s_u
        Gram_pu{i} = sdpvar(l_zeta*s_p);
        Gram_nu{i} = sdpvar(l_zeta*s_p);
        zeta_pu{i} = temp'*Gram_pu{i}*temp;
        zeta_nu{i} = temp'*Gram_nu{i}*temp;
        p_psatz = p_psatz + zeta_pu{i,1}*C.ineq(i+kx*2*s_x) + zeta_nu{i,1}*C.ineq(i+kx*2*s_x);
        cons_psatz = [cons_psatz; (Gram_pu{i} >= 0):'Gram_zeta_pu'];
        cons_psatz = [cons_psatz; (Gram_nu{i} >= 0):'Gram_zeta_mu'];
    end
    zeta_pu = reshape(zeta_pu, m, T-1);
    zeta_nu = reshape(zeta_nu, m, T-1);
else
    ku = 0;
    Gram_pu = [];
    Gram_nu = [];
end

if eps(3) ~= 0
    for i = 1:s_w
        Gram_pw{i} = sdpvar(l_zeta*s_p);
        Gram_nw{i} = sdpvar(l_zeta*s_p);
        zeta_pw{i} = temp'*Gram_pw{i}*temp;
        zeta_nw{i} = temp'*Gram_nw{i}*temp;
        p_psatz = p_psatz + zeta_pw{i,1}*C.ineq(i+kx*2*s_x+ku*2*s_u) + zeta_nw{i,1}*C.ineq(i+kx*2*s_x+ku*2*s_u);
        cons_psatz = [cons_psatz; (Gram_pw{i} >= 0):'Gram_zeta_pw'];
        cons_psatz = [cons_psatz; (Gram_nw{i} >= 0):'Gram_zeta_mw'];
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
    Coef{i} = sdpvar(l_mu,s_p,s_p,'full');
    for j = 1:s_p
        for k = 1:s_p
            mat_mu{i}(j,k) = vect_mu'*Coef{i}(:,j,k);
        end
    end
     p_psatz = p_psatz + mat_mu{i}*C.eq(i);
end
mat_mu = reshape(mat_mu, n, T-1);

% p = p_psatz and eq.(47b-e), not following exactly but leads to same result
cons_psatz = [cons_psatz; (coefficients(p - p_psatz, vars) == 0):'eq'];
if eps(1) ~= 0
    Amu = 0;
    for i = 1:n
        for j = 1:n
            Amu = Amu - mat_mu{j,1}*A(j,i);
        end
        cons_psatz = [cons_psatz; (coefficients(zeta_px{i,1} - zeta_nx{i,1} + Amu,vars) == 0):'eq_A_k=1'];
    end
    
    for k = 2:T-1     % k for sample, i for row, j for column
        for i = 1:n
            Amu = 0;
            for j = 1:n
                Amu = Amu - mat_mu{j,k}*A(j,i);
            end
            Amu = Amu + mat_mu{i,k-1};
            cons_psatz = [cons_psatz; (coefficients(zeta_px{i,k} - zeta_nx{i,k} + Amu,vars) == 0):'eq_A_k=2~T-1'];
        end
    end
    
    for i = 1:n
        cons_psatz = [cons_psatz; (coefficients(zeta_px{i,end} - zeta_nx{i,end} + mat_mu{i,end},vars) == 0):'eq_A_k=T'];
    end
end

if eps(2) ~= 0
    for k = 1:T-1     % k for sample, i for row, j for column
        for i = 1:m
            Bmu = 0;
            for j = 1:n
                Bmu = Bmu - mat_mu{j,k}*B(j,i);
            end
            cons_psatz = [cons_psatz; (coefficients(zeta_pu{i,k} - zeta_nu{i,k} + Bmu,vars) == 0):'eq_B'];
        end
    end
end

if eps(3) ~= 0
    for k = 1:T-1     % k for sample, i for row, j for column
        for i = 1:n
            cons_psatz = [cons_psatz; (coefficients(zeta_pw{i,k} - zeta_nw{i,k} + mat_mu{i,k},vars) == 0):'eq_w'];
        end
    end
end

end
