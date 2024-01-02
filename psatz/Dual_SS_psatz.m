function [p_psatz, cons_psatz, Gram, Coef] = Dual_SS_psatz(p, C, d, vars, A, B)
%% Enforce Putinar Psatz using eq.(8,18)
%    p:   nonnegative polynomial
%    C:   eq and ineq constraints
%    d:   degree of p
% vars:   variables A,B,DX
%    A:   used for partial info
n_g = length(C.ineq);       % # of inequality
n_h = length(C.eq);         % # of equality

n_gu = length(C.inequ);


n_var = length(vars);       % # of variables
d_g = degree(C.ineq);       % degree of inequality
d_h = degree(C.eq);         % degree of equality
Gram = cell(1+n_g, 1);      % gram coefficients of sigma0, sigmai
Coef = cell(n_h, 1);        % coefficients of phi
p_psatz = 0;
n = n_g/2-n_h;              % note that n_g = 2nT, n_h = n(T-1)
T = n_g/2/n;
m = n_gu / (2 * (T-1));
zeta_p = zeros(n_g/2,1,'like',sdpvar);       % zeta corresponding to dx
zeta_n = zeros(n_g/2,1,'like',sdpvar);
mu = zeros(n_h,1,'like',sdpvar);


psi_p =  zeros(n_gu/2,1,'like',sdpvar);       % zeta corresponding to dx
psi_n =  zeros(n_gu/2,1,'like',sdpvar);       % zeta corresponding to dx

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

for i = 1:n_g/2     % zeta+
    Gram{i+1} = sdpvar(l_zeta);
    zeta_p(i) = vect_zeta'*Gram{i+1}*vect_zeta;
    p_psatz = p_psatz + zeta_p(i)*C.ineq(i);
    cons_psatz = [cons_psatz; (Gram{i+1} >= 0):'Gram_zeta_p'];
end
zeta_p = reshape(zeta_p, n, T);

for i = 1:n_g/2     % zeta-
    Gram{i+1+n*T} = sdpvar(l_zeta);
    zeta_n(i) = vect_zeta'*Gram{i+1+n*T}*vect_zeta;
    p_psatz = p_psatz + zeta_n(i)*C.ineq(i+n*T);
    cons_psatz = [cons_psatz; (Gram{i+1+n*T} >= 0):'Gram_zeta_n'];
end
zeta_n = reshape(zeta_n, n, T);


use_psi = isfield(C, 'inequ') && ~isempty(C.inequ);
if use_psi
    for i = 1:n_gu/2     % psi+
        Gram{i+1+2*n*T} = sdpvar(l_zeta);
        psi_p(i) = vect_zeta'*Gram{i+1}*vect_zeta;
        p_psatz = p_psatz + psi_p(i)*C.inequ(i);
        cons_psatz = [cons_psatz; (Gram{i+1} >= 0):'Gram_zeta_p'];
    end
    psi_p= reshape(psi_n, m, T-1);
    
    for i = 1:n_gu/2     % psi-
        Gram{i+1+m*(T-1) + 2*n*T} = sdpvar(l_zeta);
        psi_n(i) = vect_zeta'*Gram{i+1+n*T}*vect_zeta;
        p_psatz = p_psatz + psi_n(i)*C.inequ(i+n*(T-1));
        cons_psatz = [cons_psatz; (Gram{i+1+n*T} >= 0):'Gram_zeta_n'];
    end
    psi_n= reshape(psi_n, m, T-1);
end




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

% p = p_psatz and eq.(18c-e)
cons_psatz = [cons_psatz; (coefficients(p - p_psatz, vars) == 0):'eq'];
cons_psatz = [cons_psatz; (coefficients(zeta_p(:,1) - zeta_n(:,1) - A'*mu(:,1), vars) == 0):'eq_k=1'];
for i = 2:T-1
    cons_psatz = [cons_psatz; (coefficients(zeta_p(:,i) - zeta_n(:,i) - (A'*mu(:,i) - mu(:,i-1)), vars) == 0):'eq_k=2~T-1'];
end

if use_psi
    for i = 1:T-1
        cons_psatz = [cons_psatz; (coefficients(psi_p(:,i) - psi_n(:,i) - (B'*mu(:,i)), vars) == 0):'eq_ku=1~T-1'];
    end
end

cons_psatz = [cons_psatz; (coefficients(zeta_p(:,end) - zeta_n(:,end) - (-mu(:,end)), vars) == 0):'eq_k=T'];
end


