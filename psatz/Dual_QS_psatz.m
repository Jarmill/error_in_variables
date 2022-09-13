function [p_psatz, cons_psatz, Gram, Coef] = Dual_QS_psatz(p, C, d, vars, A)
%% Enforce Scherer Psatz using eq.(13,14,33)
%    p:   nonnegative polynomial
%    C:   eq and ineq constraints
%    d:   degree of p
% vars:   variables A,B,DX
%    A:   used for partial info
n_g = length(C.ineq);       % # of inequality
n_h = length(C.eq);         % # of equality
n_var = length(vars);       % # of variables
d_g = degree(C.ineq);       % degree of inequality
d_h = degree(C.eq);         % degree of equality
Gram = cell(1+n_g, 1);      % gram cofficients of zeta, -Q
Coef = cell(n_h, 1);        % coefficients of mu
p_psatz = 0;
n = n_g/2-n_h;
T = n_g/2/n;
s_p = length(p);            % size of polynomial matrix
zeta_p = cell(n_g/2,1);     % zeta corresponding to dx
zeta_n = cell(n_g/2,1);
mat_mu = cell(n_h,1);

% generate -Q
[pow0, ~] = momentPowers(0, n_var, d);   
vect_0 = recovermonoms(pow0, vars);     % basis v(x)
l_0 = length(vect_0);                  
temp = kron(vect_0, eye(s_p));          % eq.(13)
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
for i = 1:n_g/2     % zeta+
    Gram{i+1} = sdpvar(l_zeta*s_p);
    zeta_p{i} = temp'*Gram{i+1}*temp;        
    p_psatz = p_psatz + zeta_p{i}*C.ineq(i);
    cons_psatz = [cons_psatz; (Gram{i+1} >= 0):'Gram_zeta_p'];
end
zeta_p = reshape(zeta_p, n, T);

for i = 1:n_g/2     % zeta-
    Gram{i+1+n*T} = sdpvar(l_zeta*s_p);
    zeta_n{i} = temp'*Gram{i+1+n*T}*temp;
    p_psatz = p_psatz + zeta_n{i}*C.ineq(i+n*T);
    cons_psatz = [cons_psatz; (Gram{i+1+n*T} >= 0):'Gram_zeta_n'];
end
zeta_n = reshape(zeta_n, n, T);

% generate mu
[pow_mu, ~] = momentPowers(0, n_var, 2*d-d_h);
vect_mu = recovermonoms(pow_mu, vars);
l_mu = length(vect_mu);
for i = 1:n_h   
    Coef{i} = sdpvar(l_mu,s_p,s_p,'full');      % l_mu matrices mu
    for j = 1:s_p
        for k = 1:s_p
            mat_mu{i}(j,k) = vect_mu'*Coef{i}(:,j,k);  
        end
    end
    p_psatz = p_psatz + mat_mu{i}*C.eq(i);
end
mat_mu = reshape(mat_mu, n, T-1);

% p = p_psatz and eq.(33d-f)
cons_psatz = [cons_psatz; (coefficients(p - p_psatz, vars) == 0):'eq'];
Amu = 0;
for i = 1:n         % k = 1
    for j = 1:n
        Amu = Amu - mat_mu{j,1}*A(j,i);
    end
    cons_psatz = [cons_psatz; (coefficients(zeta_p{i,1} - zeta_n{i,1} + Amu,vars) == 0):'eq_k=1'];
end

for k = 2:T-1       % k for sample, i for row, j for column
    for i = 1:n
        Amu = 0;
        for j = 1:n
            Amu = Amu - mat_mu{j,k}*A(j,i);
        end
        Amu = Amu + mat_mu{i,k-1};
        cons_psatz = [cons_psatz; (coefficients(zeta_p{i,k} - zeta_n{i,k} + Amu,vars) == 0):'eq_k=2~T-1'];
    end
end

for i = 1:n         % k = T
    cons_psatz = [cons_psatz; (coefficients(zeta_p{i,end} - zeta_n{i,end} + mat_mu{i,end},vars) == 0):'eq_k=T'];
end
end

