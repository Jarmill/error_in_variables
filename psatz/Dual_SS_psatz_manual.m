function [cons_psatz, cmu, Gram] = Dual_SS_psatz_manual(sim, cp, d)
% extract size
n = size(sim.X_noise,1);        % dim of state
m = size(sim.U,1);              % dim of input
T = size(sim.X_noise,2);        % # of samples
n_g = n*T;                      % # of inequalities
n_h = n*(T-1);                  % # of equalities
eps = sim.epsilon;
n_var = n^2+n*m;
s1 = nchoosek(n_var+d,d);
s2 = nchoosek(n_var+2*d,2*d);
c_p = zeros(s2,n_g,'like',sdpvar);       % coefficient corresponding to zeta
c_n = zeros(s2,n_g,'like',sdpvar);

% coefficient of -Q
Gram = cell(1+2*n_g, 1);
Gram{1} = sdpvar(s1);
cons_psatz = [(Gram{1} >= 0):'Gram0'];
c0 = compute_coeff(Gram{1});    

% coefficient of zeta
for i = 1:n_g     % zeta+
    Gram{i+1} = sdpvar(s1);
    c_p(:,i) = compute_coeff(Gram{i+1});
    cons_psatz = [cons_psatz; (Gram{i+1} >= 0):'Gram_zeta_p'];
end

for i = 1:n_g     % zeta-
    Gram{i+1+n_g} = sdpvar(s1);
    c_n(:,i) = compute_coeff(Gram{i+1+n_g});
    cons_psatz = [cons_psatz; (Gram{i+1+n_g} >= 0):'Gram_zeta_n'];
end

c_all = sum(c_p+c_n,2);    % coefficients of zeta_plus + zeta_minus

% coefficient of h
ch = zeros(s1,n_h);
ch(1,:) = sim.X_noise(1+n:end);                          % X(k+1)
ch(2:1+n^2,:) = -kron(sim.X_noise(:,1:end-1),eye(n));    % coefficients of  AX(k)
ch(2+n^2:1+n^2+n*m,:) = -kron(sim.U,eye(n));             % coefficients of  BU(k)

% coefficient of mu*h
cmuh = 0;                        % coefficients of muh
cmu = sdpvar(s1,n_h,'full');     % coefficients of mu
for i = 1:n_h
    cmuh = cmuh + compute_coeff(cmu(:,i),ch(:,i));
end

% p = p_psatz and eq.(18c-e)
cons_psatz = [cons_psatz; (cp == c0 + eps*c_all + cmuh):'eq'];
cA = zeros(s1,n^2,'like',sdpvar);     % coefficients of A
cA(2:n^2+1,:) = eye(n^2);
cAmu = zeros(s2,n_h,'like',sdpvar);   % coefficients of Amu

for j = 1:n    
    for i = 1:n
        cAmu(:,j) = cAmu(:,j) + compute_coeff(cA(:,n*(j-1)+i),cmu(:,i));
    end
    cons_psatz = [cons_psatz; (c_p(:,j) - c_n(:,j) == cAmu(:,j)):'eq_k=1'];
end

for j = 1+n:n_h
    for i = 1:n
        if mod(j,n) == 0
            k = (n-1)*n;
        else
            k = (mod(j,n)-1)*n;
        end
        cAmu(:,j) = cAmu(:,j) + compute_coeff(cA(:,k+i),cmu(:,n*ceil(j/n)-n+i));
    end
    cons_psatz = [cons_psatz; (c_p(:,j) - c_n(:,j) == cAmu(:,j) - [cmu(:,j-n);zeros(s2-s1,1)]):'eq_k=2~T-1'];
end

cons_psatz = [cons_psatz; (c_p(:,end-n+1:end) - c_n(:,end-n+1:end) == -[cmu(:,end-n+1:end);zeros(s2-s1,n)]):'eq_k=T'];
end