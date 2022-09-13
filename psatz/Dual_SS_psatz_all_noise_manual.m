function [cons_psatz, cmu, Gram] = Dual_SS_psatz_all_noise_manual(sim, cp, d, n, m, T)
%% define size/index
eps = sim.epsilon;
n_var = n^2+n*m;
s1 = nchoosek(n_var+d,d);
s2 = nchoosek(n_var+2*d,2*d);
s_x = n*T;
s_u = m*(T-1);
s_w = n*(T-1);
Gram_px = cell(s_x,1);                   % Gram matrices corresponding to dx,du,dw
Gram_nx = cell(s_x,1);
Gram_pu = cell(s_u,1);
Gram_nu = cell(s_u,1);
Gram_pw = cell(s_w,1);
Gram_nw = cell(s_w,1);
c_px = zeros(s2,s_x,'like',sdpvar);      % coefficient corresponding to dx,du,dw
c_nx = zeros(s2,s_x,'like',sdpvar);
c_pu = zeros(s2,s_u,'like',sdpvar);
c_nu = zeros(s2,s_u,'like',sdpvar);
c_pw = zeros(s2,s_w,'like',sdpvar);
c_nw = zeros(s2,s_w,'like',sdpvar);

% coefficient of -Q
Gram{1} = sdpvar(s1);
cons_psatz = [(Gram{1} >= 0):'Gram0'];
c0 = compute_coeff(Gram{1});     % coefficients of simga0

% coefficient of zeta
c_all = 0;
if eps(1) ~= 0
    for i = 1:s_x    % zeta_x
        Gram_px{i} = sdpvar(s1);
        Gram_nx{i} = sdpvar(s1);
        c_px(:,i) = compute_coeff(Gram_px{i});
        c_nx(:,i) = compute_coeff(Gram_nx{i});
        cons_psatz = [cons_psatz; (Gram_px{i} >= 0):'Gram_zeta_px'];
        cons_psatz = [cons_psatz; (Gram_nx{i} >= 0):'Gram_zeta_mx'];
    end
    c_all = c_all + eps(1)*sum(c_px+c_nx,2);
else
    Gram_px = [];
    Gram_nx = [];
end

if eps(2) ~= 0
    for i = 1:s_u    % zeta_u
        Gram_pu{i} = sdpvar(s1);
        Gram_nu{i} = sdpvar(s1);
        c_pu(:,i) = compute_coeff(Gram_pu{i});
        c_nu(:,i) = compute_coeff(Gram_nu{i});
        cons_psatz = [cons_psatz; (Gram_pu{i} >= 0):'Gram_zeta_pu'];
        cons_psatz = [cons_psatz; (Gram_nu{i} >= 0):'Gram_zeta_mu'];
    end
    c_all = c_all + eps(2)*sum(c_pu+c_nu,2);
else
    Gram_pu = [];
    Gram_nu = [];
end

if eps(3) ~= 0
    for i = 1:s_w    % zeta_w
        Gram_pw{i} = sdpvar(s1);
        Gram_nw{i} = sdpvar(s1);
        c_pw(:,i) = compute_coeff(Gram_pw{i});
        c_nw(:,i) = compute_coeff(Gram_nw{i});
        cons_psatz = [cons_psatz; (Gram_pw{i} >= 0):'Gram_zeta_pw'];
        cons_psatz = [cons_psatz; (Gram_nw{i} >= 0):'Gram_zeta_mw'];
    end
    c_all = c_all + eps(3)*sum(c_pw+c_nw,2);
else
    Gram_pw = [];
    Gram_nw = [];
end
Gram = {Gram;Gram_px;Gram_nx;Gram_pu;Gram_nu;Gram_pw;Gram_nw};

% coefficient of h
h = zeros(s1,s_w);
h(1,:) = sim.X_noise(1+n:end);
h(2:1+n^2,:) = -kron(sim.X_noise(:,1:end-1),eye(n));
h(2+n^2:1+n^2+n*m,:) = -kron(sim.U_noise,eye(n));

% coefficient of mu*h
cmuh = 0;                        % coefficients of muh
cmu = sdpvar(s1,s_w,'full');     % coefficients of mu
for i = 1:s_w
    cmuh = cmuh + compute_coeff(cmu(:,i),h(:,i));
end

% p = p_psatz and eq.(18c-e)
cons_psatz = [cons_psatz; (cp == c0 + c_all + cmuh):'eq'];
if eps(1) ~= 0
    cA = zeros(s1,n^2,'like',sdpvar);     % coefficients of A
    cA(2:n^2+1,:) = eye(n^2);
    cAmu = zeros(s2,s_w,'like',sdpvar);   % coefficients of Amu
    for j = 1:n
        for i = 1:n
            cAmu(:,j) = cAmu(:,j) + compute_coeff(cA(:,n*(j-1)+i),cmu(:,i));
        end
        cons_psatz = [cons_psatz; (c_px(:,j) - c_nx(:,j) == cAmu(:,j)):'eq_A_k=1'];
    end
    
    for j = 1+n:s_x-n
        for i = 1:n
            if mod(j,n) == 0
                k = (n-1)*n;
            else
                k = (mod(j,n)-1)*n;
            end
            cAmu(:,j) = cAmu(:,j) + compute_coeff(cA(:,k+i),cmu(:,n*ceil(j/n)-n+i));
        end
        cons_psatz = [cons_psatz; (c_px(:,j) - c_nx(:,j) == cAmu(:,j) - [cmu(:,j-n);zeros(s2-s1,1)]):'eq_A_k=2~T-1'];
    end
    
    cons_psatz = [cons_psatz; (c_px(:,end-n+1:end) - c_nx(:,end-n+1:end) == -[cmu(:,end-n+1:end);zeros(s2-s1,n)]):'eq_A_k=T'];
end

if eps(2) ~= 0
    cB = zeros(s1,n*m,'like',sdpvar);     % coefficients of B
    cB(2+n^2:end,:) = eye(n*m);
    cBmu = zeros(s2,s_w,'like',sdpvar);   % coefficients of Bmu
    for j = 1:s_u
        for i = 1:n
            if mod(j,m) == 0
                k = (m-1)*n;
            else
                k = (mod(j,m)-1)*n;
            end
            cBmu(:,j) = cBmu(:,j) + compute_coeff(cB(:,k+i),cmu(:,n*ceil(j/m)-n+i));
        end
        cons_psatz = [cons_psatz; (c_pu(:,j) - c_nu(:,j) == cBmu(:,j)):'Bmut'];
    end
end

if eps(3) ~= 0
    for j = 1:s_w
        cons_psatz = [cons_psatz; (c_pw(:,j) - c_nw(:,j) == [cmu(:,j);zeros(s2-s1,1)]):'wmut'];
    end
end
end