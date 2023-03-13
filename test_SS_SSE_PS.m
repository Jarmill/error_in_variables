clear all
clc
tic
load data_eps0.05-0.14.mat
type = 'no_prior';
obj = '[]';
d = 1;      % order of psatz
opts = sdpsettings('solver', 'mosek', 'verbose', 0);
tol = 1e-6;     % tolerance
n1 = 4;     % noise level     [0.05,0.08,0.11,0.14]
n2 = 100;     % # of samples
success_eps = zeros(n1,n2);
success_T = zeros(n1,n2);
range = 8:2:14;
Acl_eig = cell(n1,n2);

sysd=struct('A', [0.6863,    0.3968;0.3456    1.0388], ...
 B = [0.4170, 0.0001; 0.7203, 0.3023]);


% for different noise
T = 8;
for i = 1:n1
    for j = 1:n2
        yalmip('clear')
        sim = data{i,j};
%         % for SS
%         out = Dual_SS(sim, d, T, type, obj);
%         sol = optimize(out.cons, out.obj, opts);
%         K = value(out.K);
        
%         % for SSE
%         out = Dual_SS_Extended(sim, d, T, type, obj);
%         sol = optimize(out.cons, out.obj, opts);
%         v = value(out.v);
%         Y = value(out.Y);
%         K = Y * diag(1./v);

        % for PS
%         sysd=[];
        out = Dual_Pos(sim, d, T, sysd, type, obj);
        sol = optimize(out.cons, out.obj, opts);
        v = value(out.v);
        Y = value(out.Y);
        K = Y * diag(1./v);

        Acl_eig{i,j} = eig(sysd.A+sysd.B*K);
        if sol.problem==0 && max(abs(Acl_eig{i,j})) <= 1 - tol
            success_eps(i,j) = success_eps(i,j) + 1;
        end
        fprintf('%d: \t\t eps %d out of %d: %d\n', i, j, n2, sol.problem)
    end
    save('ps_experiment_eps.mat', 'success_eps');
    fprintf('eps %d out of %d\n', i, n1)
end

% % for different # of samples
% for i = 1:n1
%     for j = 1:n2
%         yalmip('clear')
%         sim = data{i,j};
%         T = range(i);
%         
%         % for SS
%         out = Dual_SS(sim, d, T, type, obj);
%         sol = optimize(out.cons, out.obj, opts);
%         K = value(out.K);
%         
% %         % for SSE
% %         out = Dual_SS_Extended(sim, d, T, type, obj);
% %         sol = optimize(out.cons, out.obj, opts);
% %         v = value(out.v);
% %         Y = value(out.Y);
% %         K = Y * diag(1./v);
%         
% %         % for PS
% %         out = Dual_Pos(sim, d, T, sysd, type, obj);
% %         sol = optimize(out.cons, out.obj, opts);
% %         v = value(out.v);
% %         Y = value(out.Y);
% %         K = Y * diag(1./v);
% 
%         Acl_eig{i,j} = eig(sysd.A+sysd.B*K);
%         if sol.problem==0 && max(abs(Acl_eig{i,j})) <= 1 - tol
%             success_T(i,j) = success_T(i,j) + 1;
%         end
%         fprintf('%d: \t\t eps %d out of %d: %d\n', i, j, n2, sol.problem)
%     end
% %     save('extended_experiment_T.mat', 'success_T');
%     fprintf('T %d out of %d\n', i, n1)
% end

toc

sum(success_eps,2)
sum(success_T,2)