Nlam = 100;
% Nlam = 12;
% Nlam = 5;
lam_list = linspace(1, 0, Nlam+1);

load('p2p_ess_result_d2.mat')
p2p_list_2 = p2p_list;
err_list(58:end) = 11; %(ran out of time to execute)
err_list_2 = err_list;


load('p2p_ess_result.mat')
[mm, ii] = min(p2p_list);
[mm, ii] = min(p2p_list);


figure(2)
clf
hold on
scatter(lam_list(err_list==0), p2p_list(err_list==0), 'filled')
scatter(lam_list(err_list_2==0), p2p_list_2(err_list_2==0), 'filled')
scatter(lam_list(ii), mm, 200, 'k')
ylim([4, max(p2p_list(err_list==0))])
set(gca, 'Yscale', 'log')
xlabel('$\lambda$', 'interpreter', 'latex', 'FontSize', 14)
ylabel('$\ell_\infty$ upper-bound', 'interpreter', 'latex', 'FontSize', 14)
title('Worst-Case Peak-to-Peak Regulation', 'Interpreter', 'latex', 'FontSize', 16)
legend({'d=1 bounds', 'd=2 bounds', sprintf('best bound %0.4f', mm)}, 'location', 'southeast')