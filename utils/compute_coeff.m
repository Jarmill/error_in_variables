function out = compute_coeff(a,b)
%%% the coefficients of a(x)b(x) is obtained by
%%% d = ca*cb'
%%% e = (d+d')./(eye(t1)+1)
%%% f = e(ind)
s = length(a);
ind = tril(true(s));   % lower triangular index
% an example 
% [1; x1; x2; x3]*[1 x1 x2 x3] =
% [1     x1     x2     x3;
%  x1    x1^2   x1x2   x1x3;
%  x2    x2x1   x2^2   x2x3;
%  x3    x3x1   x3x2   x3^2]
% graded lexicographic order corresponding to the lower triangular index (we use this)
% 1 x1 x2 x3 x1^2 x1x2 x1x3 x2^2 x2x3 x3^2

% graded reverse lexicographic order (used in polynomial function)
% 1 x1 x2 x3 x1^2 x1x2 x2^2 x1x3 x2x3 x3^2
if nargin < 2
    temp = 2*a./(eye(s)+1);
    out = temp(ind);
else
    temp = a*b';
    temp = (temp + temp')./(eye(s)+1);
    out = temp(ind);
end
end



