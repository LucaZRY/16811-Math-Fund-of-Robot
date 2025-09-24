function p = newton_eval(x, a, xstar)
% newton_eval: evaluate Newton interpolating polynomial at xstar
%
% Input:
%   x - vector of nodes
%   a - Newton coefficients from divided_diff_coeffs
%   xstar - evaluation point (scalar or vector)
% Output:
%   p - interpolated value(s)
%
% Example:
%   xstar = 1.5;
%   p = newton_eval(x, a, xstar);

    n = length(a);
    p = a(n) * ones(size(xstar));
    for k = n-1:-1:1
        p = a(k) + (xstar - x(k)) .* p;
    end
end
