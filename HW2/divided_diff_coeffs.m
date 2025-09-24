function a = divided_diff_coeffs(x, y)
% divided_diff_coeffs: compute Newton divided difference coefficients
% 
% Input:
%   x - vector of distinct nodes (length n+1)
%   y - vector of function values at x (same length as x)
% Output:
%   a - vector of Newton coefficients for interpolation polynomial
%
% Example:
%   x = [0 1 2];
%   y = [1 3 2];
%   a = divided_diff_coeffs(x, y);

    n = length(x);
    a = y(:); % ensure column vector
    for k = 2:n
        for i = n:-1:k
            a(i) = (a(i) - a(i-1)) / (x(i) - x(i-k+1));
        end
    end
end

