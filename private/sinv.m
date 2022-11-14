function sinv = sinv(s, epsilon)
% Take the psdueoinverse of a diagonal matrix S of singular values
% with s = diag(S).

% Default to system epsilon.
if nargin == 1
    epsilon = eps;
end

% Singular values should be positive or zero.
assert(all(s >= 0));

% Invert non-zero singular values.
sinv = s;
sinv(s > epsilon) = 1./s(s > epsilon);

end