function P = com_mat(m, n, fullorsparse)

arguments
    m
    n
    fullorsparse = 'sparse'
end

% determine permutation applied by K
A = reshape(1:m*n, m, n);
v = reshape(A', 1, []);

% apply this permutation to the rows (i.e. to each column) of identity matrix
P = speye(m*n);
P = P(v,:);

if strcmp(fullorsparse, 'full')
    P = full(P);
end

end