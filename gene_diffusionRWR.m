function Q = gene_diffusionRWR(A, maxiter)
tic
n = size(A,1);

% Add self-edge to isolated nodes
A = A + diag(sum(A) == 0);

% Normalize the adjacency matrix
renorm = @(M) bsxfun(@rdivide, M, sum(M));
P = renorm(A);
fprintf('Normalized network ... \n');

% Personalized PageRank
Q = ones(n,1)/n;
for i = 1 : maxiter
    %fprintf('Performing diffussion iteration #%d ...\n', i);
    Q_new = P * Q;
    delta = norm(Q - Q_new, 'fro');
    fprintf('Iter #%d. Frobenius norm: %f\n', i, delta);
    Q = Q_new;
    if delta < 1e-6
        fprintf('Converged.\n');
        break;
    end
end
end