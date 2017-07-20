function Q = diffusionRWR(A, maxiter, restartProb)
    tic
	n = size(A,1);
    
	% Add self-edge to isolated nodes
	A = A + diag(sum(A) == 0);

	% Normalize the adjacency matrix
	renorm = @(M) bsxfun(@rdivide, M, sum(M));
	P = renorm(A);
    fprintf('Normalized network ... \n');
    
	% Personalized PageRank
	restart = eye(n);
	Q = eye(n);
	for i = 1 : maxiter                                                                                                                                                 
        %fprintf('Performing diffussion iteration #%d ...\n', i);
		Q_new = (1 - restartProb) * P * Q + restartProb * restart;
		delta = norm(Q - Q_new, 'fro');
		fprintf('Iter #%d. Frobenius norm: %f\n', i, delta);
		Q = Q_new;
		if delta < 1e-6
			fprintf('Converged.\n');			
			break;
        end
    end   
    toc
end