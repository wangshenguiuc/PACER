%function DCA_embedding(NetList)


clear;
maxiter = 20;
restartProb = 0.80;
nnode = 18362;
NetList = {'test.network'};
pathway_file = '';
pathway = dlmread(pathway_file);
for i = 1 : length(NetList)
    
    tic
    fprintf('Processing network %d ...\n', i);
    netID = char(NetList(i));
    ppi_net = loadNetwork(netID, nnode);
    
    
    min(ppi_net(:))
    
    npathway = max(pathway(:,1));
    
    path_net = sparse(pathway(:,1),pathway(:,2),1,npathway,nnode);
    
    all_net = [[ppi_net,path_net'];[path_net,zeros(npathway,npathway)]];
    
    QA = diffusionRWR(all_net, maxiter, restartProb);
    
    size(QA)
    fprintf('RWR done!\n');
    
    alpha = 1 /  (nnode*nnode);
    tic
    QA = log(QA + alpha) - log(alpha);
    toc
    
    for dim=[50,100,500]
        fprintf('svd ...\n');
        [U, S] = svds(QA, dim);
        
        fprintf('computing scores ...\n');
        
        US = U * sqrt(S);
        fprintf('Writing file ...\n');
        dlmwrite([netID,'_net_',num2str(dim),'_',num2str(restartProb),'.U'],U,'delimiter','\t');
        dlmwrite([netID,'_net_',num2str(dim),'_',num2str(restartProb),'.US'],US,'delimiter','\t');
        
    end
    
end