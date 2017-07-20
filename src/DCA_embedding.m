%function DCA_embedding(NetList)
addpath('../Data/');

clear;
maxiter = 20;
restartProb = 0.5;
nnode = 3;
NetList = {'Network.txt'};
dim = 50;
for i = 1 : length(NetList)    
    tic
    fprintf('Processing network %d ...\n', i);
    netID = char(NetList(i));
    ppi_net = loadNetwork(netID, nnode);
    ppi_net = ppi_net*1000;
    
    pathway = dlmread('Pathway_property.txt');
    npathway = max(pathway(:,1));
    
    path_net = sparse(pathway(:,1),pathway(:,2),1,npathway,nnode);
    
    all_net = [[ppi_net,path_net'];[path_net,zeros(npathway,npathway)]];
    
    QA = diffusionRWR(all_net, maxiter, restartProb);
    
    size(QA)
    fprintf('RWR done!\n');
    nnode = size(all_net,1);
    alpha = 1 /  (nnode);
    tic
    QA = log(QA + alpha) - log(alpha);
    toc   
    
    fprintf('svd ...\n');
    [U, S] = svds(QA, dim);
    
    fprintf('Writing file ...\n');
    dlmwrite(['../result/',netID,'_net_',num2str(dim),'_',num2str(restartProb),'.U'],U,'delimiter','\t');
    
end