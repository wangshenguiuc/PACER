function A = loadNetwork(fileID, nnode)
	ppi = load(fileID);
    sub_mat = ppi(:,1:2);
    [~,sub_ind] = unique(sub_mat,'rows');
    ppi = ppi(sub_ind,:);
	p1 = ppi(:,1);
	p2 = ppi(:,2);
if size(ppi,2) == 2
A = sparse(p1, p2, 1, nnode, nnode);

else
	s = ppi(:,3);
	s = s ./ 1000;

	A = sparse(p1, p2, s, nnode, nnode);
end
A = max(A,A');
end
