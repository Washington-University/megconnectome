
function csvv = hcp_cs2csvv(CSf,W,ndim)

if not(ndim==3)
    error('ndim must be = 3');
end

[U,S,V] = svd(CSf);

[d1,d2] = size(W);
nvox = d1/ndim;
ord_ind = reshape((repmat([1:nvox]',1,ndim) + repmat([0 nvox 2*nvox],nvox,1))',1,[]);
W = W(ord_ind,:);

sv=sparse(diag(S));
iL = 1:d1;
jL = reshape(padarray(1:nvox,2,'replicate','post'),1,[]);

csvv = sparse(d1,d1);
for i=1:d2
    L = sparse(iL,jL,W*U(:,i),d1,nvox);
    R = sparse(jL,iL,(V(:,i)')*(W'),nvox,d1);
    csvv = csvv + sv(i).*(L*R);
end

return