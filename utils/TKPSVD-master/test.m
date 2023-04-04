
N = 1000; P = 5; L = 200;
A = rand(P*L, P*L);
% [B,sigmas]=tkpsvd(A,[P P L L]);
[B,sigmas]=tkpsvd(A,[L L   P P]);


% A= \sum_{j=1}^R \sigmas_j A^{dj} \otimes ... \otimes A^{1j}
    

Atil = zeros(size(A));
for i=1:size(B,2)
    Atil = Atil + sigmas(i) * kron(B{2,i}, B{1,i});
end

norm(A-Atil,'fro')/norm(A,'fro')
