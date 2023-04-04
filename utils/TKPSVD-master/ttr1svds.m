function [U,S,V,sigmas]=ttr1svds(A,R)
% [U,S,V,sigmas]=ttr1svds(A,R)
% ----------------------------
% Sparse Tensor Train rank-1 singular value decomposition. Decomposes an arbitrary tensor A into
% a linear combination of orthonormal rank-1 terms. Returns the orthogonal
% vectors U,V and singular values S from each of the SVDs in the TTr1 tree.
% Use the function getAtilde.m to obtain a low rank approximation using the
% U,S,V obtained from this function. Uses svds instead of svd.
%
% U         =   cell, contains the U vectors of each of the SVDs in the
%               tree,
%
% S         =   cell, contains the singular values S of each of the SVDs
%               in the tree,
%
% V         =   cell, contains the V vectors of each of the SVDs in the
%               tree,
%
% sigmas    =   vector, contains the final singular values in the linear
%               combination of rank-1 terms.
%
% A         =   array, d-way array,
%
% R         =   scalar, desired number of terms.
%
% Reference
% ---------
%
% A Constructive Algorithm for Decomposing a Tensor into a Finite Sum of Orthonormal Rank-1 Terms
% http://arxiv.org/abs/1407.1593
%
% 2014, Kim Batselier, Haotian Liu, Ngai Wong
n=size(A);
r=zeros(1,length(n)-1);
for i=1:length(n)-1
    r(i) = min(n(i),R);
end

totalsvd=1;
svdsperlevel=zeros(1,length(r));   % add length 1 for the first level
svdsperlevel(1)=1;  % first level
for i=2:length(r)
    svdsperlevel(i)=prod(r(1:i-1));
    totalsvd=totalsvd+svdsperlevel(i); 
end
nleaf=prod(r);

U=cell(1,totalsvd);
S=cell(1,totalsvd);
V=cell(1,totalsvd);

[Ut,St,Vt]=svds(sparse(reshape(A,[n(1),prod(n(2:end))])),R);
U{1}=Ut;
S{1}=diag(St);
V{1}=Vt;
counter=2; % this counter keeps track during the iterations which V{i} we are computing. This is a linear counter that counts breadth-first
whichvcounter=1;    % this counter keeps track during the iterations of which V we are taking the svd
sigmas=kron(S{1}, ones(nleaf/length(S{1}),1));

tol=prod(n)*eps(max(S{whichvcounter}));


for i=1:length(r)-1           % outer loop over the levels
    Slevel=[];
    for j=1:prod(r(1:i))      % inner loop over the number of eigs for this level 
        if rem(j,r(i)) == 0
            col=r(i);
        else
            col=rem(j,r(i));
        end
        if ~isempty(S{whichvcounter}) && S{whichvcounter}(col) > tol
			[Ut,St,Vt]=svds(sparse(reshape(V{whichvcounter}(:,col),[n(i+1),prod(n(i+2:end))])),R);
			U{counter}=Ut;
			S{counter}=diag(St);
			V{counter}=Vt;
			Slevel=[Slevel;S{counter}];            
        else
            Slevel=[Slevel;zeros(R,1)];
        end
        counter=counter+1;
        if rem(j,length(S{whichvcounter}))==0
            whichvcounter =  whichvcounter+1;
        end
    end
    Slevel=kron(Slevel, ones(nleaf/length(Slevel),1));
    sigmas=sigmas.*Slevel;
end

end