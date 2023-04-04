function [M0,id]=alignEMmatrices(M_true, M0)
% Align the signatures of the extracted endmembers M0
% to the correspoding materials in M_true using a least 
% angle criteria

[L,P] = size(M_true);

id = zeros(P,1);
for k = 1:P
    for l = 1:P
        s(l) = 180*acos( (M_true(:,k).')*M0(:,l) /(norm(M_true(:,k))*norm(M0(:,l))) )/pi;
    end
    [~, id(k)] = min(s);
end
M0 = M0(:,id);

