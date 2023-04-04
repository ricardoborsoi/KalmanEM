
function [A, P_0, Q, psi_0, sigmaR]=estimateLinSysParamEM(P_t,G_t,psi_t,yt,psi_0_old,M0)
% =========================================================================
% Estimates the parameters of the following dnamical system
%       psi_t = psi_{t-1} + q_{t-1}
%         y_t = H * psi_t + r_{t-1}
% 
% comprised of H = A (kron) Il, P_0, psi_0, R = sigma_r^2 I and Q
% ... based on the result of the RTS smoother: P_t, Psi_t, G_t   obtained 
%     from t=0,..,T, or equivalently from t=1,..,T+1
% 
% The inputs are:
%    P_t : cell array of dimension T+1, comprising the estimates from t=0
%            up to t=T, indexed as i=1,..,T+1
%    G_t : cell array of dimension T+1, comprising the estimates from t=0
%            up to t=T, indexed as i=1,..,T+1
%    psi_t : cell array of dimension T+1, comprising the estimates from t=0
%            up to t=T, indexed as i=1,..,T+1
%    yt : cell array of dimension T, comprising the observations from t=1
%            up to t=T
%    psi_0_old : estimative of the initial state returned by EM at the
%            previous iteration
%    L and P : (constants) number of bands and number of endmembers
%    M0 : reference endmember matrix
%    
% 
% The outputs are:
%    [A, P_0, Q, psi_0, sigmaR]
%    NOTES: 
%      a) we are currently not computing R due to the large problem size
%      b) the optimization w.r.t A (abundances) considers R to be a
%         multiple of the identity matrix
%         
% 
%     Authors: Ricardo Borsoi and Tales Imbiriba
%     date of current version: 26/08/2019
% =========================================================================

% addpath('utils/TKPSVD-master')

% Variables we need from the previous EM iteration:
% psi_0_old
% R_old
% H_old

% get the constant values
[L,P] = size(M0);

% order observed image as a vector if necessary
if size(yt{1},2) > 1
    for i=1:length(yt)
        yt{i} = yt{i}(:);
    end
end

% get the variable dimensions
if length(yt) == length(P_t) % discard one sample if the covariances do not include t=0
    T  = length(yt)-1;
else
    T  = length(yt); 
end
LN = size(yt{1},1);
K  = size(psi_t{1},1);


% introduce a synthetic "yt{0}" observation
temp = cell(1,T+1);
for t=2:T+1
    temp{t} = yt{t-1};
end
yt = temp;


% Estimate matrices depending on the smoother variables
Sigma1 = zeros(K,K);
Sigma2 = zeros(K,K);
Sigma3 = zeros(LN,K);
Sigma4 = zeros(K,K);
% Sigma5 = zeros(L,L);

for t=2:T+1
    Sigma1 = Sigma1 + (1/T) * (P_t{t} + psi_t{t} * psi_t{t}');
    Sigma2 = Sigma2 + (1/T) * (P_t{t-1} + psi_t{t-1} * psi_t{t-1}');
    Sigma3 = Sigma3 + (1/T) * (yt{t} * psi_t{t}');
    Sigma4 = Sigma4 + (1/T) * (P_t{t} * G_t{t-1}' + psi_t{t} * psi_t{t-1}');
%     Sigma5 = Sigma5 + (1/T) * (yt{t} * yt{t}');
end

% Compute the solutions to the EM problem ---------------------------------
P_0   = P_t{1} + (psi_t{1}-psi_0_old) * (psi_t{1}-psi_0_old)';
Q     = Sigma1 - Sigma4 - Sigma4' + Sigma2;
psi_0 = psi_t{1};
A     = estimateA(Sigma3,Sigma1,L,P,M0);
% A = estimateA_fmincon(Sigma3,Sigma1,L,P,M0);
% R     = Sigma5 - H_old * Sigma3' - Sigma3 * H_old' + H_old * Sigma1 * H_old';
H = kron(A',eye(L)) * diag(M0(:));
sigmaR = sqrt(((1/T)*norm(cell2mat(yt),'fro')^2 - 2*trace(Sigma3'*H) + trace(Sigma1*H'*H))/(LN));



end

function A = estimateA_fmincon(Sigma3,Sigma1,L,P,M0)

Sigma1 = diag(M0(:)) * Sigma1 * diag(M0(:));
Sigma3 = Sigma3 * diag(M0(:));

N = size(Sigma3,1)/L;

afun = @(A)( trace(kron(A*A',eye(L))*Sigma1) ...
           -2*trace(Sigma3'*kron(A',eye(L))) );

A0 = ones(P,N)/P;

options = optimoptions('fmincon','Display','iter','Algorithm','interior-point',...
    'HessianApproximation','lbfgs');
LB = zeros(size(A0));

% X = fmincon(FUN,X0,A,B,Aeq,Beq,LB,UB,NONLCON,OPTIONS)
A = fmincon(afun,A0,[],[],[],[],LB,[],[],options);


end


function A = estimateA(Sigma3,Sigma1,L,P,M0)
% This estimates H(theta) = A \otimes I by assuming that matrix R is a
% scaled multiple of the identity. It is based on the kronecker
% decomosition of the other matrices involved (Sigma3 and Sigma1).

% get the number of pixels
N = size(Sigma3,1)/L;

% incorporate M0 in the Sigmas
Sigma1 = diag(M0(:)) * Sigma1 * diag(M0(:));
Sigma3 = Sigma3 * diag(M0(:));


% Warning: the function tkpsvd returns a "reversed" decomposition!
% Decompose Sigma1 ------------------------------------
[B,sigmas]=tkpsvd(Sigma1,[L L   P P]);
K1 = size(B,2); % number of terms of the decomposition
D = cell(K1,1);
C = cell(K1,1);

for i=1:K1
    C{i} = sqrt(sigmas(i)) * B{2,i};
    D{i} = sqrt(sigmas(i)) * B{1,i};
end

% Decompose Sigma3 ------------------------------------
[B,sigmas]=tkpsvd(Sigma3,[L L   N P]);
K2 = size(B,2); % number of terms of the decomposition
Dtil = cell(K2,1);
Ctil = cell(K2,1);

for i=1:K2
    Ctil{i} = sqrt(sigmas(i)) * B{2,i};
    Dtil{i} = sqrt(sigmas(i)) * B{1,i};
end

% Invert linear system -------------------------------
Mx = zeros(P,P);
My = zeros(P,N);
for i=1:K1
    Mx = Mx + trace(D{i}) * (C{i} + C{i}');
end

for i=1:K2
    My = My + 2 * trace(Dtil{i}') * (Ctil{i}');
end

A = Mx\My;
end

