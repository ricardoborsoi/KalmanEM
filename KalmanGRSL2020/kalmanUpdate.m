function [psi, P, v] = kalmanUpdate(psi, y, B, P, Q, sigma_r)
%function [psi, P] = kalmanUpdate(psi, y, B, P, Q, R)
%[psi, P] = kalmanUpdate(psi, y, B, P, Q, R)
%  
%   Perform the kalman filter update for the state-space model:
%   psi = psi + q
%   y = B*psi + r
%
%   return updated psi and P
%
%   IMPORTANT: matrix R must be a sparse matrix for large problems
% 
%   Authors: Tales Imbiriba and Ricardo A. Borsoi.
%   Aug 22, 2019.

% prediction:
% psi = psi;
P = P + Q;

% % update:
% v = y - B*psi;
% S = B*P*B' + R;
% K = (P*B')/S;
% psi = psi + K*v;
% P = P - K*S*K';

% update:
v = y - B*psi;
% S = B*P*B' + R;
% K = (P*B')/S;

sigma2_r_inv = sigma_r^(-2);
BtRinvB = sigma2_r_inv * (B' * B);
BtSinv = sigma2_r_inv * B' - sigma2_r_inv*BtRinvB*((inv(P) + BtRinvB)\B');

K = P*BtSinv;
psi = psi + K*v;
% P = P - K*S*K';
P = P - (K*B)*P*(B'*K') - sigma_r^2 * (K*K');


end

