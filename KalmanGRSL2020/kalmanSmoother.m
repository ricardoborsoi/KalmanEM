function [psis_s, Ps_s] = kalmanSmoother(psis, Ps, Q)
%function [psis_s, Ps_s] = kalmanSmoother(psis, Ps, Q)
%
% The inputs are:  
%   psis is a cell with all the psis (Kalman states) obtained using the 
%       Kalman filter and ordered from t=0 to t=T.
%   
%   Ps is a cell with the the evolution of the state covariance matrix
%   ordered from t=0 to t=T.
%
% The outputs are:
%   psis_s: a cell containing the ordered (t=0 to t=T) smoother state 
%           vectors.
%   Ps_s: a cell containing the ordered state covariance matrices. 
%
%   Authors: Tales Imbiriba and Ricardo Borsoi
%   date of current version: 27/08/2019
% =========================================================================


T = length(psis);

psis_s = cell(1,T);
Ps_s = cell(1,T);

psis_s{T} = psis{T};
Ps_s{T} = Ps{T};

for t=T-1:-1:1
    
    psi_pred = psis{t};
    P_pred = Ps{t} + Q;
    
    G = Ps{t}/P_pred;
    psis_s{t} = psis{t} + G * (psis_s{t+1} - psi_pred);
    Ps_s{t} = Ps{t} + G * (Ps_s{t+1} - P_pred) * G';
    
end

end

