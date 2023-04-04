function [psis_smoother, Ps_smoother, psis, Ps, A] = KalmanUnmixing(Y, Q, sigma_r, P0, M0, A0, EM_iters,flag_abundanceSimplexConstraint)
%function [psis_smoother, Ps_smoother, psis, Ps ] = KalmanUnmixing(Y, Q, R, P0, M0, A0, EM_iters)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

%Y is a cell of images. 

% If true, project the abundances to the unit simplex at each EM iteration
if nargin < 8
    flag_abundanceSimplexConstraint = false;
end

% get constants
T = length(Y);
[L, Nmembers] = size(M0);

% initial state vector
psi_0 = ones(L, Nmembers);
psi_0 = psi_0(:);

% T+1 positios to also save psi_{t=0}
psis = cell(1,T+1);
Ps = cell(1,T+1);

psis{1} = psi_0;
Ps{1} = P0;

for k=1:EM_iters
    
    % run kalman and smoother using the current set os parameters \theta
    disp('Running Kalman+smoother...')
    
    v_norms = zeros(T,1);
    fprintf('Iteration: ')
    for t=1:T
%         disp(t);
        fprintf('%d ', t)
        if t==1
            psi = psi_0;
            P   = P0;                    
        end
        y = Y{t}(:);
        H_theta = kron(A0', eye(L));
        B = H_theta * diag(M0(:));
        [psi, P, v] = kalmanUpdate(psi, y, B, P, Q, sigma_r);
        psis{t+1} = psi;
        Ps{t+1} = P;    
        v_norms(t) = norm(v);
    end

    [psis_smoother, Ps_smoother] = kalmanSmoother(psis, Ps, Q);



    % Use EM to estimate the system parameters \theta based on the smoother
    % results
    disp('Running EM...')
    
    G_t = cell(1,T);
    for t=1:T
        G_t{t} = Ps{t}/(Ps{t} + Q);    
    end
    [A0, P0, Q, psi_0, sigma_r]=estimateLinSysParamEM(Ps_smoother,G_t,psis_smoother,Y,psi_0,M0);
    
    if flag_abundanceSimplexConstraint == true
%         A0(A0<=0) = 0; 1e-6;
        for ii=1:size(A0,2)
            A0(:,ii) = ProjectOntoSimplex(A0(:,ii), 1);
        end
    end
    
end
A = A0;
end


