function [A_kalman,M_kalman,Y_hat]=adaptor_KalmanEM(Y_cell, A0, M0, nr, nc)
% codes for the paper:
% Kalman fitler and expectation maximization for multitemporal hyperspectral unmixing
% Ricardo A Borsoi, Tales Imbiriba, Pau Closas, José Carlos M Bermudez, Cédric Richard
% IEEE Geoscience and Remote Sensing Letters, 2020
% 
% Inputs:
%   Y_cell : Cell array with imahes at each time instant, L * N
%   A0 : initial reference abundances
%   M0 : reference EM matrix
%   nr, nc : spatial dimension of the images
% Outputs:
%   A_kalman: cell array with abundances
%   M_kalman: cell array with endmembers for each time instant
%   Y_hat: cell array with reconstructed images
% 
% Ricardo Borsoi, 2023

T = length(Y_cell);
[L,P] = size(M0);
Y = Y_cell;

if nr*nc > 500%5000
    flag_decimated = true;
    decimateFact = ceil(sqrt(nr*nc / 5000));
    Y_to_kalman = Y;
    for t=1:T
        Y_to_kalman{t} = reshape(Y_to_kalman{t}', [nr,nc,L]);
        Y_to_kalman{t} = Y_to_kalman{t}(1:decimateFact:end, 1:decimateFact:end, :);
        [ss1,ss2,~] = size(Y_to_kalman{t});
        Y_to_kalman{t} = reshape(Y_to_kalman{t}, [ss1*ss2,L])';
    end
    A0 = reshape(A0', [nr,nc,P]);
    A0 = A0(1:decimateFact:end, 1:decimateFact:end, :);
    [ss1,ss2,~] = size(A0);
    A0 = reshape(A0, [ss1*ss2,P])';
    
else
    flag_decimated = false;
    Y_to_kalman = Y;
end


Q = 0.1*eye(P*L);
%R = sigma*eye(L*N);
sigma_r = 0.01;
P0 = 1*eye(L*P);

M0_kalman = M0;

EM_iters = 5;

% If true, project the abundances to the unit simplex at each EM iteration
flag_abundanceSimplexConstraint = false;

% Reestimate Kalman's abundances after unmixing using FCLS if true:
flag_kalman_reestimate_A = true;

[psis_smoother, Ps_smoother, psis_kalman, Ps_kalman, A_kalman] = KalmanUnmixing(Y_to_kalman, Q, sigma_r, P0, M0_kalman, A0, EM_iters, flag_abundanceSimplexConstraint);

M_kalman = cell(1,T);
for t=1:T
    M_kalman{t} = reshape(psis_smoother{t+1}, L,P) .* M0_kalman;
end



% Re-estimate the abundances using FCLS, or keep them constant
if flag_kalman_reestimate_A == true
    temp = A_kalman;            
    A_kalman = cell(1,T);
    for t=1:T
        % Fully Constrained Least Squares Unmixing (FCLSU) with Kalman results
        alpha_reg = 1e-8;
        if flag_decimated == true % if decimated, A_kalman will have a smaller dimension
            A_kalman{t} = FCLSU(Y{t},M_kalman{t})';
        else
            A_kalman{t} = FCLSU_reg(Y{t},M_kalman{t},temp,alpha_reg)';
        end
    end
else
    temp = A_kalman;
    A_kalman = cell(1,T);
    for t=1:T
        A_kalman{t} = temp;
    end
end


Y_hat = cell(T,1);
for t=1:T
    Y_hat{t} = M_kalman{t} * A_kalman{t};
    Y_hat{t} = reshape(Y_hat{t}', [nr,nc,L]);
end

for t=1:T
    A_kalman{t} = reshape(A_kalman{t}', [nr,nc,P]);
end




