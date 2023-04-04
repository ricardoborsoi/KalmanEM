% =========================================================================
% Run methods for synthetic example Tahoe
% 
% Kalman fitler and expectation maximization for multitemporal hyperspectral unmixing
% Ricardo A Borsoi, Tales Imbiriba, Pau Closas, José Carlos M Bermudez, Cédric Richard
% IEEE Geoscience and Remote Sensing Letters, 2020
% =========================================================================

clear
rng(1)
addpath(genpath('utils'))
addpath(genpath('KalmanGRSL2020'))
clc

% load dataset
load('../v1/data/rd_tip_with_outliers_t5.mat')
load('../v1/data/M1.mat')

% get constants
T = 6;
[L,P] = size(M0);
nr = H;
nc = W;
N = nr*nc;

% reorder image
Y_time = cell(T,1);
for t=1:T
    Y{t} = reshape(Y{t}',H,W,L);
    Y_time{t} = reshape(Y{t}, [N,L])';
end



%% extract endmembers

Mth=M0;

% Y_concatenated = zeros(L,N*T);
% for t=1:T
%     Y_concatenated(:,(t-1)*N+1:t*N) = Y_time{t};
% end
% 
% M0 = vca(Y_concatenated,'Endmembers',P);
% 
% % Sort M0 with respect to real/desired EM signatures to ease the comparison 
% % of estimated abundance maps
% [M0,id] = alignEMmatrices(Mth, M0);




%%
disp('FCLSU...')
A_FCLS = cell(1,P);

tic
Mvca = cell(1,T);
for t=1:T
    % extract the EMs
    Mvca{t} = vca(Y_time{t},'Endmembers',P);
    % align the EMs to the real ones to compute the metrics
    [Mvca{t},id] = alignEMmatrices(Mth, Mvca{t});
    
    % Fully Constrained Least Squares Unmixing (FCLSU)
    A_FCLS{t} = FCLSU(Y_time{t},Mvca{t})';
end
time_FCLS = toc;

% compute the reconstructed image
Y_hat_fcls = cell(T,1);
for t=1:T
    Y_hat_fcls{t} = Mvca{t}*A_FCLS{t};
end


%%
% run Kalman+EM

A0 = A_FCLS{1};
tic
[A_kalman,M_kalman,Y_hat_kalman] = adaptor_KalmanEM(Y_time, A0, M0, nr, nc);
time_klaman = toc;

for t=1:T
    A_kalman{t} = reshape(A_kalman{t}, [N,P])';
end

% reorder reconstructed image
for t=1:T
    Y_hat_kalman{t} = reshape(Y_hat_kalman{t}, [N,L])';
end



%%
% compute image resconstruction error and display times

RMSE_Y_FLCS   = NRMSE_Y(Y_time, Y_hat_fcls);
RMSE_Y_Kalman = NRMSE_Y(Y_time, Y_hat_kalman);

fprintf('RMSE_Y:\n')
fprintf('FLCS.......... %f\n', RMSE_Y_FLCS)
fprintf('Kalman+EM..... %f\n', RMSE_Y_Kalman)

fprintf('TIMES:\n')
fprintf('FLCS.......... %f\n', time_FCLS)
fprintf('Kalman+EM..... %f\n', time_klaman)


%% Plot images

fh = figure;
[ha, pos] = tight_subplot(1, T, 0.01, 0.1, 0.1);
for t=1:T
    Y_tmp = reshape(Y_time{t}',nr,nc,L);
    axes(ha(t));
    imagesc(3*Y_tmp(:,:,[32 20 8])), set(gca,'ytick',[],'xtick',[])
end



%% Plotting abundances

fontSize = 12;

for pp=1:P
    fh = figure;
    [ha, pos] = tight_subplot(2, T, 0.01, 0.1, 0.1);
    for t=1:T
        
        A_cube = reshape(A_FCLS{t}',nr,nc,P);
        axes(ha(0*T + t));
        imagesc(A_cube(:,:,pp),[0 1]), set(gca,'ytick',[],'xtick',[])

        A_cube = reshape(A_kalman{t}',nr,nc,P);
        axes(ha(3*T + t));
        imagesc(A_cube(:,:,pp),[0 1]), set(gca,'ytick',[],'xtick',[])
    end
    axes(ha(0*T + 1)); ylabel('FCLS','interpreter','latex','fontsize',fontSize)
    axes(ha(1*T + 1)); ylabel('Kalman','interpreter','latex','fontsize',fontSize)
    colormap jet
end






