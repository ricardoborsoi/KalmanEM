% =========================================================================
% Run methods for synthetic example 2
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
load synth_dataset_ex2.mat

% get constants
[L,nr,nc,T] = size(Y);
P = size(M,2);
N = nr*nc;

% reorder image
Y_time = cell(T,1);
for t=1:T
    Y_time{t} = reshape(Y(:,:,:,t), [L,N]);
end


%% extract endmembers

try 
    load('endmembers_vca_synth_ex2.mat','M0')
    Mth=M;
catch
    Y_concatenated = zeros(L,N*T);
    for t=1:1 %; T
        Y_concatenated(:,(t-1)*N+1:t*N) = Y_time{t};
    end

    M0 = vca(Y_concatenated,'Endmembers',P);
    Mth=M;

    % Sort M0 with respect to real/desired EM signatures to ease the comparison 
    % of estimated abundance maps
    [M0,id] = alignEMmatrices(Mth, M0);

    save('endmembers_vca_synth_ex2.mat','M0')
end


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
% compute metrics and show stuff

RMSE_A_FLCS   = NRMSE_A(A, A_FCLS);
RMSE_A_Kalman = NRMSE_A(A, A_kalman');

fprintf('RMSE_A:\n')
fprintf('FLCS.......... %f\n', RMSE_A_FLCS)
fprintf('Kalman+EM..... %f\n', RMSE_A_Kalman)

RMSE_M_FLCS   = NRMSE_M(M_nt, Mvca);
RMSE_M_Kalman = NRMSE_M(M_nt, M_kalman);

fprintf('RMSE_M:\n')
fprintf('FLCS.......... %f\n', RMSE_M_FLCS)
fprintf('Kalman+EM..... %f\n', RMSE_M_Kalman)

SAM_M_FLCS   = SAM_M(M_nt, Mvca);
SAM_M_Kalman = SAM_M(M_nt, M_kalman);

fprintf('SAM_M:\n')
fprintf('FLCS.......... %f\n', SAM_M_FLCS)
fprintf('Kalman+EM..... %f\n', SAM_M_Kalman)

RMSE_Y_FLCS   = NRMSE_Y(Y_time, Y_hat_fcls);
RMSE_Y_Kalman = NRMSE_Y(Y_time, Y_hat_kalman);

fprintf('RMSE_Y:\n')
fprintf('FLCS.......... %f\n', RMSE_Y_FLCS)
fprintf('Kalman+EM..... %f\n', RMSE_Y_Kalman)

fprintf('TIMES:\n')
fprintf('FLCS.......... %f\n', time_FCLS)
fprintf('Kalman+EM..... %f\n', time_klaman)


%%
% fh = figure;
% [ha, pos] = tight_subplot(1, T, 0.01, 0.1, 0.1);
% for t=1:T
%     Y_cube = reshape(Y{t}',H,W,L);
%     axes(ha(t));
%     imagesc(3*Y_cube(:,:,[32 20 8])), set(gca,'ytick',[],'xtick',[])
% end




%% Plotting abundances

fontSize = 12;

for pp=1:P
    fh = figure;
    [ha, pos] = tight_subplot(3, T, 0.01, 0.1, 0.1);
    for t=1:T
        
        A_cube = permute(A(:,:,:,t),[2,3,1]);
        axes(ha(0*T + t));
        imagesc(A_cube(:,:,pp),[0 1]), set(gca,'ytick',[],'xtick',[])
        
        A_cube = reshape(A_FCLS{t}',nr,nc,P);
        axes(ha(1*T + t));
        imagesc(A_cube(:,:,pp),[0 1]), set(gca,'ytick',[],'xtick',[])

        A_cube = reshape(A_kalman{t}',nr,nc,P);
        axes(ha(2*T + t));
        imagesc(A_cube(:,:,pp),[0 1]), set(gca,'ytick',[],'xtick',[])
    end
    axes(ha(0*T + 1)); ylabel('True','interpreter','latex','fontsize',fontSize)
    axes(ha(1*T + 1)); ylabel('FCLS','interpreter','latex','fontsize',fontSize)
    axes(ha(2*T + 1)); ylabel('Kalman','interpreter','latex','fontsize',fontSize)
    colormap jet
    set(gcf,'PaperType','A3')
end










