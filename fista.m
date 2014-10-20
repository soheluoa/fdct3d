%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% FISTA DEMO:
%%% Author: Sohel Bhuiyan
%%% FISTA can be used as to solve non-linear Optimization Problem. 
%%% Cost Function to Optimize: J = ||y - Ax||_2^2 + lambda |x|_1^1
%%% with respect to x the model parameter. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [q_factor, data_recon, misfit] = fista(d_obs, T, original_data, beta, PARAM)

[Nx, Ny] = size(d_obs);
dt = 0.002;
lambda = [0.5 0.3];  %%%%%% This coefficient will be changed based on the problem of interpolation and nature of seismic data.
init = ones(Nx, Ny);
T_real = T;
Nit = 90; 
PARAM.T = T;

%%%%%%  To Compute the Maximum eigen Value of Operator %%%%%%%%%%%%%%%%%%%
%alpha = 1.1*power_eig(ones(Nx,Ny),PARAM);
alpha = 1.1;

%%%%%%%%%%%%%%%%%End of Data Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Random Sampling %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j=1:1:length(lambda);
       
    [model, misfit] = solver(d_obs, PARAM, lambda(j), alpha, Nit);
    
    %%%[model, misfit] = solver_amp(d_obs, PARAM, beta(j), Nit);
    
    PARAM.T = ones(Nx,Ny);
    
    data_recon = operator_curve(model, PARAM, 1);    
    
    bias_factor = norm(original_data'*data_recon + data_recon'*original_data)/(2*norm(data_recon'*data_recon));
    
    q_factor(j) = 10*log10(norm(original_data(:))^2/norm(bias_factor.*data_recon(:) - original_data(:))^2);
    
    PARAM.T = T_real;
    
    j
    
end

% coeff = fwa2(data_recon, 'q', 'directional');
% C = coeff;
%         
% % Scale J = 3 %%%%
% makecellFig(C, 182, 3, 3, 0.5 );
% 
% % Scale J = 4 %%%%
% makecellFig(C, 686, 2.1, 4, 0.5 );
% 
% % Scale J = 5 %%%%
% makecellFig(C, 1126, 1.0, 5, 0.5 );

% z_axis = -800:6.25/2:800-6.25/2;
% t_axis = (0:1:Nx-1)*dt;
% 
% figure,
% subplot(1,3,1);
% imagesc(z_axis, t_axis, n_final_data.final_data); xlabel('Offset (m)'); ylabel('Time (s)'); 
% subplot(1,3,2);
% imagesc(z_axis, t_axis, d_obs); xlabel('Offset (m)'); ylabel('Time (s)');
% subplot(1,3,3);
% imagesc(z_axis, t_axis, data_recon); xlabel('Offset (m)'); ylabel('Time (s)');
% 

% nfftt = 2^nextpow2(size(d_obs,1));
% nfftx = 2^nextpow2(size(d_obs,2));
% freq = [0:1:nfftt-1]./(2*dt * nfftt);
% k_x = linspace(-0.5, 0.5, nfftx);
% 
% DATA = (abs(fft2(d_obs, nfftt, nfftx)));
% DATA = fftshift(DATA(1:nfftt/2+1,:),2); DATA2 = (abs(fft2(data_recon, nfftt, nfftx)));
% DATA2 = fftshift(DATA2(1:nfftt/2+1,:),2); DATA3 = (abs(fft2(original_data, nfftt, nfftx)));
% DATA3 = fftshift(DATA3(1:nfftt/2+1,:),2); 

% 
% figure,
% subplot(1,3,1);
% imagesc(k_x, freq(1:nfftt/2+1), DATA); xlabel('Normalized wavenumber'); ylabel('Frequency (Hz)');
% subplot(1,3,2);
% imagesc(k_x, freq(1:nfftt/2+1), DATA2); xlabel('Normalized wavenumber'); ylabel('Frequency (Hz)'); 
% subplot(1,3,3);
% imagesc(k_x, freq(1:nfftt/2+1), DATA3); xlabel('Normalized wavenumber'); ylabel('Frequency (Hz)'); 

%%%%%%% Q-factor %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%q_factor = 10*log10(norm(original_data(:))^2/norm(data_recon(:) - original_data(:))^2);

end