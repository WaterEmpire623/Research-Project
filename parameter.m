function [para]=parameter()
f = 3e9;                        % frequency (9 GHz)
c = 3e8;                        % speed of light
para.lambda = c / f;            % wavelength
%% System parameters
para.d              = c / f / 2;% inter-element spacing for virtual channel representation
para.d_min          = c / f / 2;% allowed minimal inter-element spacing in practice
para.user_num       = 4;        % number of users
para.ant_num        = 4;        % number of movable antennas
para.path_num       = 15;       % number of channel paths
para.sigma_2        = 1;        % noise power
para.power          = 1;        % transmit power
para.alpha          = 0.01;     % precoder regularized factor, 1 for MMSE, 0.01 for ZF
%% Realization
para.monte_carlo    = 2;     % Realization times