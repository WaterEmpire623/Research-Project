function [beta,phi,theta] =  generate_DOA(para)
L  = para.path_num;
K  = para.user_num;

beta = normrnd(0,1,L,K) + 1i * normrnd(0,1,L,K); % path gains
phi = unifrnd(-1,1,L,K); % BS path DoA-azimuth
theta = unifrnd(-1,1,L,K); % BS path DoA-elevation

% beta  = zeros(L,K); % path gains
% phi   = zeros(L,K); % BS path DoA-azimuth
% theta = zeros(L,K); % BS path DoA-elevation
% 
% for i=1:K
%     beta(:,:,i)  = normrnd(0,1,L,K) + 1i * normrnd(0,1,L,K); 
%     phi(:,:,i)   = unifrnd(-1,1,L,K); 
%     theta(:,:,i) = unifrnd(-1,1,L,K); 
% end

end