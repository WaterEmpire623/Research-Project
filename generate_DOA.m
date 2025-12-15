function [beta,phi,theta] =  generate_DOA(para)
L  = para.path_num;
K  = para.user_num;

beta = normrnd(0,1,L,K) + 1i * normrnd(0,1,L,K); % path gains
phi = unifrnd(-1,1,L,K); % BS path DoA-azimuth
theta = unifrnd(-1,1,L,K); % BS path DoA-elevation

end