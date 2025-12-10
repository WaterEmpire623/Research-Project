function [h] = dictionary_channel(para,beta,phi,theta,CASE)
lambda = para.lambda;
L = para.path_num;
d = para.d;
NT = CASE;
h=zeros(NT,1);

% combination of channel path
for l=1:L
    h=h+sqrt(1/L)*beta(l)*PW(theta(l),phi(l),d,lambda,sqrt(NT),sqrt(NT));
end

end