clear; close all; clc;
addpath("functions");
% rng(0) % the random number generator start from the same starting point
%% System parameters
[para] = parameter(); % System parameters stored in function

Movable_region = [9,16,25,36,49,64,81,100]; % available region for antenna move

CASE = Movable_region;
XX = length(CASE);
% L  = para.path_num;
% K  = para.user_num;
%% Determination
all_rate  = zeros(para.monte_carlo,XX);
mean_rate = zeros(1,XX);
%% Main
for mon = 1:para.monte_carlo
    % disp(mon);
    [beta,phi,theta] = generate_DOA(para); % virtual DoAs of l path of the k user channel
    for cse = 1:XX
        [h] = dictionary_channel(para,beta,phi,theta,CASE(cse));
        H=zeros(CASE(cse),K);
        for k = 1:K
            H(:,k) = dictionary_channel(para,beta(:,k),phi(:,k),theta(:,k),CASE(cse));
        end
        %---------- you should write your algorithm here ----------
        
    end
end
%% Average rate
for i = 1:XX
    mean_rate(i) = mean(all_rate(:,i));      
end
%% Plot