function [R,D,Gamma,S,T,U] = Multi_Dim_RW(lambda1,k1,lambda2,k2,c,m,theta,f,K,R,D,Gamma,S,T,U)
%*********************************************************************************
%% function "Multi_Dim_RW" simulates a c-dimensional random walk with each step 
%% at coordinate i being S*I(U=i)-T, where S, U, T are independent
%% T ~ Erlang(k1,lambda1), S ~ Erlang(k2,lambda2), U ~ Unif{1,...,c}
%% until function value of f >= K
%% inputs:
    % lambda1,k1 = interarrival distribution parameters, i.e. T~Erlang(k1,lambda1)
    % lambda2,k2 = service time distribution parameters, i.e. S~Erlang(k2,lambda2)
    % c = number of servers
    % m = layer height
    % theta = exponential tilting parameter
    % f = function for termination condition
    % K = function value threshold for termination condition
    % R = previously sampled c-dimensional random walk
    % D = previously sampled downward milestone index sequence
    % Gamma = previously sampled upward milestone index sequence
    % S = previously sampled service times sequence
    % T = previously sampled interarrival times sequence
    % U = previously sampld service assignment sequence
%% outputs:
    % R = updated c-dimensional random walk
    % D = updated index sequence of downward milestone events
    % Gamma = updated index sequence of upward milestone events
    % S = updated sequence of service times
    % T = updated sequence of interarrival times
    % U = updated sequence of service assignments in RA policy
%*********************************************************************************
    
    while(f(R,D,Gamma,S,T,U) < K)
        % call function "Global_Max"
        [R_temp,D_temp,Gamma_temp,S_temp,T_temp,U_temp,Delta,M0] = Global_Max(lambda1,k1,lambda2,k2,c,m,theta,0);
        if(prod(M0<=m)) % accept the sub-path to the end of previous path
            % update RW path
            R = [R,repmat(R(:,end),[1 Delta])+R_temp(:,2:Delta+1)];
            Gamma = [Gamma,D(end)+Gamma_temp(2:end)];
            D = [D,D(end)+D_temp(2:end)];
            S = [S,S_temp];
            T = [T,T_temp];
            U = [U,U_temp];
        end
    end


end

