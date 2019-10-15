function [N,R,S,T,U,D,Gamma] = CoalDetect(lambda1,k1,lambda2,k2,c,m,theta) 
%*********************************************************************************
%% function "CoalDetect" simulates a stopping time N when the multi-dimensional 
%% process reversed process achieves its running time maximum
%% inputs:
    % lambda1,k1 = interarrival distribution parameters, i.e. T~Erlang(k1,lambda1)
    % lambda2,k2 = service time distribution parameters, i.e. S~Erlang(k2,lambda2)
    % c = number of servers
    % m = layer height
    % theta = exponential tilting parameter
%% outputs:
    % N = coalescence time (or stopping time)
    % R = c-dimensional reversed random walk
    % S = sequence of service times
    % T = sequence of interarrival times
    % U = sequence of service assignments in RA policy
    % D = index sequence of downward milestone events
    % Gamma = index sequence of upward milestone events
%*********************************************************************************
    [R,D,Gamma,S,T,U,~,~] = Global_Max(lambda1,k1,lambda2,k2,c,m,theta,1); 
    N = 0; % record last running time max checking point
    while(1)
        f = @(R,D,Gamma,S,T,U) (D(end-1) >= N);
        K = 1;
        [R,D,Gamma,S,T,U] = Multi_Dim_RW(lambda1,k1,lambda2,k2,c,m,theta,f,K,R,D,Gamma,S,T,U);
        % find running time maximum
        M = max(R(:,N+1:end)')';
        if(sum(M==R(:,N+1))==c)
            break;
        else
            N = N + 1;
        end
    end

end


