function [S_fifo,T_arrival,t,R,D,Gamma,S,T,U] = Reorder_Service(c,S,T,U,R,D,Gamma,m,theta,kappa,lambda1,k1,lambda2,k2)
%*********************************************************************************
%% function "Reorder_Service" reconstructs the mlti-server queue under RA discipline
%% and reorder the service times in the order of service initiations
%% inputs:
    % c = number of servers
    % N = coalescence time
    % S = sequence of service times
    % T = sequence of interarrival times
    % U = sequence of server assignments
    % R = c-dimensional reversed random walk
    % D = index sequence of downward milestone events
    % Gamma = index sequence of upward milestone events
    % m = layer height
    % theta = exponential tilting parameter
    % kappa = index of inspection time
    % lambda1,k1 = interarrival distribution parameters, i.e. T~Erlang(k1,lambda1)
    % lambda2,k2 = service time distribution parameters, i.e. S~Erlang(k2,lambda2)
%% outputs:
    % S_fifo = reordered sequence of service times by initiation time
    % T_arrival = sequence of arrival times from -N to -1
    % t = inspection time
    % R = updated c-dimensional reversed random walk
    % D = updated index sequence of downward milestone events
    % Gamma = updated index sequence of upward milestone events
    % S = updated sequence of service times
    % T = updated sequence of interarrival times
    % U = updated sequence of server assignments
%*********************************************************************************

    %-- simulate further backward in time to get the correct service time
    % following vector store the first index for each single server queue
    % that it first empties before the inspection time t_{\kappa}
    single_server_empty_index = zeros(c,1); 
    curr = kappa;
    while(prod(single_server_empty_index)==0)
        if(curr<=D(end-1))
            V = max(R(:,curr+1:end)')' - R(:,curr+1);
            if(sum(V==0)>0)
                single_server_empty_index(V==0) = single_server_empty_index(V==0)+curr*(single_server_empty_index(V==0)==0);
            end
        else
            % simulate further backwards
            f = @(R,D,Gamma,S,T,U) (D(end-1)>=curr);
            K = 1;
            [R,D,Gamma,S,T,U] = Multi_Dim_RW(lambda1,k1,lambda2,k2,c,m,theta,f,K,R,D,Gamma,S,T,U);
        end
        curr = curr + 1;
    end
    
    index = [];
    for i = 1:c
        index = [index,find(U(1:single_server_empty_index(i))==i)];
    end
    index = sort(index);
    
    T_fifo = -cumsum(T(1:index(end)));
    T_arrival = fliplr(T_fifo(index));
    S_target = fliplr(S(index));
    U_target = fliplr(U(index));
    S_InitTime = zeros(1,length(index)); % store the start time of each service time
    % consider for each queue
    Delay = zeros(c,1); % record delay time of last customer with server i before 0
    Arr = zeros(c,1); % record arrival time of last customer with server i before 0
    Svc = zeros(c,1); % record service time of last customer with server i before 0
    Svc_InitTime = zeros(c,1); % record service initiation time of last customer with server i before 0
    for i = 1:c
        id = find(U_target==i);
        if(~isempty(id))
            % first arrival of each single server queue i does not wait
            W = 0;
            S_InitTime(id(1)) = T_arrival(id(1)) + W;
            for k = 2:length(id)
                W = max(W+S_target(id(k-1))-(T_arrival(id(k))-T_arrival(id(k-1))),0);
                S_InitTime(id(k)) = T_arrival(id(k)) + W;
            end
            Delay(i) = W;
            Arr(i) = T_arrival(id(end));
            Svc(i) = S_target(id(end));
            Svc_InitTime(i) = S_InitTime(id(end));
        else
            Arr(i) = -Inf;
        end
    end
    % construct RA queues forward in time until enough service times are obtained for the FIFO queue until 0
    lastInitTime = max(S_InitTime);
    if(lastInitTime>=0)
        % for each single server queue, simulate forward in time until
        % service initiation time is greater than lastInitTime
        
        % for the first customer arriving after time 0
        T_arrival = [T_arrival,-1/lambda1*sum(log(rand(randi(k1),1)))]; % Te ~ Ae
        U_target = [U_target,randi(c)];
        S_target = [S_target,-1/lambda2*sum(log(rand(k2,1)))];
        Delay(U_target(end)) = max(Delay(U_target(end))+Svc(U_target(end))-(T_arrival(end)-Arr(U_target(end))),0)*(Arr(U_target(end)~=-Inf));
        Arr(U_target(end)) = T_arrival(end);
        Svc(U_target(end)) = S_target(end);
        S_InitTime = [S_InitTime,Arr(U_target(end))+Delay(U_target(end))];
        Svc_InitTime(U_target(end)) = S_InitTime(end);
        
        while(sum(Svc_InitTime>=lastInitTime)<c)
            T_arrival = [T_arrival,T_arrival(end)-1/lambda1*sum(log(rand(k1,1)))];
            U_target = [U_target,randi(c)];
            S_target = [S_target,-1/lambda2*sum(log(rand(k2,1)))];
            Delay(U_target(end)) = max(Delay(U_target(end))+Svc(U_target(end))-(T_arrival(end)-Arr(U_target(end))),0)*(Arr(U_target(end)~=-Inf));
            Arr(U_target(end)) = T_arrival(end);
            Svc(U_target(end)) = S_target(end);
            S_InitTime = [S_InitTime,T_arrival(end)+Delay(U_target(end))];
            Svc_InitTime(U_target(end)) = S_InitTime(end);
        end
    end
    
    % reorder the service times by initiation time under RA
    [~,Order] = sort(S_InitTime);
    S_fifo = S_target(Order);
    
    % only return arrival times and service times of customers arriving
    % between inspection time and 0
    S_fifo = S_fifo(T_arrival<0);
    S_fifo = S_fifo(end-(kappa-1):end);
    T_arrival = T_arrival(T_arrival<0);
    T_arrival = T_arrival(end-(kappa-1):end);
    t = -sum(T(1:kappa));
end

