function [S_fifo,T_arrival,t] = Reorder_Service(c,N,S,T,U,lambda1,k1,lambda2,k2)
%*********************************************************************************
%% function "Reorder_Service" reconstructs the mlti-server queue under RA discipline
%% and reorder the service times in the order of service initiations
%% inputs:
    % c = number of servers
    % N = coalescence time
    % S = sequence of service times
    % T = sequence of interarrival times
    % U = sequence of server assignments
    % lambda1,k1 = interarrival distribution parameters, i.e. T~Erlang(k1,lambda1)
    % lambda2,k2 = service time distribution parameters, i.e. S~Erlang(k2,lambda2)
%% outputs:
    % S_fifo = reordered sequence of service times by initiation time
    % T_arrival = sequence of arrival times from -N to -1
    % t = arrival time of N-th customer (counting backwards in time)
%*********************************************************************************
    
    t = -sum(T(1:N)); % current time stamp (stopping time)
    T_fifo = fliplr(T(2:N)); % forward in time
    T_arrival = t + [0,cumsum(T_fifo)]; % arrival times forward in time before 0
    S_target = fliplr(S(1:N)); % service times of interest and forward in time
    U_target = fliplr(U(1:N)); % server assignment of interest and forward in time
    S_InitTime = zeros(1,N); % vector to store the start time of each service
    % consider for each queue
    Delay = zeros(c,1); % record delay time of last customer with server i before 0
    Arr = zeros(c,1); % record arrival time of last customer with server i before 0
    Svc = zeros(c,1); % record service time of last customer with server i before 0
    Svc_InitTime = zeros(c,1); % record service initiation time of last customer with server i before 0
    for i = 1:c
        id = find(U_target==i); %find arrivals directed to server i
        if(~isempty(id))
            W = 0; %waiting time of first arrival of queue i
            S_InitTime(id(1)) = T_arrival(id(1)); % start service immediately
            for k=2:length(id)
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
    % construct RA queues forward in time until enough service times are
    % obtained for the FIFO queue
    lastInitTime = max(S_InitTime);
    if(lastInitTime>=0)
        % for each single server queue, simulate forward in time until
        % service initiation time is greater than lastInitTime
        
        % for the first customer arriving after time 0 (T ~ Ae)
        T_arrival = [T_arrival,-1/lambda1*sum(log(rand(randi(k1),1)))];
        U_target = [U_target,randi(c)];
        S_target = [S_target,-1/lambda2*sum(log(rand(k2,1)))];
        Delay(U_target(end)) = max(Delay(U_target(end))+Svc(U_target(end))-(T_arrival(end)-Arr(U_target(end))),0)*(Arr(U_target(end))~=-Inf);
        Arr(U_target(end)) = T_arrival(end);
        Svc(U_target(end)) = S_target(end);
        S_InitTime = [S_InitTime,Arr(U_target(end))+Delay(U_target(end))];
        Svc_InitTime(U_target(end)) = S_InitTime(end);
        
        while(sum(Svc_InitTime>=lastInitTime)<c)
            T_arrival = [T_arrival,T_arrival(end)-1/lambda1*sum(log(rand(k1,1)))];
            U_target = [U_target,randi(c)];
            S_target = [S_target,-1/lambda2*sum(log(rand(k2,1)))];
            Delay(U_target(end)) = max(Delay(U_target(end))+Svc(U_target(end))-(T_arrival(end)-Arr(U_target(end))),0)*(Arr(U_target(end))~=-Inf);
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
    % before time 0
    T_arrival = T_arrival(1:N);
    S_fifo = S_fifo(1:N);
end

