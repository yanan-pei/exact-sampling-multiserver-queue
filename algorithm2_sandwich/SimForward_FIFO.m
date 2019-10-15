function [Queue,RS,NumCus] = SimForward_FIFO(c,S_fifo,T_arrival,tau,W)
%% Function Description 
    % this function reconstruct the multi-server queue under FIFO
    % forward in time and gives system status at time 0
%% Inputs:
    % c = number of servers
    % S_fifo = reordered sequence of service times by initiation time
    % T_arrival = sequence of arrival times from -N to -1
    % tau = coalesence time (<0)
    % W = workload vector at coalescence time
%% Outputs:
    % Queue = service requirements of customers waiting in queue, in actual waiting order
    % RS = ordered remainng service time of c servers
    % NumCus = number of customers in the system (both waiting & being served) at time 0
    % TimeElapse = time length since last arrival before 0
    
    % Initialize
    t = tau;
    Queue = []; % service requirements of customers waiting in queue, in actual waiting order
    RS = W; % ordered remaininig service time of c servers
    k = 1; % counter
%     NET = [T_arrival(1)-t,+Inf]; % 2-entry vector to store time to next arrival and time to next service completion
%     
    % evolve the system forward in time
    while(1)
        if(k > length(T_arrival))
            NET = [+Inf,RS(1)];
        else
            NET = [T_arrival(k)-t,RS(1)];
        end
        
        [NTm,NTp] = min(NET); % time to next event and next event type
        t = t + NTm; % update time stamp
        NET = NET - NTm; % udpate time until next events
        RS = RS - NTm; % update sorted remaining service times of c servers
        
        if(t < 0) 
            if(NTp==1) % next event is an arrival
                if(RS(c)==+Inf) % some server is idle
                    RS(c) = S_fifo(k); % assign the new arrival to an idle server
                    RS = sort(RS); % sort RS ascendingly
                else % all servers are busy
                    Queue = [Queue,S_fifo(k)];
                end
                k = k + 1; % udpate index of arrivals
            else % next event is a service completion
                if(isempty(Queue)) % no one is waiting to be served
                    RS(1) = +Inf; % indicate the server is idle
                else % someone is waiting to be served in queue
                    RS(1) = Queue(1); % add the first waiting customer to idle server
                    Queue = Queue(2:end);
                end
                RS = sort(RS);
            end
        elseif(t==+Inf) % system is idle at time 0
            RS = +Inf*ones(c,1);
            break;
        else % some server is busy at time 0
            RS = RS + t; % ordered remaining service at time 0
            break;
        end

    end
    
    NumCus = length(Queue)+sum(RS~=+Inf);
    

end

