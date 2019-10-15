function [Coal,tau,W] = CoalDetect(S_fifo,T_arrival,t,c,Wu,Wl,kappa) 
%*********************************************************************************
%% function "CoalDetect" evolves two bounding FIFO GI/GI/c queues forward in time
%% using the same interarrival times and service times until coalescence or time 0
%% inputs:
    % S_fifo = reordered sequence of service times by initiation time
    % T_arrival = sequence of arrival times from -kappa to -1
    % t = arrival time of (-kappa)-th customer
    % c = number of servers
    % Wu = ordered remaining workload vector of upper bound process
    % Wl = ordered remaining workload vector of lower bound process
    % kappa = inspection index (before 0)
%% outputs:
    % Coal = indicator that tells whether two bouding processes coalesce
    %        before 0 (Coal = 1) or not (Coal = 0)
    % tau = coalensence time of two bounding processes (-Inf if Coal = 0)
    % W = ordered workload vector of the target process if Coal = 1 ([] o.w.)
%**********************************************************************************

    k = 1;
    Coal = 0;
    tau = -Inf;
    W = [];
    while(1)
        if(k > kappa)
            NET = [Inf;Wu(1);Wl(1)];
        else
            NET = [T_arrival(k)-t;Wu(1);Wl(1)];
        end
        [NTm,NTp] = min(NET);
        t = t + NTm;
        if(t>0)
            break;
        end
        Wu = Wu - NTm;
        Wl = Wl - NTm;
        if(NTp==1) % an arrival comes next
            if(Wu(c)==Inf)
                Wu(c) = S_fifo(k);
            else
                Wu(1) = Wu(1) + S_fifo(k);
            end
            Wu = sort(Wu);
            if(Wl(c)==Inf)
                Wl(c) = S_fifo(k);
            else
                Wl(1) = Wl(1) + S_fifo(k);
            end
            Wl = sort(Wl);
            k = k + 1;
        end
        id_u = Wu==0;
        Wu(id_u) = Inf;
        Wu = sort(Wu);
        id_l = Wl==0;
        Wl(id_l) = Inf;
        Wl = sort(Wl);
        if(sum(Wu==Wl)==c)
            tau = t;
            Coal = 1;
            W = Wu;
            break;
        end
    end
    

end


