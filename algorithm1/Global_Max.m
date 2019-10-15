function [R,D,Gamma,S,T,U,Delta,M0] = Global_Max(lambda1,k1,lambda2,k2,c,m,theta,flag)
%*********************************************************************************
%% function "Global_Max" finds the global maximum of a c-dimensional random walk
%% with each step at coordinate i being S*I(U=i)-T, S, U, T independent
%% T ~ Erlang(k1,lambda1), S ~ Erlang(k2,lambda2), U ~ Unif{1,...,c}
%% inputs:
    % lambda1,k1 = interarrival distribution parameters, i.e. T~Erlang(k1,lambda1)
    % lambda2,k2 = service time distribution parameters, i.e. S~Erlang(k2,lambda2)
    % c = number of servers
    % m = layer height
    % theta = exponential tilting parameter
    % flag = 1 if the first equilibrium interarrival has not been sampled
    % flag = 0 otherwise
%% outputs:
    % R = c-dimensional reversed random walk until Delta (matrix c*(n+1))
    % D = index sequence of downward milestone events
    % Gamma = index sequence of upward milestone events
    % S = sequence of service times
    % T = sequence of interarrival times
    % U = sequence of service assignements in RA policy
    % Delta = the first downward milestone that has its corresponding
    %         upward milestone being +Inf
    % M0 = global maximum of the c-dimensional random walk (coordinate-wise)
%*********************************************************************************

%--- initialization ---%
R = zeros(c,1);
n = 0;
D = 0;
Gamma = Inf;
Lv = zeros(c,1);
servers = (1:c)';
S = []; T = []; U = [];

%--- downward and upward patches sampling ---%
while(1)
    % downward patch
    while(sum(R(:,n+1)>=Lv-m))
        u = randi(c); % server assignment
        s = -1/lambda2*sum(log(rand(k2,1)));
        if(flag==1)
            t = -1/lambda1*sum(log(rand(randi(k1),1)));
            flag = 0;
        else
            t = -1/lambda1*sum(log(rand(k1,1)));
        end
        R = [R,R(:,n+1)+s*(u==servers)-t*ones(c,1)];
        n = n+1;
        S = [S,s];
        T = [T,t];
        U = [U,u];
    end
    D = [D,n]; % update downward milestone sequence
    Lv = R(:,n+1); % update level for next downward milestone search
    
    % upward patch
    R2 = zeros(c,1);
    n2 = 0;
    S2 = []; T2 = []; U2 = [];
    index = randi(c); % choose the direction to do exponential tilting
    while(prod(R2(:,n2+1)<=m))
        if(rand < (lambda2/(lambda2-theta))^k2*(lambda1/(lambda1+theta))^k1/c)
            % assign to the tilted direction
            u2 = index;
            s2 = -1/(lambda2-theta)*sum(log(rand(k2,1)));
        else
            % assign to other directions
            u2 = randi(c-1);
            u2 = u2+(u2>=index);
            s2 = -1/lambda2*sum(log(rand(k2,1)));
        end
        
%         % simulate service times independently
%         if(rand < (lambda2/(lambda2-theta))^k2*(lambda1/(lambda1+theta))^k1/c)
%             s2 = -1/(lambda2-theta)*sum(log(rand(k2,1)));
%         else
%             s2 = -1/lambda2*sum(log(rand(k2,1)));
%         end
        
        t2 = -1/(lambda1+theta)*sum(log(rand(k1,1)));
        R2 = [R2,R2(:,n2+1) + s2*(u2==servers) - t2*ones(c,1)];
        n2 = n2 + 1;
        S2 = [S2,s2];
        T2 = [T2,t2];
        U2 = [U2,u2];
    end 
    % generate I ~ Ber(Tm < Inf)
    I = (rand<=c/sum(exp(theta*R2(:,n2+1))));
    if(I==1)
        % accept and update
        R = [R,repmat(R(:,n+1),[1 n2])+R2(:,2:n2+1)];
        U = [U,U2];
        S = [S,S2];
        T = [T,T2];
        n = n + n2;
        Gamma = [Gamma,n];
    else
        Gamma = [Gamma,+Inf];
        Delta = n;
        M0 = max(R')';
        break;
    end
end
    
end