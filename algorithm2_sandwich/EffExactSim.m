function [tau,Queue,RS,NumCus,TimeElapse,NumArrSimulated,T_inspect]...
    = EffExactSim(lambda1,k1,lambda2,k2,c)
%*********************************************************************************
%% function "EffExactSim_Phasetype" does perfect sampling for Erlang/Erlang/c queues
%% inputs:
    % lambda1,k1 = interarrival distribution parameters, i.e. T~Erlang(k1,lambda1)
    % lambda2,k2 = service time distribution parameters, i.e. S~Erlang(k2,lambda2)
    % c = number of servers
    % NOTE: require lambda1/k1 < c*lambda2/k2
%% outputs:
    % tau = coalensence time of two bounding processes
    % Queue = vector of service requirements of customers waiting in queue at time 0
    % RS = remaining service time of c servers at time 0 (+Inf if anyone is idle)
    % NumCus = number of customers in the system at time 0
    % TimeElapse = elapsed time from last arrival before time 0 to time 0
    % NumArrSimulated = number of interarrival times simulated backward in time
    % T_inspect = first successful inspection time that yields coalescence
%**********************************************************************************

%--- find proper theta* and layer height m ---%
% use Newton's method to find exponential tilting parameter theta
% solve theta for equation
% lambda1^k1*lambda2^k2+(c-1)*lambda^k1*(lambda2-theta)^k2=c*(lambda2-theta)^k2*(lambda+theta)^k1
theta = lambda2/2;
theta_new = +Inf;
while(abs(theta-theta_new)>=10^(-6))
    theta_new = theta;
    numerator = c*(lambda2-theta)^k2*(lambda1+theta)^k1-(c-1)*lambda1^k1*(lambda2-theta)^k2-lambda1^k1*lambda2^k2;
    denominator = -c*k2*(lambda2-theta)^(k2-1)*(lambda1+theta)^k1+c*k1*(lambda2-theta)^k2*(lambda1+theta)^(k1-1)+(c-1)*k2*lambda1^k1*(lambda2-theta)^(k2-1);
    theta = theta - numerator/denominator;
end
% find layer height m to satisfy m > log(c)/theta
m = ceil(log(c)/theta)+1;

%--- construct two bounding processes from inspection time ---%   
[R,D,Gamma,S,T,U,Delta,M0] = Global_Max(lambda1,k1,lambda2,k2,c,m,theta,1);
kappa = 10; % start coupling from the (-kappa)-th arrival
f = @(R,D,Gamma,S,T,U) (D(end-1)>=kappa);
K = 1;
[R,D,Gamma,S,T,U] = Multi_Dim_RW(lambda1,k1,lambda2,k2,c,m,theta,f,K,R,D,Gamma,S,T,U);
V = max(R(:,kappa+1:end)')' - R(:,kappa+1); % workload vector of stationary RA GI/GI/c at t_{-\kappa}
Wu = V; 
Wu(Wu==0) = +Inf;
Wu = sort(Wu); % ordered workload vector of upper bounding process
Wl = +Inf*ones(c,1); % workload vector of lower bound process
    
%--- reorder the service times in the order of service initiations ---%
[S_fifo,T_arrival,t,R,D,Gamma,S,T,U] = Reorder_Service(c,S,T,U,R,D,Gamma,m,theta,kappa,lambda1,k1,lambda2,k2);
        
%--- coupling two bounding systems and deteck coalescence before time 0 ---%
[Coal,tau,W] = CoalDetect(S_fifo,T_arrival,t,c,Wu,Wl,kappa);

%--- double kappa and check coalescence for new inspection time if not met ---%
while(~Coal)
    kappa = 2*kappa;
    f = @(R,D,Gamma,S,T,U) (D(end-1)>=kappa);
    K = 1;
    [R,D,Gamma,S,T,U] = Multi_Dim_RW(lambda1,k1,lambda2,k2,c,m,theta,f,K,R,D,Gamma,S,T,U);
    M = max(R(:,kappa+1:end)')';
    V = M - R(:,kappa+1); % workload vector of stationary RA GI/GI/c queue at inspection time
    Wu = V;
    Wu(Wu==0) = +Inf;
    Wu = sort(Wu); % ordered workload vector of upper bound process
    Wl = Inf*ones(c,1); % workload vector of lower bound process
    [S_fifo,T_arrival,t,R,D,Gamma,S,T,U] = Reorder_Service(c,S,T,U,R,D,Gamma,m,theta,kappa,lambda1,k1,lambda2,k2);
    % coupling two bounding systems and deteck coalescence before time 0
    [Coal,tau,W] = CoalDetect(S_fifo,T_arrival,t,c,Wu,Wl,kappa);      
end
T_inspect = -sum(T(1:kappa));

%--- simulate forward from coalescence time ---%
id = T_arrival>=tau;
T_arrival = T_arrival(id);
S_fifo = S_fifo(id);
[Queue,RS,NumCus] = SimForward_FIFO(c,S_fifo,T_arrival,tau,W);
TimeElapse = T(1);
NumArrSimulated = length(T);
    
end

