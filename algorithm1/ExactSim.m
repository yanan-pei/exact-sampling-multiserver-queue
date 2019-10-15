function [Queue,RS,NumCus,t] = ExactSim(lambda1,k1,lambda2,k2,c)
%*********************************************************************************
%% function "ExactSim_Phasetype" does perfect sampling for Erlang/Erlang/c queues
%% inputs:
    % lambda1,k1 = interarrival distribution parameters, i.e. T~Erlang(k1,lambda1)
    % lambda2,k2 = service time distribution parameters, i.e. S~Erlang(k2,lambda2)
    % c = number of servers
    % NOTE: require lambda1/k1 < c*lambda2/k2
%% outputs:
    % Queue = vector of service requirements of customers waiting in queue at time 0
    % RS = remaining service time of c servers at time 0 (+Inf if anyone is idle)
    % NumCus = number of customers in the system at time 0
    % TimeElapse = elapsed time from last arrival before time 0 to time 0
    % t = first time the RA system is empty
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

%--- simulate the process backwards in time until hitting 0 ---%
[N,R,S,T,U,D,Gamma] = CoalDetect(lambda1,k1,lambda2,k2,c,m,theta);

%--- reorder the service times in the order of service initiations ---%
[S_fifo,T_arrival,t] = Reorder_Service(c,N,S,T,U,lambda1,k1,lambda2,k2);

%--- simulate the queue forward under FIFO until time 0 ---%
[Queue,RS,NumCus,TimeElapse] = SimForward_FIFO(c,S_fifo,T_arrival,t,N);


end