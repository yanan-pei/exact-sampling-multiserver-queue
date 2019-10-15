# exact-sampling-multiserver-queue
Here are the MATLAB codes for article   
EXACT SAMPLING FOR SOME MULTI-DIMENSIONAL QUEUEING MODELS WITH RENEWAL INPUT  
by Jose Blanchet, Yanan Pei, Karl Sigman (Oct 14, 2019)  
   
The codes perform exact simulation for Erlang (k1, λ)/ Erlang (k2, μ)/c queue:  
Folder "algorithm1" contains codes for Algorithm 1 (Section 4.1);  
Folder "algorithm2_sandwich" contains codes for the more effcient Algorithm 2 (Section 4.2).  
  
To run the algorithms, you first set  
lambda1 = λ, lambda2 = μ, k1 = k1, k2 = k2, c = c.  
To run Algorithm 1:  
[Queue,RS,NumCus,t] = ExactSim(lambda1,k1,lambda2,k2,c)  
To run Algorithm 2:  
[tau,Queue,RS,NumCus,TimeElapse,NumArrSimulated,T_inspect] = EffExactSim(lambda1,k1,lambda2,k2,c),   
  
where   
tau = coalensence time of two bounding processes,  
Queue = vector of service requirements of customers waiting in queue at time 0,  
RS = remaining service time of c servers at time 0 (+Inf if anyone is idle),  
NumCus = number of customers in the system at time 0,  
TimeElapse = elapsed time from last arrival before time 0 to time 0,  
NumArrSimulated = number of interarrival times simulated backward in time,  
T_inspect = first successful inspection time that yields coalescence,  
t = first time the RA system is empty.  
