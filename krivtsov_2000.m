clear all

T = 100; % time period of interest
num = 100; % number of intervals
dt = T/num; % interval
k = 100000; % number of iterations

beta = 1.5; % shape parameter
eta = 20; % scale parameter

H = zeros(9, num);

for q = 0.1:0.1:0.9 % restoration parameter
    
    N = zeros(k, num);
    
    for i = 1 : k % index of iteration (system)
        
        S = 0; % initialize actual age
        n = 1; % initialize index of consecutive failure
        
        while ~(S > T)
            
            if n == 1
                x = eta*(-log(1-rand))^(1/beta); % inverse of TTFF
                S = S + x; % update actual age
                A = q * x; % update virtual age
            end
            
            if S > T % if the first failure happen outside the period
                break; % of interest, then just break
            end 
            
            N(i, ceil(S/dt)) = N(i, ceil(S/dt)) + 1;
            % record the failure in corresponding interval
            
            n = n + 1; % update the index of failure
            
            x = (A^beta - eta^beta * log(1-rand))^(1/beta) - A;
            % inverse of conditional TBF
            A = q * x + A; % update virtual age
            S = A/q; % update actual age
            
        end      
    end
    
    for i = 1 : k
        for j = 2 : num
            N(i,j) = N(i,j) + N(i,j-1);
        end % accumulate the number of failures for each system
    end
    
    Nm = mean(N,1); % calculate the estimate of GRP
    H(round(q*10),1:num) = Nm(1:num);
    % record the estimate of GRP with different restoration parameter
    
end

%% plot figure

for i = 1 : 9
    plot(1:num,H(i,1:num));
    hold on
end
