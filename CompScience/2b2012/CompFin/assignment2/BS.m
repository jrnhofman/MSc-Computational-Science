function c = BSsq(S0,K,r,sigma,T)
%BS Summary of this function goes here
%   Detailed explanation goes here
    d1 = ((log(S0/K) + (r-0.5*sigma^2)*T)/(sigma*sqrt(T)) + sigma*sqrt(T));
    d2 = d1 - sigma * sqrt(T);
    c = S0*cdf('norm',d1,0,1) - K*exp(-r*T)*cdf('norm',d2,0,1);
end

