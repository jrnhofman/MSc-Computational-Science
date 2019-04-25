function c = BSsq(S0,K,r,sigma,T)
%BS Summary of this function goes here
%   Detailed explanation goes here
    d1 = ((log(S0/K) + (r-0.5*sigma^2)*T)/(sigma*sqrt(T)) + sigma*sqrt(T));
    d2 = d1 - sigma * sqrt(T);
    c = -2*K*S0*cdf('norm',d1,0,1) + K*K*exp(-r*T)*cdf('norm',d2,0,1) + S0^2*exp(r*T+sigma^2*T)*cdf('norm',d1+sigma*sqrt(T),0,1);
end

