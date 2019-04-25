function delta = deltabinary(S0,K,r,sigma,T)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    d1 = ((log(S0/K) + (r-0.5*sigma^2)*T)/(sigma*sqrt(T)) + sigma*sqrt(T));
    d2 = d1 - sigma * sqrt(T);
    delta = exp(-r*T)*(1/(sqrt(2*pi)*sigma*S0*sqrt(T)))*exp(-d2^2/2);

end

