function c = asiangeom(S0,K,r,sigma,T,N)
%BS Summary of this function goes here
%   Detailed explanation goes here
    sigmaz = sigma*sqrt((2*N+1)/(6*(N+1)))
    rho = ((r-0.5*sigma^2) + sigmaz^2)/2
    d1 = ((log(S0/K) + (rho+0.5*sigmaz^2)*T)/(sigmaz*sqrt(T)))
    d2 = d1 - sigmaz*sqrt(T)
    c = exp(-r*T)*(S0*exp(rho*T)*cdf('norm',d1,0,1) - K*cdf('norm',d2,0,1));
end

