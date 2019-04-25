function [mean,error] = MC(n,S0,K,r,sigma)

var = 0;
mean = 0;
term = S0*exp(r-sigma^2/2);
for i=1:n
    payoff = max(0,term*exp(sigma*normrnd(0,1))-K);
    var = var + (i-1)*(payoff-mean)^2/i; 
    mean = mean + (payoff-mean)/i;
end

var = var/(n-1);
std = sqrt(var);
error = 1.96*std/sqrt(n);
mean = mean*exp(-0.06);

end

