function [S,call,delta,loan,final_debt] = hedge_sim(price,strike,interest,delta_t,sigma,sigma2)
%function [S] = hedge_sim(price,strike,interest,delta_t,sigma,sigma2)
% Number of steps %
n = 1/delta_t;

% Creation of matrices %
S = linspace(0,0,n+1);
call = linspace(0,0,n+1);
delta = linspace(0,0,n);
loan = linspace(0,0,n);

% Initial stock price, call price and delta %
S(1) = price;
call(1) = blsprice(S(1),strike,interest,(n+1-1)/n,sigma,0);
delta(1) = blsdelta(S(1),strike,interest,(n+1-1)/n,sigma,0);

% Initial loan: paid for stock - received for call %
loan(1) = delta(1) * S(1) - call(1);

acc_interest = exp(interest * delta_t);
for i=2:n
    
    % Calculate new stock price, new call price and new delta %
    S(i) = S(i-1) + interest * S(i-1) * delta_t + sigma2 * S(i-1) * normrnd(0,sqrt(delta_t));
    call(i) = blsprice(S(i),strike,interest,1 + 1/n - i/n,sigma,0);
    delta(i) = blsdelta(S(i),strike,interest,1 + 1/n - i/n,sigma,0);
    
    % New loan given by old loan with interest + paying new stock -
    % selling old stock %
    loan(i) = loan(i-1) * acc_interest + (delta(i) - delta(i-1)) * S(i);
    
end

% Maturity of the contract, calculate stock price and call price %
S(n+1) = S(n) + sigma * S(n+1-1) * delta_t + sigma2 * S(n+1-1) * normrnd(0,sqrt(delta_t));
call(n+1) = blsprice(S(n+1),strike,interest,(n+1-(n+1))/n,sigma,0); 
    
% Final debt is previous debt with interest - selling old stock + value of call %
final_debt = loan(n) * acc_interest - delta(n) * S(n+1) + max(S(n+1)-strike,0);

end

