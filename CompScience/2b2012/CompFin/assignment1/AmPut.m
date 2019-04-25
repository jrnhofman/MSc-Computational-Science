function [Am_Put,option_price] = AmPut(n,S0,sigma,strike_price,r)

delta_t = 1/(n-1);
u = exp(sigma*sqrt(delta_t));
d = exp(-sigma*sqrt(delta_t));
p = (exp(r*delta_t) - d)/(u - d);

stock_price = zeros(n);
option_price = zeros(n);

% bereken de aandelenprijzen op tijdstip T

for i = 1: 1: n
    for j=1: 1: i
        stock_price(i,j) = S0*d^((i-1)-(j-1))*u^(j-1);
    end    
end

% bereken de put optie prijs op tijdstip T

for j=1: 1: n
    option_price(n,j) = max(strike_price-stock_price(n,j),0);
end

% reken terug met de formule max(f_p,K-S_i)

for i = n-1: -1 : 1
    for j = 1: 1 : i
        option_price(i,j) = max(exp(-r*delta_t)*(...
            p*option_price(i+1,j+1)+...
            (1-p)*option_price(i+1,j)),strike_price-stock_price(i,j));
    end
end

Am_Put = option_price(1,1);
end