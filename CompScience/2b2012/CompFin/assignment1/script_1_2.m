% Script assignment 1.2

% Perform n tests for different step sizes, calculate mean and std of debt
% at end (NOTE: THIS IS SLOW!!)
% 
% n = 1000;
% debt = linspace(0,0,n);
% stock_price = linspace(0,0,n);
% for steps=[104,365]%,104,156,365]
%     for j=1:n
%         j;
%         [S,~,~,~,final_debt] = hedge_sim(100,99,0.06,1./steps,0.20,0.20);
%         debt(j) = final_debt;
%         stock_price(j) = S(steps+1);
%     end
%     0.20
%     0.20
%     steps_year = steps
%     mean_debt = mean(debt)
%     std_debt = std(debt)
% end
% 
% n = 1000;
% debt = linspace(0,0,n);
% stock_price = linspace(0,0,n);
% for steps=[104,365]%,104,156,365]
%     for j=1:n
%         j;
%         [S,~,~,~,final_debt] = hedge_sim(100,99,0.06,1./steps,0.20,0.19);
%         debt(j) = final_debt;
%         stock_price(j) = S(steps+1);
%     end
%     0.20
%     0.19
%     steps_year = steps
%     mean_debt = mean(debt)
%     std_debt = std(debt)
% end
% 
n = 1000;
debt = linspace(0,0,n);
stock_price = linspace(0,0,n);
for steps=[365]%,104,156,365]
    for j=1:n
        j
        [S,~,~,~,final_debt] = hedge_sim(100,99,0.06,1./steps,0.20,0.21);
        debt(j) = final_debt;
        stock_price(j) = S(steps+1);
    end
    0.20
    0.21
    steps_year = steps
    mean_debt = mean(debt)
    std_debt = std(debt)
end



%Plot normal distribution
% steps = 365;
% n = 10000
% x = linspace(3.8,5.4,25)
% stock_price = linspace(0,0,n)
% for j=1:n
%     j
%     S = hedge_sim(100,99,0.06,1./steps,0.20,0.20);
%     stock_price(j) = S(steps+1);
% end
% h2=figure;
% hist(log(stock_price),x)
% hold
% mu = (0.06-0.5*0.20^2)+log(100)
% sd = 0.20
% ix = -3*sd+mu:1e-3:3*sd+mu; %covers more than 99% of the curve
% iy = 666.66*pdf('normal', ix, mu, sd);
% plot(ix,iy)
% 
% saveTightFigure(h2,'dist.pdf')


% Plot stock price, call price, delta and loan per week (as example) %
steps =365;
[S,call,delta,loan,final_debt] = hedge_sim(100,99,0.06,1./steps,0.20,0.20);
h1=figure;
subplot(2,2,1)
plot(0:steps,S)
xlabel('Day')
ylabel('Stock price')

subplot(2,2,2)
plot(0:steps,call)
xlabel('Day')
ylabel('Call price')

subplot(2,2,3)
plot(0:steps-1,delta)
xlabel('Day')
ylabel('Delta')

subplot(2,2,4)
plot(0:steps-1,loan)
xlabel('Day')
ylabel('Debt')
