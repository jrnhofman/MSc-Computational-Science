%Script Assignment 1 Part 1%

%Part A%
for i=1:5
    EurCall = EuCall(50,100,0.05*i,99,0.06);
    EurPut = EuPut(50,100,0.05*i,99,0.06);
end
EurCall = EuCall(50,100,0.50,99,0.06);
EurPut = EuPut(50,100,0.50,99,0.06);

%Part B%
n = 1000;
A=zeros(1,(n-10)/3);
B=zeros(1,n);
C=zeros(1,n);
D=zeros(1,n);
bcall=zeros(1,n);
bput=zeros(1,n);
j=1
for i=10:3:n
    i
    A(j)=EuCall(i,100,0.20,99,0.06);
    %B(i)=EuPut(i,100,0.20,99,0.06);
    %C(i)=AmCall(i,100,0.20,99,0.06);
    %D(i)=AmPut(i,100,0.20,99,0.06);
    %[bcall(i),bput(i)]=blsprice(100,99,0.06,1,0.20,0);
    j = j+1;
end
h=figure;
plot(10:3:n,A)
%plot(1:n,A,'r',1:n,bcall,'--r',1:n,B,'b',1:n,bput,'--b',1:n,C,'black',1:n,D,'g')
xlabel('Number of levels')
ylabel('European Call Option Price')
axis([10 1000 11.5 11.6])
%leg1=legend('Eur. call','BS Eur. call','Eur. put','BS Eur. put','Am. call','Am. put');
saveTightFigure(h,'closeup.pdf');

% %Part C%
% diff = EuCall(100,100,0.20,99,0.06)-EuPut(100,100,0.20,99,0.06)-100+99*exp(-0.06)
% 
% %Part D%
% option_price = zeros(1,100);
% delta1 = zeros(1,100);
% delta2 = zeros(1,100);
% delta3 = zeros(1,100);
% delta4 = zeros(1,100);
% for i=2:100
%     [test,option_price] = EuCall(i,100,0.20,99,0.06);
%     delta1(i) = (option_price(2,2)-option_price(2,1))/(100*exp(0.20*sqrt(1./(i-1))) - 100*exp(-0.20*sqrt(1./(i-1))));
%     [test,option_price] = EuPut(i,100,0.20,99,0.06);
%     delta2(i) = (option_price(2,2)-option_price(2,1))/(100*exp(0.20*sqrt(1./(i-1))) - 100*exp(-0.20*sqrt(1./(i-1))));
%     [test,option_price] = AmCall(i,100,0.20,99,0.06);
%     delta3(i) = (option_price(2,2)-option_price(2,1))/(100*exp(0.20*sqrt(1./(i-1))) - 100*exp(-0.20*sqrt(1./(i-1))));
%     [test,option_price] = AmPut(i,100,0.20,99,0.06);  
%     delta4(i) = (option_price(2,2)-option_price(2,1))/(100*exp(0.20*sqrt(1./(i-1))) - 100*exp(-0.20*sqrt(1./(i-1))));
% end
% h1=figure;
% plot(1:100,delta1,'r',1:100,delta2,'b',1:100,delta3,'black',1:100,delta4,'g')
% xlabel('Number of levels')
% ylabel('delta at first level')
% leg2=legend('Eur. call','Eur. put','Am. call','Am. put');
% saveTightFigure(h1,'delta.pdf');
%     
% %Part E%
% for i=1:5
%     AmeCall = AmCall(50,100,0.05*i,99,0.06)
%     AmePut = AmPut(50,100,0.05*i,99,0.06)
% end
% AmeCall = AmCall(50,100,0.50,99,0.06)
% AmePut = AmPut(50,100,0.50,99,0.06)