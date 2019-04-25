function price = CN(S0,K,r,sigma,deltat,deltax)
    
    Smin = 0.1*K;
    Smax = 10*K;
    M1 = log(Smin);
    M2 = log(Smax);
            
    t = round((log(S0)-M1)/deltax);
    M1 = log(S0) - t*deltax;
    N = round((M2-M1)/deltax - 1);
    M2 = M1 + (N+1)*deltax;
    n = 1/deltat;

    a1 = -(r-0.5*sigma*sigma)*deltat/(4*deltax) - 0.25*sigma*sigma*deltat/(deltax*deltax)
    a0 = 1 + 0.5*sigma*sigma*deltat/(deltax*deltax) + r*deltat/2
    am1 = (r-0.5*sigma*sigma)*deltat/(4*deltax) - 0.25*sigma*sigma*deltat/(deltax*deltax)
    
    b1 = (r-0.5*sigma*sigma)*deltat/(4*deltax) + 0.25*sigma*sigma*deltat/(deltax*deltax)
    b0 = 1 - 0.5*sigma*sigma*deltat/(deltax*deltax) - 0.5*r*deltat
    bm1 = -(r-0.5*sigma*sigma)*deltat/(4*deltax) + 0.25*sigma*sigma*deltat/(deltax*deltax)
    
    V = zeros(n+1,N+2);
    for i=1:n+1
        V(i,1) = 0;
        V(i,N+2) = blsprice(exp(M2),K,r,(i-1)*deltat,sigma);
        %V(i,N+2) = exp(M2)-K;
        %V(i,N+2) = (exp(M2)-K)*exp((i-1)*deltat);
        %V(i,N+2) = exp(r*(i-1)*deltat)*exp(M2)-K
    end
    for j=2:N+1
        V(1,j) = max((exp(M1+(j-1)*deltax)-K),0);
    end
    
    A = zeros(N);
    A(1,1) = a0;
    A(1,2) = a1;
    A(N,N-1) = am1;
    A(N,N) = a0;
    for i=2:N-1
        A(i,i-1) = am1;
        A(i,i) = a0;
        A(i,i+1) = a1;
    end
    
    c = linspace(0,0,N);

    for j=2:n+1
        c(1) = b1*V(j-1,3) + b0*V(j-1,2) + bm1*V(j-1,1) - am1*V(j,1);
        c(N) = b1*V(j-1,N+2) + b0*V(j-1,N+1) + bm1*V(j-1,N) - a1*V(j,N+2);
        for i=2:N-1
            c(i) = b1*V(j-1,i+2) + b0*V(j-1,i+1) + bm1*V(j-1,i);
        end
        y = A\c';
        for i=2:N+1
            V(j,i) = y(i-1); 
        end
    end
   
    price = V(n+1,t+1);
    blackscholes = blsprice(S0,K,r,1,sigma);
    diff = (price-blackscholes)/blackscholes * 100
    
    
%     price = V(n+1,:);
%     a = exp(linspace(M1,M2,N+2));
%     exact = blsprice(a,110,0.04,1,0.30);
%     
%     figure;
%     plot(a,price,a,exact,'-r','LineWidth',2)
%     xlabel('S_0')
%     ylabel('V(S_0)')
%     
%     figure;
%     subplot(1,2,1);
%     plot(a,abs(price-exact),'-r','LineWidth',2)
%     xlabel('S_0')
%     ylabel('Absolute error')
%     subplot(1,2,2);
%     plot(a,abs(price-exact)./exact,'-r','LineWidth',2)
%     xlabel('S_0')
%     ylabel('Relative error')
%     axis([0 450 0 0.05])
%     
%     diff = norm(price-exact)
%     
%     delta = linspace(0,0,N);
%     for i=1:N
%         delta(i) = (price(i+2) - price(i))/(2*deltax) / a(i+1);
%     end
%     
%     figure;
%     plot(a(2:end-1),delta,'-r','LineWidth',2)
%     xlabel('S_0')
%     ylabel('\Delta')
%     

end

