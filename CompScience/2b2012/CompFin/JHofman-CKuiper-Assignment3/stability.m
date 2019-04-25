function [coeff] = stability(deltax,deltat,r,sigma)
    coeff =((1-r*deltat)-2*sigma^2*deltat/deltax^2)^2+4*((r-0.5*sigma^2)*deltat/(2*deltax))^2;


end

