function [a,b] = galerkin_method_flap(psi, p,q)
global gamma mu theta0 theta1c theta1s Omega lambda
%GALERKIN_METHOD_FLAP
%Original source: Professor Elisa Capello (elisa.capello@polito.it )

    a = zeros(5,5);
    b = zeros(5,1);

    xc = (gamma/8)*(1+(4/3)*mu*sin(psi));
    xk = 1+(gamma/8)*((4/3)*mu*cos(psi)+(mu^2)*sin(2*psi));
    L_theta = (gamma/8)*(1+(8/3)*mu*sin(psi)+2*(mu^2)*(sin(psi))^2);
    L_lambda = (-gamma/8)*((4/3)+2*mu*sin(psi));
    L_q = (gamma/8)*(cos(psi)+(2/3)*mu*sin(2*psi))-(2*sin(psi));
    L_p = (gamma/8)*sin(psi)*(1+(4/3)*mu*sin(psi))+(2*cos(psi));
    theta = theta0+theta1c*cos(psi)+theta1s*sin(psi);
    L_0 = (L_theta*theta+L_lambda*lambda+L_q*(q/Omega)+L_p*(p/Omega));
    
    a(1,1) = a(1,1)+xk;
    a(1,2) = a(1,2)+(-cos(psi)-xc*sin(psi)+xk*cos(psi));
    a(1,3) = a(1,3)+(-sin(psi)+xc*cos(psi)+xk*sin(psi));
    a(1,4) = a(1,4)+(-4*cos(2*psi)-2*xc*sin(2*psi)+xk*cos(2*psi));
    a(1,5) = a(1,5)+(-4*sin(2*psi)+2*xc*cos(2*psi)+xk*sin(2*psi));
    a(2,1) = a(2,1)+xk*cos(psi);
    a(2,2) = a(2,2)+(-cos(psi)-xc*sin(psi)+xk*cos(psi))*cos(psi);
    a(2,3) = a(2,3)+(-sin(psi)+xc*cos(psi)+xk*sin(psi))*cos(psi);
    a(2,4) = a(2,4)+(-4*cos(2*psi)-2*xc*sin(2*psi)+xk*cos(2*psi))*cos(psi);
    a(2,5) = a(2,5)+(-4*sin(2*psi)+2*xc*cos(2*psi)+xk*sin(2*psi))*cos(psi);
    a(3,1) = a(3,1)+xk*sin(psi);
    a(3,2) = a(3,2)+(-cos(psi)-xc*sin(psi)+xk*cos(psi))*sin(psi);
    a(3,3) = a(3,3)+(-sin(psi)+xc*cos(psi)+xk*sin(psi))*sin(psi);
    a(3,4) = a(3,4)+(-4*cos(2*psi)-2*xc*sin(2*psi)+xk*cos(2*psi))*sin(psi);
    a(3,5) = a(3,5)+(-4*sin(2*psi)+2*xc*cos(2*psi)+xk*sin(2*psi))*sin(psi);    
    a(4,1) = a(4,1)+xk*cos(2*psi);
    a(4,2) = a(4,2)+(-cos(psi)-xc*sin(psi)+xk*cos(psi))*cos(2*psi);
    a(4,3) = a(4,3)+(-sin(psi)+xc*cos(psi)+xk*sin(psi))*cos(2*psi);
    a(4,4) = a(4,4)+(-4*cos(2*psi)-2*xc*sin(2*psi)+xk*cos(2*psi))*cos(2*psi);
    a(4,5) = a(4,5)+(-4*sin(2*psi)+2*xc*cos(2*psi)+xk*sin(2*psi))*cos(2*psi);    
    a(5,1) = a(5,1)+xk*sin(2*psi);
    a(5,2) = a(5,2)+(-cos(psi)-xc*sin(psi)+xk*cos(psi))*sin(2*psi);
    a(5,3) = a(5,3)+(-sin(psi)+xc*cos(psi)+xk*sin(psi))*sin(2*psi);
    a(5,4) = a(5,4)+(-4*cos(2*psi)-2*xc*sin(2*psi)+xk*cos(2*psi))*sin(2*psi);
    a(5,5) = a(5,5)+(-4*sin(2*psi)+2*xc*cos(2*psi)+xk*sin(2*psi))*sin(2*psi); 
    
    b(1) = b(1)+L_0;
    b(2) = b(2)+L_0*cos(psi);
    b(3) = b(3)+L_0*sin(psi);
    b(4) = b(4)+L_0*cos(2*psi);
    b(5) = b(5)+L_0*sin(2*psi);

end