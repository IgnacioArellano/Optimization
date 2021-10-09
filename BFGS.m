function dfp

%%%%función objetivo
U = @(x)  100*(sqrt(x(1).^2  + (x(2)+1).^2) -1).^2 +   90*(sqrt(x(1).^2  + (x(2)+1).^2) -1).^2 - (20.*x(1) + 40.*x(2));

%%%% Derivadas parciales para gradiente
Dx1= @(x,delta)   (U([x(1)+delta,x(2)])- U([x(1)-delta,x(2)]))/ (2*delta);
Dx2= @(x,delta)   (U([x(1),x(2)+delta])- U([x(1),x(2)-delta]))/ (2*delta);



%%%% Condiciones iniciales
delx=1e-3; % para derivada
eps=1e-10; % condicion de paro 1D
xi = [-1 1];
x = xi;
A=eye(length(x));
fprintf('Initial function value = %7.4f\n',U(x))
fprintf('No. x-vector f(x) Deriv \n')
fprintf('__________________________________________\n')
for i = 1:100
    if i==1
        grad_prev(1) = Dx1(x,delx);
        grad_prev(2) = Dx2(x,delx);
        search = -grad_prev;
        
        [alpha,fx_prev] = seccion_dorada(x,search ,xi,eps,U);
        if norm(grad_prev)<eps
            break;
        end
        x_c = x + alpha*search;
        fprintf('%3d %8.3f %8.3f % 8.3f %8.3f \n',i, x_c ,fx_prev,norm(grad_prev))
    else
        
        deltax = x_c-x;
        grad_c(1) = Dx1(x_c,delx);
        grad_c(2) = Dx2(x_c,delx);
        deltag = grad_c-grad_prev;
        term1 = (deltag'*deltag)/(deltag*deltax');
        term2 = (grad_prev'*grad_prev)/(grad_prev*search');
        A = A + term1 + term2;
        search = -inv(A)*grad_c';
        [alpha,fx_c] = seccion_dorada(x_c,search' ,xi,eps,U);
        
        if abs(fx_c-fx_prev)<eps || norm(grad_c)<eps
            break;
        end
         search= search';
        fx_prev = fx_c;
        grad_prev = grad_c;
        x=x_c;
        x_c = x+alpha*search;
        fprintf('%3d %8.3f %8.3f % 8.3f %8.3f \n',i,x_c,fx_c,norm(grad_c))
    end
end
fprintf('__________________________________________\n')
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end


function [alpha1,falpha1]=seccion_dorada (x,Si,lims,ep1,U)

tau = 0.381967;
alpha1 = lims(1)*(1-tau) + lims(2)*tau;
alpha2 = lims(1)*tau + lims(2)*(1-tau);
falpha1 = U(x+alpha1*Si);
falpha2 = U(x+alpha2*Si);

for i= 1:1000
    if falpha1 > falpha2
        lims(1) = alpha1;
        alpha1 = alpha2;
        falpha1 = falpha2;
        alpha2 = tau*lims(1) + (1-tau)*lims(2);
        falpha2 = U(x+alpha2*Si);
    else
        lims(2) = alpha2;
        alpha2 = alpha1;
        falpha2 = falpha1;
        alpha1 = tau*lims(2) + (1-tau)*lims(1);
        falpha1 = U(x+alpha1*Si);
    end
    if abs(U(x+alpha1*Si)- U(x+alpha2*Si)) < ep1
        break;
    end
end
end