function SR_1

%%%%función objetivo
U = @(x)  100*(sqrt(x(1).^2  + (x(2)+1).^2) -1).^2 +   90*(sqrt(x(1).^2  + (x(2)+1).^2) -1).^2 - (20.*x(1) + 40.*x(2));
U2 = @(x1,x2)  100*(sqrt(x1.^2  + (x2+1).^2) -1).^2 +   90*(sqrt(x1.^2  + (x2+1).^2) -1).^2 - (20.*x1 + 40.*x2);
[x1,x2]=meshgrid(-1:0.1:1,-1:0.1:1);
f=U2(x1,x2);
% condicion de paro 2D
delta = 1e-3;
ep= 1e-3;
xl = [-1;1];
x = xl;
fx_prev=U(x);
%%%% Derivadas parciales para gradiente
Dx1= @(x,delta)   (U([x(1)+delta,x(2)])- U([x(1)-delta,x(2)]))/ (2*delta);
Dx2= @(x,delta)   (U([x(1),x(2)+delta])- U([x(1),x(2)-delta]))/ (2*delta);

%%%% Condiciones iniciales
% Gradiente
grad(1) = Dx1(x,delta);
grad(2) = Dx2(x,delta);
grad = [grad(1);grad(2)];
% Hessiana
H = eye(2);
% Seccion dorada
Si = -grad;
[alpha,~] = seccion_dorada(x,Si,xl,ep,U);
figure('color',[1 1 1]);
hold on
contour(x1,x2,f,50)
plot(0.5,0,'ro','MarkerFaceColor','r', 'MarkerSize',10)
plot(x(1),x(2),'go','MarkerFaceColor','g')
fprintf('Initial function value = %7.4f\n ',fx_prev)
fprintf(' No. x-vector f(x) norma \n')
fprintf('__________________________________________\n')
for i = 1:100
    xi = x - alpha*H*grad;
    %Gradiente
    grad2(1) = Dx1(xi,delta);
    grad2(2) = Dx2(xi,delta);
    grad2 = [grad2(1);grad2(2)];
    Pk = xi - x;
    qk = grad2 - grad;
    H = H + (Pk-H*qk)*transpose((Pk-H*qk))/(transpose(qk)*(Pk-H*qk));
    Si = -H*grad;
    [alpha,falpha] = seccion_dorada(xi,Si,xl,ep,U);
    plot(xi(1),xi(2),'bo','MarkerFaceColor','b')
    plot([x(1),xi(1)],[x(2),xi(2)],'b')
    if abs(falpha-fx_prev) < ep
       break;
    end
    x = xi;
    grad = grad2;
    fx_prev = falpha;
    fprintf('%3d %8.3f %8.3f % 8.3f %8.3f \n',i,x(1),x(2),U(x),norm(grad))
end
fprintf('__________________________________________\n')
end
function [alpha1,falpha1]=seccion_dorada(x,Si,lims,ep1,U)
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