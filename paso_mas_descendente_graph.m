function paso_mas_descendente_graph

%%%%función objetivo
U = @(x)  100*(sqrt(x(1).^2  + (x(2)+1).^2) -1).^2 +   90*(sqrt(x(1).^2  + (x(2)+1).^2) -1).^2 - (20.*x(1) + 40.*x(2));
U2 = @(x1,x2)  100*(sqrt(x1.^2  + (x2+1).^2) -1).^2 +   90*(sqrt(x1.^2  + (x2+1).^2) -1).^2 - (20.*x1 + 40.*x2);
[x1,x2]=meshgrid(-1:0.1:1,-1:0.1:1);

f=U2(x1,x2);
figure('color',[1 1 1])
surf(x1,x2,f), shading interp
hold on
plot3(-1,1,U([-1,1]),'ko', 'MarkerSize',10,'markerfacecolor',[0 0 0])
plot3(0.5,0,U([0.5,0]),'ro', 'MarkerSize',10,'markerfacecolor',[1 0 0])

%%%% Derivadas parciales para gradiente
Dx1= @(x,delta)   (U([x(1)+delta,x(2)])- U([x(1)-delta,x(2)]))/ (2*delta);
Dx2= @(x,delta)   (U([x(1),x(2)+delta])- U([x(1),x(2)-delta]))/ (2*delta);

%%%% Condiciones iniciales
delta=1e-3; % para derivada
ep1=1e-5; % condicion de paro 1D
ep2=1e-5; % condicion de paro 2D
xi = [-1 1];
x = xi;
fx_prev=U(x);

for i = 1:3000
    grad(1) = Dx1(x,delta);
    grad(2) = Dx2(x,delta);
    Si = -grad;
    [alpha,falpha] = seccion_dorada(x,Si,xi,ep1,U);
    for j= 0:alpha/20:alpha
            r=U(x+j*Si);
            xt=x+j*Si;
            plot3(xt(1),xt(2),r,'k.', 'MarkerSize',10,'markerfacecolor',[0 0 0])
    end
        
    if abs(falpha-fx_prev)<ep2 || norm(grad)<ep2
        break;
    end
    fx_prev = falpha;
    x = x + alpha*Si;
    fprintf('%3d %8.3f %8.3f % 8.3f %8.3f \n',i,x,falpha,norm(grad))
end

fprintf('__________________________________________\n')
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
