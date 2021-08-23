function compare_DG_CG
%%%%función objetivo
U = @(x)  100*(sqrt(x(1).^2  + (x(2)+1).^2) -1).^2 +   90*(sqrt(x(1).^2  + (x(2)+1).^2) -1).^2 - (20.*x(1) + 40.*x(2));
U2 = @(x1,x2)  100*(sqrt(x1.^2  + (x2+1).^2) -1).^2 +   90*(sqrt(x1.^2  + (x2+1).^2) -1).^2 - (20.*x1 + 40.*x2);
%%%% Derivadas parciales para gradiente
Dx1= @(x,delta)   (U([x(1)+delta,x(2)])- U([x(1)-delta,x(2)]))/ (2*delta);
Dx2= @(x,delta)   (U([x(1),x(2)+delta])- U([x(1),x(2)-delta]))/ (2*delta);


[x1,x2]=meshgrid(-1:0.1:1,-1:0.1:1);
f=U2(x1,x2);

xi = [-1 1];
x = xi;
xp=xi;
xpg=xi;
ep1=1e-5; % condicion de paro 1D
delta=1e-3; % para derivada

figure('color',[1 1 1]);
hold on
contour(x1,x2,f,50)
plot(0.5,0,'ro','MarkerFaceColor','r', 'MarkerSize',10)
plot(xi(1),xi(2),'go','MarkerFaceColor','g')
for i = 1:50
   %%%%%%%%%%% Descent Gradient
   grad(1) = Dx1(x,delta);
   grad(2) = Dx2(x,delta); 
   Si = -grad;
   [alpha, ~] = seccion_dorada(x,Si,xi,ep1,U);
    xp=x;
    x = x + alpha*Si;
    plot(x(1),x(2),'go','MarkerFaceColor','g')
    plot([xp(1),x(1)],[xp(2),x(2)],'g')
    %%%%%%%%%%% Conjugate Gradient
     if i==1
         grad_prev=grad;
         search_prev=Si;
         xg=x;
     else
        grad2(1) = Dx1(xg,delta);
        grad2(2) = Dx2(xg,delta);
        search = -grad2 +((norm(grad2)^2)/(norm(grad_prev)^2))*search_prev;
        [alpha, ~] = seccion_dorada(xg,search,xi,ep1,U);
        grad_prev = grad2;
        search_prev = search;
        xpg=xg;
        xg = xg + alpha*search;
        plot(xg(1),xg(2),'bo','MarkerFaceColor','b')
        plot([xpg(1),xg(1)],[xpg(2),xg(2)],'b')
     end 
    
end

xg
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
