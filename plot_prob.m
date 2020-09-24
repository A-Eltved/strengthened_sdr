function [minVal,h] = plot_prob(prob)
% SYNTAX: h = plot_prob(prob) 
%
% DESCRIPTION:
% Plots 2D problems




if prob.data.n ~= 2
    disp('Can only plot for 2D')
    return
end
r = prob.data.r;
R = prob.data.R;
a = prob.data.a;
b = prob.data.b;
c = prob.data.c;

H = prob.data.H;
g = prob.data.g;

np = 1000;
xlin = linspace(-R-0.1,R+0.1,np);
[X1,X2] = meshgrid(xlin,xlin);
B = zeros(size(X1));
VAL = B;

for i = 1:np
    for j = 1:np
        xij = [X1(i,j);X2(i,j)];
        [B(i,j),in] = in_set(xij,r,R,a,b,c);
        if B(i,j) == 1
            VAL(i,j) = xij'*H*xij + 2*g'*xij;
%         else
%             B(i,j) = 1.3;
        end
    end
end

VAL(B==0) = max(max(VAL))+1;
[minVec,idy] = min(VAL);
[minVal,idx] = min(minVec);
idy = idy(idx);
xmin = X1(idy,idx);
ymin = X2(idy,idx);
fprintf('Minimum numerical value: %.4e, \t at (x1,x2)=(%.2f,%.2f)\n',minVal,xmin,ymin)


h = figure;
%title_str = sprintf('b = %.2f, a = (%.2f, %.2f)',b,a(1),a(2));
%title(title_str)

contourf(X1,X2,VAL,10);
xlabel('$x_1$','Interpreter','latex')
ylabel('$x_2$','Interpreter','latex')
axis equal
colorbar

hold on
plot(prob.data.xhat(1),prob.data.xhat(2),'k*')
plot(xmin,ymin,'gx')
hold off
legend('$x^T H x + 2 g^T x$','$\hat{x}$','Numerical minimum','Interpreter','latex')
end

function [bool,in] = in_set(x,r,R,a,b,c)
in = zeros(3,1);
in(1) = (norm(x,2) <= R);
in(2) = (norm(x,2) >= r);
in(3) = (norm(x - c,2) <= b'*x - a);
if sum(in) == 3
    bool = 1;
else
    bool = 0;
end
end
