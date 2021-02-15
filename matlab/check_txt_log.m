


M = dlmread('./../log/coputed_alpha.txt');


n2 = length(M(1,:));
n = (n2+4)/3;

points = M(1:2,:);
W = M(3,:);

path = M(1:2,1:n);

extra_points = M(1:2,n+1:end);

figure(1)
plot(path(1,:),path(2,:),'r.','MarkerSize',10)
axis equal
grid on
xlabel('x')
ylabel('y')

hold on
plot(extra_points(1,:),extra_points(2,:),'bo','MarkerSize',4)
hold off

%%


ws = [min(M(1,:)) max(M(1,:)) min(M(2,:)) max(M(2,:))] + [-1 1 -1 1]*2;

dxy = sqrt( (ws(2)-ws(1)) * (ws(4)-ws(3)) )/30;
xv = ws(1):dxy:ws(2);
yv = ws(3):dxy:ws(4);


%% Compute and plot alpha function

[X,Y] = meshgrid(xv,yv);
ALPHA = X*0;

for i = 1:1:length(xv)
    x = xv(i);
    for j = 1:1:length(yv)
        y = yv(j);
         
        ALPHA(j,i) = alpha_func(x,y,points,W);

    end
end


figure(1)

h = surf(X,Y,ALPHA);
h.LineStyle = 'none';
colorbar
hold on

plot(path(1,:),path(2,:),'r.','MarkerSize',25)
plot(extra_points(1,:),extra_points(2,:),'w.','MarkerSize',20)

hold off

axis equal
axis([ws min(ALPHA(:)) max(ALPHA(:))])

figure(1000)
C = contour(X,Y,ALPHA,[0 0]);
close(1000)
C(:,1) = [];
hold on
plot(C(1,:),C(2,:),'k','LineWidth',2)
hold off




%% Compute field

delta = 0.001;
Kf = 4;

Phi_X = [];
Phi_Y = [];

for i = 1:1:length(xv)
    x = xv(i);
    for j = 1:1:length(yv)
        y = yv(j);
        
        a = alpha_func(x,y,points,W);
        grad_a = ([alpha_func(x+delta,y,points,W); alpha_func(x,y+delta,points,W)] - [a; a])/delta;
        
        grad_a_H = [0 -1; 1 0]*grad_a;
        
        P = 0.5*a^2;
        
        G = (2/pi)*atan(Kf*a);
        H = sqrt(1-G^2);
        
        
        field = -G*grad_a/(norm(grad_a)+1e-6) + H*grad_a_H/(norm(grad_a)+1e-6);
        
        
        Phi_X(j,i) = field(1);
        Phi_Y(j,i) = field(2);
    end
end

figure(2)
quiver(X,Y,Phi_X,Phi_Y,'b')
hold on
plot(path(1,:),path(2,:),'r.','LineWidth',1.0,'MarkerSize',25)
plot(C(1,:),C(2,:),'k','LineWidth',2)
hold off
axis equal
axis(ws)
