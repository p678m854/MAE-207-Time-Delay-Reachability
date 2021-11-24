%% Main Document
% Input document for the functions:

rho = 0; % rho < |0.2|
A = [-2 0; 0, -0.9+rho];
Ad = [-1, 0; -1, -1+0.5*rho];
B = [-0.5; 1];
h = 0.7;
utilde = 1;


[P, F, delta]=fridman_method1(A, Ad, B, h, utilde) 

%% Plot the Ellipse

[V,D] = eig(P);
u = V(:,1);
v = V(:,2);
l1 = D(1,1);
l2 = D(2,2);

pts = [];

delta = .1;

for alpha = -1/sqrt(l1)-delta:delta:1/sqrt(l1)+delta
    beta = sqrt((1 - alpha^2 * l1)/l2);
    pts(:,end+1) = alpha*u + beta*v;
end
for alpha = 1/sqrt(l1)+delta:-delta:-1/sqrt(l1)-delta
    beta = -sqrt((1 - alpha^2 * l1)/l2);
    pts(:,end+1) = alpha*u + beta*v;
end

plot(pts(1,:), pts(2,:))
xlabel('x1')
ylabel('x2')
title(['Bounding Ellipse for h = ',num2str(h)])
axis('equal')