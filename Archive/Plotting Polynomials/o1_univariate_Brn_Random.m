function [] = o2_univariate()

% Given a set of two dimensional control points b_{i}, plot the points, 
% and plot the associated curve. This code is for non-evenly distributed
% control points.

b = [...
    1       1;
    2       7;
    8       6;
    12      2;
    ]

b = [
    0       0.05;
    0.5     -0.25   ;
    1       0.45;
]

% plot the control points
figure(1)
hold on
scatter(b(:,1),b(:,2),...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor',[0 .75 .75])



% plot the associated polynomial by taking values over t \in [0,1]
data = [];
t = deCasteljau(b,0.25)

for i = 0:0.01:2
    
    j = deCasteljau(b,i);
    data = [data ;  j];
end

figure(1)
hold on 

plot(data(:,1),data(:,2),'-')
hold off

end


function [y] = deCasteljau(b,t)


while size(b,1) > 1
    rows = size(b,1);
    
    new_b = zeros(rows-1,2);
    for i = 1:1:rows-1
        new_b(i,:) = (1-t).*b(i,:) + t.*b(i+1,:);
    end
    
    b = new_b;
    
end

y = b;
end