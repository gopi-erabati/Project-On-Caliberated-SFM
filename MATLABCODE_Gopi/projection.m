function [ x ] = projection( P,X )
% this function to computer projections of 3d points onto 2d using camera
% parameters

X = [X;ones(1,size(X,2))];

x = P * X;

 x = x ./repmat(x(3,:),3,1) ;
 
%  x = x(1:2,:);
 
%  figure(1)
% plot(x(1,:),x(2,:),'o')
% figure(2)
% plot(X(1,:),X(2,:),'o')

end

