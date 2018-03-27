%% lab1
clc; close all; clear all;

%read the 3D object
ptCloud = pcread('hulk.ply');
%ptCloud = pcdownsample(ptCloud, 'gridAverage' , 1); % to downsample
showPointCloud(ptCloud);
axis([-5 ptCloud.XLimits(2) -5 ptCloud.YLimits(2) -5 ptCloud.ZLimits(2)]);
xlabel('X'); ylabel('Y'); zlabel('Z');
title('3D scene and projections onto two image planes');
hold on;

% take the 3dPoints
X = ptCloud.Location;

%% simulate sterovision system

%% camera 1
fx1 = 1 ; fy1 = 1; fc1 = 50; u01 = 0; v01 = 0; % pixel density and focal length
kx1 = fc1 * fx1 ; ky1 = fc1 * fy1;
K1 = [kx1 0 u01;
    0 ky1 v01;
    0 0 1]; % intrinsic parameters
normCam1 = [eye(3) zeros(3,1)]; %normalsised camera matrix
R1 = eye(3);
t1 = [ 0*fc1 0 0]'; % to put camera 1 at origin Extrinsic parmeters

%Camera Matric of camera1
P1 = K1 * normCam1 * [R1 t1; zeros(1,3) 1];

% to draw camera1 on plot
Rc1 = [1 0 0 ; 0 1 0; 0 0 1];
tc1 = [0 0 0]';
Draw3DCamera([Rc1 tc1], 'C1');

%% camera 2
fx2 = 1; fy2 = 1; fc2 = 50; u02 = 0; v02 = 0; % pixel density and focal length
kx2 = fc2 * fx2 ; ky2 = fc2 * fy2;
K2 = [kx2 0 u02;
    0 ky2 v02;
    0 0 1]; % intrinsic parameters
normCam2 = [eye(3) zeros(3,1)]; %normalsised camera matrix
R2 = [cos(0) -sin(0) 0;
    sin(0) cos(0) 0;
    0 0 1];
txx = 400; % trsnalatoin along x
t2 = [ -txx 0 0]'; %to trsnalate camera by 400 Extrinsic parameters

% the above '-' is becasue "THAT IS TRANSLATION OF WORLD COOROINATE ORIGIN
% WITHRESPECT TO CAMERA2"

%Camera Matric of camera2
P2 = K2 * normCam2 * [R2 t2; zeros(1,3) 1];

% to draw camera2 on plot
Rc1 = [cos(pi/4) sin(pi/4) 0 ; -sin(pi/4) cos(pi/4) 0; 0 0 1];
tc1 = [txx 0 0]';
Draw3DCamera([Rc1 tc1], 'C2', 'b');


%% Compute resulting images

%% camera1

% projections of scene onto camera planes
x1 = projection(P1, X');

C1 = camstruct('f', 50); % to construct camera structure to be used by imagept2plane function
% convert image pts to a plane
pt1 = imagept2plane(C1, x1(1:2,:), [0 0 fc1], [0 0 1]);

%define a plane at focal length fc1
transparency = 0.3;    %mostly clear
% XL = [10 200];
% YL = [10 200];
min1 = min(pt1(1:2,:),[],2);
max1 = max(pt1(1:2,:),[],2);
XL = [min1(1,1)-20 max1(1,1)+20];
YL = [min1(2,1)-20 max1(2,1)+20];
zpos = 50;
patch(XL([1 2 2 1 1]), YL([1 1 2 2 1]), zpos * ones(1,5), [0 .2 .7], 'FaceAlpha', transparency);


% to plot the points on image plane
plot3(pt1(1,:),pt1(2,:), pt1(3,:),'.','MarkerSize',10, 'LineWidth', 1,'Color','r');

% UNCOMMENT below to show the projection lines
% for i = 1:size(X,1)
%     plot3([X(i,1) 0],  [X(i,2) 0],[X(i,3) 0], 'LineWidth', 0.00002, 'Color', 'g');
% end

%% camera 2

%projections of scene onto camera planes
x2 = projection(P2, X');

C2 = camstruct('f', 50, 'P', [txx;0;0]); % to construct camera structure to be used by imagept2plane function
% convert image points to a  a plane
pt2 = imagept2plane(C2, x2(1:2,:), [0 0 fc1], [0 0 1]);

%define a plane at focal length 1
transparency = 0.3;    %mostly clear
% XL = [280 450];
% YL = [20 200];
min1 = min(pt2(1:2,:),[],2);
max1 = max(pt2(1:2,:),[],2);
XL = [min1(1,1)-20 max1(1,1)+20];
YL = [min1(2,1)-20 max1(2,1)+20];
zpos = 50;
patch(XL([1 2 2 1 1]), YL([1 1 2 2 1]), zpos * ones(1,5), [0 .2 .7], 'FaceAlpha', transparency);

% to plot the points on image plane
plot3(pt2(1,:),pt2(2,:), pt2(3,:),'.','MarkerSize',10, 'LineWidth', 1,'Color','b');

% UNCOMMENT below to show the projection lines
% for i = 1:size(X,1)
%     plot3([X(i,1) txx],  [X(i,2) 0],[X(i,3) 0], 'LineWidth', 0.00002, 'Color' , 'c');
% end

hold off;

%% to create images with the 2D Points

%create images
image1 = createImage(pt1(1:2,:));
image2 = createImage(pt2(1:2,:));

%% Essential matrix from known parameters (GROUND TRUTH)

tx = [-txx 0 0]; % trsnalation of world cordinate with respect to camera2
tx = star(tx); % skew symmetric matrix

R = [cos(0) -sin(0) 0;
    sin(0) cos(0) 0;
    0 0 1];

% Essential Matrix
E = tx * R;

disp('************************************************************');
disp('The computed Essential matrix is');
disp(E)
disp('************************************************************');

%% estimate R and t from Essential matrix

[U D V] = svd(E); % compute SVD of E

U =  det(U) * U ;
V = det(V) * V;

if det(U) > 0 && det(V) > 0
    
    % translation is propertional to this
    t = [U(1,3) U(2,3) U(3,3)]';
    
    % Rotation
    D = [0 1 0
        -1 0 0
        0 0 1];
    
    Ra =  U * D * V';
    
    Rb = U * D' * V';
    
    
    % Camera configuration 1
    P1_est = [eye(3) zeros(3,1)];
    
    %Camera configurations 2
    P2A_est = [Ra t];
    P2B_est = [Ra -t];
    P2C_est = [Rb t];
    P2D_est = [Rb -t];
    
    Ht = [eye(3) zeros(3,1); -2*V(1,3) -2*V(2,3) -2*V(3,3) -1]; % to twist the cameras
    
    % to check the cheirality constraint
    
    % triangulate one point
    Xone = f_intersection(P1_est, P2A_est, x1(1:2,1), x2(1:2,1));
    
    c1 = Xone(3) * Xone(4);
    temp = P2A_est * Xone ;
    c2 = temp(3) * Xone(4);
    
    if c1 > 0 && c2 > 0
        R_Est = Ra;
        t_Est = t;
    elseif c1 < 0 && c2 < 0
        R_Est = Ra;
        t_Est = -t;
    else
        temp1 = Ht * Xone;
        if (Xone(3)*temp1(4)) > 0
            R_Est = Rb;
            t_Est = t;
        else
            R_Est = Rb;
            t_Est = -t;
        end
    end
    
    %estimated essential matrix
    E_Est = star(t_Est) * R_Est;
    
    disp('************************************************************');
    disp('The estimated rotation is');
    disp(R_Est);
    disp('The estimated translation is');
    disp(t_Est)
    disp('The estimated essential matrix is');
    disp(E_Est);
    disp('************************************************************');
else
    disp('Deteriminant of U or V of svd(E) is not greater than zero');
end

%% Compute 3D scene

% camera matrix 1 at origin
P1_Est = K1 * normCam1 * [eye(3) zeros(3,1); zeros(1,3) 1];
% P1_Est = [eye(3) zeros(3,1)];

% compute camera matrix 2 with intrinsic parameters of camera 2 and
% estimated Rotataion and Translation matrices for extrinsic parameters

P2_Est = K2 * normCam1 * [R_Est t_Est; zeros(1,3) 1]; % true comibination of R and t
P2_Est1 = K2 * normCam1 * [Ra t; zeros(1,3) 1];
P2_Est2 = K2 * normCam1 * [Ra -t; zeros(1,3) 1];
P2_Est3 = K2 * normCam1 * [Rb t; zeros(1,3) 1];
P2_Est4 = K2 * normCam1 * [Rb -t; zeros(1,3) 1];

% compute 3D scene by triangulation
X_Est = f_intersection(P1_Est, P2_Est, x1(1:2,:), x2(1:2,:));
X_Est1 = f_intersection(P1_Est, P2_Est1, x1(1:2,:), x2(1:2,:));
X_Est2 = f_intersection(P1_Est, P2_Est2, x1(1:2,:), x2(1:2,:));
X_Est3 = f_intersection(P1_Est, P2_Est3, x1(1:2,:), x2(1:2,:));
X_Est4 = f_intersection(P1_Est, P2_Est4, x1(1:2,:), x2(1:2,:));

figure
plot3(X_Est(1,:)', X_Est(2,:)', X_Est(3,:)','.','MarkerSize',10, 'LineWidth', 1,'Color','b');
title('reconstruction from estimated true r and T from ground truth essential matrix');


figure
subplot(2,2,1)
plot3(X_Est1(1,:)', X_Est1(2,:)', X_Est1(3,:)','.','MarkerSize',10, 'LineWidth', 1,'Color','r');
% title('Reconstruction from estimated R and t');

subplot(2,2,2)
plot3(X_Est2(1,:)', X_Est2(2,:)', X_Est2(3,:)','.','MarkerSize',10, 'LineWidth', 1,'Color','r');

subplot(2,2,3)
plot3(X_Est3(1,:)', X_Est3(2,:)', X_Est3(3,:)','.','MarkerSize',10, 'LineWidth', 1,'Color','r');

subplot(2,2,4)
plot3(X_Est4(1,:)', X_Est4(2,:)', X_Est4(3,:)','.','MarkerSize',10, 'LineWidth', 1,'Color','r');

mtit('Four reconstructions from four camera posiibilities (ground truth essential matrix)');

%% Estimate the essential matrix from a set of 2D correspondences

E_Est_fromPoints = eight_point_algo(x1(1:2,:), x2(1:2,:), 10 ,10);

disp('************************************************************');
disp('The estimated essential matrix from point correspondances (Eight Point Algo) is');
disp(E_Est_fromPoints);
disp('************************************************************');


%% estimate R and t from Essential matrix (estimated from point correspondences , see above)

[U D V] = svd(E_Est_fromPoints); % compute SVD of E

U =  det(U) * U ;
V = det(V) * V;

if det(U) > 0 && det(V) > 0
    
    % translation is propertional to this
    t = [U(1,3) U(2,3) U(3,3)]';
    
    % Rotation
    D = [0 1 0
        -1 0 0
        0 0 1];
    
    Ra =  U * D * V';
    
    Rb = U * D' * V';
    
    
    % Camera configuration 1
    P1_est = [eye(3) zeros(3,1)];
    
    %Camera configurations 2
    P2A_est = [Ra t];
    P2B_est = [Ra -t];
    P2C_est = [Rb t];
    P2D_est = [Rb -t];
    
    Ht = [eye(3) zeros(3,1); -2*V(1,3) -2*V(2,3) -2*V(3,3) -1]; % to twist the cameras
    
    % to check the cheirality constraint
    
    % triangulate one point
    Xone = f_intersection(P1_est, P2A_est, x1(1:2,1), x2(1:2,1));
    
    c1 = Xone(3) * Xone(4);
    temp = P2A_est * Xone ;
    c2 = temp(3) * Xone(4);
    
    if c1 > 0 && c2 > 0
        R_Est = Ra;
        t_Est = t;
    elseif c1 < 0 && c2 < 0
        R_Est = Ra;
        t_Est = -t;
    else
        temp1 = Ht * Xone;
        if (Xone(3)*temp1(4)) > 0
            R_Est = Rb;
            t_Est = t;
        else
            R_Est = Rb;
            t_Est = -t;
        end
    end
    
    %estimated essential matrix
    E_Est = star(t_Est) * R_Est;
    
    disp('************************************************************');
    disp('The estimated rotation is');
    disp(R_Est);
    disp('The estimated translation is');
    disp(t_Est)
    disp('The estimated essential matrix is');
    disp(E_Est);
    disp('************************************************************');
else
    disp('Deteriminant of U or V of svd(E) is not greater than zero');
end

%% Compute 3D scene

% camera matrix 1 at origin
P1_Est = K1 * normCam1 * [eye(3) zeros(3,1); zeros(1,3) 1];
% P1_Est = [eye(3) zeros(3,1)];

% compute camera matrix 2 with intrinsic parameters of camera 2 and
% estimated Rotataion and Translation matrices for extrinsic parameters

P2_Est = K2 * normCam1 * [R_Est t_Est; zeros(1,3) 1]; % true comibination of R and t
P2_Est1 = K2 * normCam1 * [Ra t; zeros(1,3) 1];
P2_Est2 = K2 * normCam1 * [Ra -t; zeros(1,3) 1];
P2_Est3 = K2 * normCam1 * [Rb t; zeros(1,3) 1];
P2_Est4 = K2 * normCam1 * [Rb -t; zeros(1,3) 1];

% compute 3D scene by triangulation
X_Est = f_intersection(P1_Est, P2_Est, x1(1:2,:), x2(1:2,:));
X_Est1 = f_intersection(P1_Est, P2_Est1, x1(1:2,:), x2(1:2,:));
X_Est2 = f_intersection(P1_Est, P2_Est2, x1(1:2,:), x2(1:2,:));
X_Est3 = f_intersection(P1_Est, P2_Est3, x1(1:2,:), x2(1:2,:));
X_Est4 = f_intersection(P1_Est, P2_Est4, x1(1:2,:), x2(1:2,:));

figure
plot3(X_Est(1,:)', X_Est(2,:)', X_Est(3,:)','.','MarkerSize',10, 'LineWidth', 1,'Color','b');
title('reconstruction from estimated true r and T from estimated essential matrix');

figure
subplot(2,2,1)
plot3(X_Est1(1,:)', X_Est1(2,:)', X_Est1(3,:)','.','MarkerSize',10, 'LineWidth', 1,'Color','r');
% title('Reconstruction from estimated R and t');

subplot(2,2,2)
plot3(X_Est2(1,:)', X_Est2(2,:)', X_Est2(3,:)','.','MarkerSize',10, 'LineWidth', 1,'Color','r');

subplot(2,2,3)
plot3(X_Est3(1,:)', X_Est3(2,:)', X_Est3(3,:)','.','MarkerSize',10, 'LineWidth', 1,'Color','r');

subplot(2,2,4)
plot3(X_Est4(1,:)', X_Est4(2,:)', X_Est4(3,:)','.','MarkerSize',10, 'LineWidth', 1,'Color','r');

mtit('Four reconstructions from four camera posiibilities (estimated essntial matrix)');
