% |---*---|---*---|---*---|
% Note : number of interface = number of cells + 1

close all
clear all
format short

%% Default input
LX = 100;
LY = 100;
NX = 30;
NY = 30;
time = 3e8;
rho = 1000; % water density in kg/m3
phi = 0.3;
k = 1e-12; % absolute permeability in m2
visco_w = 1e-3;
visco_o = 10e-3;
nw = 3.5;
no = 2;
Swc = 0.2;
Sor = 0.1;
krwe = 0.7;
kroe = 0.8;
CFL = 0.3;
grid = [1 NX*NY];
PI = [100 100];
PW = [10000 1000];
homogenous = 0; % 1 for homogenous reservoir, otherwise heteroegenous

Sw_dummy = linspace(Swc,1-Sor,NX*NY);

q = zeros(NY,NX);

%% Initialization

dx = LX/NX; %grid size
x = linspace(dx/2,LX-dx/2,NX); %location of grid center
xi = linspace(0,LX,NX+1); %location of interfaces

dy = LY/NY; %grid size
y = linspace(dy/2,LY-dy/2,NY); %location of grid center
yi = linspace(0,LY,NY+1); %location of interfaces

[xx yy] = meshgrid(x,y);

[~,~,lambdaw,lambdao,dlambdaw,dlambdao] = relativepermeability2D(Sw_dummy,Swc,Sor,no,nw,k,kroe,krwe,visco_w,visco_o);
[~,~,dfwds] = fractionalflow2D(lambdaw,lambdao,dlambdaw,dlambdao);
dt = 0.5*dx;

%% IMPES
Sw = ones(NY,NX)*Swc;
Sw(1,1) = 1-Sor;

if homogenous ~= 1
    k =  k.*(10.^rand(NY,NX));
end

A = zeros(NX*NY,NX*NY);

t = 0;
m = 1;

while t < time
[~,~,lambdaw,lambdao,dlambdaw,dlambdao] = relativepermeability2D(Sw,Swc,Sor,no,nw,k,kroe,krwe,visco_w,visco_o);
[fw,fo,~] = fractionalflow2D(lambdaw,lambdao,dlambdaw,dlambdao);
[UX,UY,P,qt] = incompressible_flow_solver2D(NX,NY,LX,LY,PW,PI,grid,dx,x,xi,dy,y,yi,lambdaw,lambdao);
q = fw.*vec2mat(qt,NX);

ux(:,:,m) = UX;
uy(:,:,m) = UY;
p(:,:,m) = P;
sw(:,:,m) = Sw;

Sw = saturation2D(Sw,dt,q,phi,dx,dy,UX,UY,fw,NX,NY); % find saturation explicitly

speed  = max(max(max(UX)),max(max(UY)))/phi*max(dfwds);

if max(max(UX)) > max(max(UY))
    dt = CFL*(dx)/speed;
else
    dt = CFL*(dy)/speed;
end

m = m + 1;
t = t + dt;
  
end

%% Saturation Plot
figure
surf(xx,yy,sw(:,:,end))
xlabel('x location')
ylabel('y location')
zlabel('Sw')
zlim([0 1])
title('2D IMPES Saturation Profile')

%% Velocity Plot
figure
surf(ux(:,:,end))
xlabel('x location')
ylabel('y location')
zlabel('UX')
title('2D IMPES Velocity-X Profile')

figure
surf(uy(:,:,end))
xlabel('x location')
ylabel('y location')
zlabel('UY')
title('2D IMPES Velocity-Y Profile')

%% Pressure Plot
figure
surf(xx,yy,vec2mat(p(:,:,end),NX))
xlabel('x location')
ylabel('y location')
zlabel('Pressure')
title('2D IMPES Pressure Profile')