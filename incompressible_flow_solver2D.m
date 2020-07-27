function [UX,UY,P,q] = incompressible_flow_solver2D(NX,NY,LX,LY,PW,PI,grid,dx,x,xi,dy,y,yi,lambdaw,lambdao);
%% Initialization

well = 1;
lambda = zeros(NY,NX);
lambdaharx = zeros(NY,NX+1);
lambdahary = zeros(NY+1,NX);
TX = zeros(NY,NX+1);
TY = zeros(NY+1,NX);

lambda = lambdaw + lambdao; % input lambda

[TX,TY,lambdaharx,lambdahary] = computetransmisibility2D(lambda,dx,NX,dy,NY);

%% Pressure Solver

A = zeros(NX*NY,NX*NY);
P = zeros(NX*NY,1);
q = zeros(NX*NY,1);
I = zeros(NY,NX);

for i = 1 : NX
    for j = 1:NY
            gridindex(j,i) = (j-1)*NX+i;
            [I,IN,IE,IS,IW] = findindex2D(j,i,NX);
        if i > 1 % there is left neighbor
            A(I,I)= A(I,I)+TX(j,i);
            A(I,IW)= -TX(j,i);
        end
        if i < NX % there is right neighbor
            A(I,I)= A(I,I)+TX(j,i+1);
            A(I,IE)= -TX(j,i+1);
        end
        if j > 1 % there is top neighbor
            A(I,I)= A(I,I)+TY(j,i);
            A(I,IN)= -TY(j,i);
        end
        if j < NY % there is bottom neighbor
            A(I,I)= A(I,I)+TY(j+1,i);
            A(I,IS)= -TY(j+1,i);
        end
    end
end

%% Insert BC
switch well
    case(0) % case no well
    for i = 1 : NX
        for j = 1 : NY
        [I,IN,IE,IS,IW] = findindex2D(j,i,NX);
        if i == 1; % left boundary
        A(I,I) = A(I,I)+ TX(j,i);
        q(I) = q(I)+TX(j,i)*PL;

        elseif i == NX; % right boundary
        A(I,I) = A(I,I)+ TX(j,i+1);
        q(I) = q(I)+TX(j,i+1)*PR;
        
        elseif j == 1
        A(I,I)= A(I,I)+TY(j,i);
        q(I) = q(I)+TY(j,i)*PT;
        
        elseif j == NY
        A(I,I)= A(I,I)+TY(j+1,i);
        q(I) = q(I)+TY(j+1,i)*PB;
        
        end
        end
    end
    
    case(1) % case with wells, constant PI and Pw
    [A,q] = wells2D(PI,PW,lambda,A,q,grid);
end

P = A\q;
P = P';
[UX,UY] = computevelocity2D(lambdaharx,lambdahary,dx,dy,P,NX,NY,i,j);

return
