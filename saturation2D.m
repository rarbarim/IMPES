function Sw = saturation2D(Sw,dt,q,phi,dx,dy,UX,UY,fw,NX,NY);

for i = 1 : NX
        for j = 1 : NY
            if i == 1 && j > 1
                Sw(j,i) = Sw(j,i) + dt*q(j,i)/phi - dt/(phi*dy)*(fw(j,i).*UY(j+1,i)-fw(j-1,i).*UY(j,i)) - dt/(phi*dx)*(fw(j,i).*UX(j,i+1)-0.*UX(j,i));
            elseif j == 1 && i > 1
                Sw(j,i) = Sw(j,i) + dt*q(j,i)/phi - dt/(phi*dy)*(fw(j,i).*UY(j+1,i)-0.*UY(j,i)) - dt/(phi*dx)*(fw(j,i).*UX(j,i+1)-fw(j,i-1).*UX(j,i));
            elseif i > 1 && j > 1
                Sw(j,i) = Sw(j,i) + dt*q(j,i)/phi - dt/(phi*dy)*(fw(j,i).*UY(j+1,i)-fw(j-1,i).*UY(j,i)) - dt/(phi*dx)*(fw(j,i).*UX(j,i+1)-fw(j,i-1).*UX(j,i));
            end
        end
end

return
    