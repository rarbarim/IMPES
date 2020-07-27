function [krw,kro,lambdaw,lambdao,dlambdaw,dlambdao] = relativepermeability2D(Sw,Swc,Sor,no,nw,k,kroe,krwe,visco_w,visco_o)
    k = k.*(Sw.^0);

    krw = krwe.*((Sw-Swc)/(1-Swc-Sor)).^nw; 
    kro = kroe.*((1-Sw-Sor)/(1-Swc-Sor)).^no;
    
    lambdaw = k.*krw./visco_w;
    lambdao = k.*kro./visco_o; 
    
    dlambdaw = nw*k.*(krwe./visco_w).*((Sw-Swc)./(1-Swc-Sor)).^(nw-1);
    dlambdao = -no*k.*(kroe./visco_o).*((1-Sw-Sor)./(1-Swc-Sor)).^(no-1);
return