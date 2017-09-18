% Written by Marco Pingaro and Paolo Venini
%
% Solution present in the paper: Three field formulation

function er_u = error_beam_l2_norm(sp,mc,el,cr,l,E,nu,f)
    nel = size(el,1);

    [gauss_w, gauss_p] = GaussQuad2D(2,2);
    npg = size(gauss_w,1);
   
    err_ux = 0.;
    err_uy = 0.;
    norm_ux = 0.;
    norm_uy = 0.;
    
    for i=1:nel
        pt = [cr(el(i,1),1),cr(el(i,1),2);
              cr(el(i,2),1), cr(el(i,2),2);
              cr(el(i,3),1), cr(el(i,3),2);
              cr(el(i,4),1), cr(el(i,4),2)];
        [DFF,JF] = jacobian(pt);        
        for j=1:npg

            
            [x,y] = map_quad(pt,gauss_p(j,:));
            
            sol_x = 2*f*(1-nu^2)/(E*l)*x*(l/2 - y);
            sol_y = f/(E*l)*( x^2 + nu/(1-nu)*(y^2-l*y) );
                  
            xi  = gauss_p(j,1);
            eta = gauss_p(j,2);
            
            psi(1) = 0.25*(1-xi)*(1-eta);
            psi(2) = 0.25*(1+xi)*(1-eta);
            psi(3) = 0.25*(1+xi)*(1+eta);
            psi(4) = 0.25*(1-xi)*(1+eta);
            
            sp_x = sp(mc(i,1))*psi(1) + sp(mc(i,3))*psi(2) + ...
                sp(mc(i,5))*psi(3) + sp(mc(i,7))*psi(4);
            
            sp_y = sp(mc(i,2))*psi(1) + sp(mc(i,4))*psi(2) + ...
                sp(mc(i,6))*psi(3) + sp(mc(i,8))*psi(4);
             
            err_ux  = err_ux + ((sp_x - sol_x)^2)*gauss_w(j)*JF(j);
            err_uy  = err_uy + ((sp_y - sol_y)^2)*gauss_w(j)*JF(j);
            
            norm_ux = norm_ux + (sol_x^2)*gauss_w(j);
            norm_uy = norm_uy + (sol_y^2)*gauss_w(j);

        end
    end
    %
    norm_u = sqrt ( norm_ux + norm_uy ) ;
    er_u = sqrt ( err_ux + err_uy )/norm_u;
    %er_u = sqrt(err_ux);
    %er_u = sqrt( err_ux/norm_ux);
    
end