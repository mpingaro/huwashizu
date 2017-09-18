% Written by Marco Pingaro and Paolo Venini
%
% Solution present in the paper: Three field formulation

function [er_exx,er_eyy,er_exy] = error_strain_l2_norm(strain,cr,el,l)

    [gauss_w, gauss_p] = GaussQuad2D(4,4);
    npg = size(gauss_w,1);
    nel = size(el,1);
    
    strain = reshape(strain',3,[]);  
    str_xx = strain(1,:);
    str_yy = strain(2,:);
    str_xy = strain(3,:);

    err_uxx = 0.;
    err_uyy = 0.;
    err_uxy = 0.;
    norm_uxx = 0.;
    norm_uyy = 0.;
    norm_uxy = 0.;
    
     for i=1:nel
        pt = [cr(el(i,1),1),cr(el(i,1),2);
              cr(el(i,2),1), cr(el(i,2),2);
              cr(el(i,3),1), cr(el(i,3),2);
              cr(el(i,4),1), cr(el(i,4),2)];
        [DFF,JF] = jacobian(pt);        
        for j=1:npg
            [x,y] = map_quad(pt,gauss_p(j,:));   
            A = 2/(1+l);
            b = 1/25;
        
            Bx = 0.5*A*pi*cos(pi*x)*sin(pi*y);
            By = 0.5*A*pi*sin(pi*x)*cos(pi*y);
    
            u_xx = b*(-2*pi*sin(2*pi*y)*sin(2*pi*x) + Bx);
            u_xy = b*(2*pi*cos(2*pi*y)*(-1+cos(2*pi*x)) + By); 
            %u_yx = b*(2*pi*cos(2*pi*x)*(1-cos(2*pi*y)) + Bx);
            u_yy = b*(2*pi*cos(2*pi*x)*sin(2*pi*y) + By);    
            
            xi  = gauss_p(j,1);
            eta = gauss_p(j,2);
            
            psi(1) = 0.25*(1-xi)*(1-eta);
            psi(2) = 0.25*(1+xi)*(1-eta);
            psi(3) = 0.25*(1+xi)*(1+eta);
            psi(4) = 0.25*(1-xi)*(1+eta);
    
            sp_xx = str_xx(el(i,1))*psi(1) + str_xx(el(i,2))*psi(2) + ...
                str_xx(el(i,3))*psi(3) + str_xx(el(i,4))*psi(4);
            
            sp_yy = str_yy(el(i,1))*psi(1) + str_yy(el(i,2))*psi(2) + ...
                str_yy(el(i,3))*psi(3) + str_yy(el(i,4))*psi(4);
            
            sp_xy = str_xy(el(i,1))*psi(1) + str_xy(el(i,2))*psi(2) + ...
                str_xy(el(i,3))*psi(3) + str_xy(el(i,4))*psi(4);           
              
            err_uxx  = err_uxx + ((sp_xx - u_xx)^2)*gauss_w(j)*JF(j);
            err_uyy  = err_uyy + ((sp_yy - u_yy)^2)*gauss_w(j)*JF(j);
            err_uxy  = err_uxy + ((sp_xy - u_xy)^2)*gauss_w(j)*JF(j);
            
            norm_uxx = norm_uxx + (u_xx^2)*gauss_w(j);
            norm_uyy = norm_uyy + (u_yy^2)*gauss_w(j);
            norm_uxy = norm_uxy + (u_xy^2)*gauss_w(j);
        end
    end
    %
    norm_uxx = sqrt ( norm_uxx ) ;
    norm_uyy = sqrt ( norm_uyy ) ;
    norm_uxy = sqrt ( norm_uyy ) ;
    
    er_exx = sqrt ( err_uxx )/norm_uxx;
    er_eyy = sqrt ( err_uyy )/norm_uyy;
    er_exy = sqrt ( err_uxy )/norm_uxy;    

end           
    
















