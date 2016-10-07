% Written by Marco Pingaro and Paolo Venini
%
% Solution present in the paper: Three field formulation

function [er_exx,er_eyy,er_exy] = error_l2_norm(strain,cr,l,m)
    np = size(cr,1);

    strain = reshape(strain',3,[]);  
    str_xx = strain(1,:);
    str_yy = strain(2,:);
    str_xy = strain(3,:);

    sol_exx = zeros(1,np);
    sol_eyy = zeros(1,np);
    sol_exy = zeros(1,np);

    for i=1:np
        pt(1) = cr(i,1);
        pt(2) = cr(i,2);

        A = 2/(1+l);
        b = 1/25;
        
        Bx = 0.5*A*pi*cos(pi*pt(1))*sin(pi*pt(2));
        By = 0.5*A*pi*sin(pi*pt(1))*cos(pi*pt(2));
    
        u_xx = b*(-2*pi*sin(2*pi*y)*sin(2*pi*x) + Bx);
        u_xy = b*(2*pi*cos(2*pi*y)*(-1+cos(2*pi*x)) + By); 
        u_yx = b*(2*pi*cos(2*pi*x)*(1-cos(2*pi*y)) + Bx);
        u_yy = b*(2*pi*cos(2*pi*x)*sin(2*pi*y) * By);

        sol_exx(i) = u_xx;
        sol_eyy(i) = u_yy;
        sol_exy(i) = 0.5*(u_xy+u_yx);
    end
    e_xx = sum ( ( str_xx(1,:) - sol_exx(1,:) ).^2 )/np;
    e_yy = sum ( ( str_yy(1,:) - sol_eyy(1,:) ).^2 )/np;
    e_xy = sum ( ( str_xy(1,:) - sol_exy(1,:) ).^2 )/np;

    norm_exx = sqrt ( sum( sol_exx.^2) );
    norm_eyy = sqrt ( sum( sol_eyy.^2) );
    norm_exy = sqrt ( sum( sol_exy.^2) );
    
    er_exx = sqrt ( e_xx )/norm_exx;
    er_eyy = sqrt ( e_yy )/norm_eyy;
    er_exy = sqrt ( e_xy )/norm_exy;

end 
