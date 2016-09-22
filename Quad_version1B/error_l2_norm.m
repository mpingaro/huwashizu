% Written by Marco Pingaro and Paolo Venini
%

function er_u = error_l2_norm(sp,cr,l,m)
    np = size(cr,1);
    sp = reshape(sp, 2, [])';
    sol_a = zeros(np,2);
    sp = sp(1:np,:);

    for i=1:np
        pt(1) = cr(i,1);
        pt(2) = cr(i,2);

        A = 2/(1+l);
        B = 0.5*A*sin(pi*pt(1))*sin(pi*pt(2));
        b = 1/25;

        sol_a(i,1) = b*(sin(2*pi*pt(2))*(-1+cos(2*pi*pt(1))) + B);
        sol_a(i,2) = b*(sin(2*pi*pt(1))*( 1-cos(2*pi*pt(2))) + B);

        %sol_a(i,1) = cos( pi*pt(1) )*sin( 2*pi*pt(2) );
        %sol_a(i,2) = sin( pi*pt(1) )*cos( pi*pt(2) );
        
        %sol_a(i,1) = sin( pi*pt(1) )*sin( pi*pt(2) );
        %sol_a(i,2) = sin( pi*pt(1) )*sin( pi*pt(2) );

    end
    %er_ux = sqrt( sum( (sol_a(:,1) - sp(:,1)).^2 ) );
    %er_uy = sqrt( sum( (sol_a(:,2) - sp(:,2)).^2 ) );

    er_ux = sum ( ( sp(:,1) - sol_a(:,1) ).^2 );
    er_uy = sum ( ( sp(:,2) - sol_a(:,2) ).^2 );

    norm_u = sqrt ( sum( sol_a(:,1).^2 + sol_a(:,2).^2 ) );
    %
    er_u = sqrt ( er_ux + er_uy ); 
end 
