function devia = deviatore(tenso);

indice = [1 3 4 6 7 9];

devia = tenso;

tracce = 1/2*[tenso(1,:)+tenso(3,:)
    tenso(1,:)+tenso(3,:) 
    tenso(4,:)+tenso(6,:) 
    tenso(4,:)+tenso(6,:)
    tenso(7,:)+tenso(9,:) 
    tenso(7,:)+tenso(9,:)];

devia(indice,:) = devia(indice,:)-tracce;

end
