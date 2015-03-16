function ene = satura(a,b)

an = [a(1,:)
    a(2,:)
    a(2,:)
    a(3,:)
    a(4,:)
    a(5,:)
    a(5,:)
    a(6,:)
    a(7,:)
    a(8,:)
    a(8,:)
    a(9,:)];

bn = [b(1,:)
    b(2,:)
    b(2,:)
    b(3,:)
    b(4,:)
    b(5,:)
    b(5,:)
    b(6,:)
    b(7,:)
    b(8,:)
    b(8,:)
    b(9,:)];

ene = sum(an.*bn);  

end