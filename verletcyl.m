function [verlet,warningsC] = verletcyl(largos,num_enzymes,enzyme_locations,enzyme_size,rc,num_tubules,cylinder_locations,tubule_size,growing,vectors_cylinders,dx)
verlet=[0,0,0];
contar=0;
for i=1:num_tubules
    for j=1:num_enzymes
            cilinder1(:,1)=cylinder_locations(i,[1,2,3]);
            vector=vectors_cylinders(1,[1,2,3])';
            cilinder1(:,2)=cylinder_locations(i,[1,2,3])+vectors_cylinders(i,[1:3])*largos(i);
            sphere1=enzyme_locations(j,:)';
            distancias(i,j)=dcylindersphere(cilinder1, tubule_size, sphere1, enzyme_size);
            if distancias(i,j)<rc
                contar=contar+1;
                verlet(contar,1)=i;%tubule
                verlet(contar,2)=j;%enzyme
                verlet(contar,3)=distancias(i,j);%distance tubule-enzyme
            end;
    end;
end;
warningsC=distancias';


end

