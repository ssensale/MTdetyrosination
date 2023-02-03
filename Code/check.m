function [warningsC,positions_temp,verletC]=check(largos,positions_temp,displacements,warningsC,num_enzymes,enzyme_locations,enzyme_size,rc,num_tubules,cylinder_locations,tubule_size,growing,vectors_cylinders,dx,verletC);
    sizeverlet=size(verletC,1);
    contar=0;
    for i=1:num_enzymes%enzyme
        for j=1:num_tubules%tubule
            if and(warningsC(i,j)-displacements(i)<rc,warningsC(i,j)>=rc)
                cilinder1(:,1)=cylinder_locations(j,[1,2,3]);
                vector=vectors_cylinders(j,[1,2,3]);
                cilinder1(:,2)=cylinder_locations(j,[1,2,3])+vector*largos(j);
                sphere1=enzyme_locations(i,:)';
                distancias(i,j)=dcylindersphere(cilinder1, tubule_size, sphere1, enzyme_size);
                positions_temp(i,:)=enzyme_locations(i,:);
                warningsC(i,j)=distancias(i,j);
                if distancias(i,j)<rc
                    contar=contar+1;
                    verletC(contar+sizeverlet,1)=j;%tubule
                    verletC(contar+sizeverlet,2)=i;%enzyme
                    verletC(contar+sizeverlet,3)=distancias(i,j);%distance tubule-enzyme
                end;
            end;
        end;
    end;
end

