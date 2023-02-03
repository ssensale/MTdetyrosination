function [newposition,check2D,section,modified,counterE] = brownian2d_fixed(counterE,vectors_cylinders,tubule_size, i, enzyme_locations,zetaANGLE,zetaHEIGHT,kT,dt,center_rad,enzyme_size,cell_rad,cylinder_locations,unbindingprobability,tubuleid,discretization,growing,modified, probabilitymodifiestubule,zeta3D,rand_values1)
    
    %Get point in the cylinder P0    
	counterE=counterE;
    cilinder1(:,1) = cylinder_locations([1,2,3]);
    cilinder1(:,2) = cylinder_locations([4,5,6]);
    AB=cilinder1(:,2)-cilinder1(:,1);
    [~,pointincylinder,~] = dcylindersphere(cilinder1, tubule_size, enzyme_locations', enzyme_size);
    
    %Use rotational Brownian motion calculation to get the angular movement
    theta = ((2* zetaANGLE * kT * dt)^(0.5)) * (1/zetaANGLE) *3*pi*rand_values1(1);%rand_values
    newposition = rotatepoint(pointincylinder, theta, vectors_cylinders([1:3]),cilinder1(:,1)')';

    %Evaluate section to determine where the enzyme is
    differencebottom = newposition - cilinder1(:,1); 
    t1 = (differencebottom(1)*vectors_cylinders(1)+differencebottom(2)*vectors_cylinders(2)+differencebottom(3)*vectors_cylinders(3));

    section=floor(discretization*t1/vectors_cylinders(4))+1; 
    if section>=discretization
        section=discretization;
    else if section<=1
            section=1;
        end;
    end;
    if section>abs(growing)
        check2D=0;
        newposition = tubule_eject_2d(enzyme_locations, newposition, tubule_size, enzyme_size, cilinder1(:, 1), cilinder1(:, 2));
		counterE=0;
        %fprintf('Enzyme %d Abandoned by Tubule %d\n',i,tubuleid)  
    else
        %fprintf('Enzyme %d Up by Tubule %d@ %d\n',i,tubuleid,section)  
        if not(ismember(section,modified))
            modification=rand();              
            if modification<probabilitymodifiestubule
                if isempty(modified)
                    modified=[section];
                else
                    modified=cat(1,modified,section); 
                end;
            end;            
        end;
        unbinding = rand();
        if unbinding < unbindingprobability
            %fprintf('Enzyme %d Left Tubule %d\n',i,tubuleid) 
            rand_values = normrnd(zeros(1, 3), ones(1, 3));
            coefficient = ((2* zeta3D * kT * dt)^0.5) * (1/zeta3D);
            random = coefficient .* rand_values;
            displacement=(random(1)^2+random(2)^2+random(3)^2)^0.5;
            new_position = enzyme_locations + random;
            check2D = 0;
            newposition = tubule_eject_2d(enzyme_locations, new_position, tubule_size, enzyme_size, cilinder1(:, 1), cilinder1(:, 2));
			counterE=0;
        else
            check2D = tubuleid;
        end
    end;
    
end

