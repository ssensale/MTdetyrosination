function [newposition,check2D,counterE] = brownian2d_mobile(counterE,vectors_cylinder,tubule_size, i, enzyme_location,zetaANGLE,zetaHEIGHT,kT,dt,center_rad,enzyme_size,cell_rad,cylinder_location,unbindingprobability,tubuleid,discretization,growing,dx,zeta3D,rand_values1)
    %Get point in the cylinder P0    
	counterE=counterE;
    cilinder1(:,1) = cylinder_location([1,2,3]);
    cilinder1(:,2) = cylinder_location([4,5,6]);
    [~,pointincylinder,~] = dcylindersphere(cilinder1, tubule_size, enzyme_location', enzyme_size);
    %Vector from position to beginning of cylinder
    Ap0=pointincylinder-cilinder1(:,1);
    %If tubule reduced and it's out now, kick it out
    t0 = (Ap0(1)*vectors_cylinder(1)+Ap0(2)*vectors_cylinder(2)+Ap0(3)*vectors_cylinder(3));
    section0=floor(discretization*t0/vectors_cylinder(4))+1; 
    if section0>=discretization
        section0=discretization;
    else if section0<=1
            section0=1;
        end;
    end;
    if section0>abs(growing)
        check2D=0;
        newposition = tubule_eject_2d(enzyme_location, pointincylinder, tubule_size, enzyme_size, cilinder1(:, 1),cilinder1(:, 2));    
    else
        %move

        %Use rotational Brownian motion calculation to get the angular movement
        theta = ((2* zetaANGLE * kT * dt)^(0.5)) * (1/zetaANGLE) *3 * pi * rand_values1(1);
        point2 = rotatepoint(pointincylinder, theta, vectors_cylinder([1:3]),cilinder1(:,1)')';

        %Move point an amount vector*number with number a random number from
        %Brownian motion in 1D. If it escapes the boundaries of the segment,
        %rebound it as usual
        distance_change = ((2* zetaHEIGHT * kT * dt)^(0.5)) * (1/zetaHEIGHT) * rand_values1(2);
        vector2 = distance_change .* vectors_cylinder([1:3]);

        newposition = point2 + vector2';

        %Vector from new_position to beginning of cylinder
        Ap1=newposition-cilinder1(:,1);
        %Vector from new_position to end of cylinder
        cilinder1(:,2)=cylinder_location([1,2,3])+vectors_cylinder([1:3])*dx*abs(growing);
        Ap2=newposition-cilinder1(:,2);
        %Their projections on the segment of the cylinder to impose rebound
        t1 = (Ap1(1)*vectors_cylinder(1)+Ap1(2)*vectors_cylinder(2)+Ap1(3)*vectors_cylinder(3));
        t2 = (Ap2(1)*vectors_cylinder(1)+Ap2(2)*vectors_cylinder(2)'+Ap2(3)*vectors_cylinder(3));

        if t1-enzyme_size<=0
          %  disp('rebound')
            amount_to_rebound=enzyme_size-t1;
            direction_to_rebound = vectors_cylinder([1:3])';
            newposition = newposition + 2*amount_to_rebound .* direction_to_rebound;  
        elseif t2+enzyme_size>=0
           % disp('outer rebound')
            amount_to_rebound=enzyme_size+t2;
            direction_to_rebound = -vectors_cylinder([1:3])';
            newposition = newposition + 2*amount_to_rebound .* direction_to_rebound;
        end; 

        differencebottom = newposition - cilinder1(:,1); 
        t1 = (differencebottom(1)*vectors_cylinder(1)+differencebottom(2)*vectors_cylinder(2)+differencebottom(3)*vectors_cylinder(3));

        section=floor(discretization*t1/vectors_cylinder(4))+1; 
        if section>=discretization
            section=discretization;
        else if section<=1
                section=1;
            end;
        end;
        unbinding = rand();
        
        if section>abs(growing)
            check2D=0;
            newposition = tubule_eject_2d(enzyme_location, pointincylinder, tubule_size, enzyme_size, cilinder1(:, 1),cylinder_location([4,5,6])');    
            counterE=0;
        else    
            if unbinding < unbindingprobability
               %fprintf('Enzyme %d Left Tubule %d\n',i,tubuleid)  
               rand_values = normrnd(zeros(1, 3), ones(1, 3));
               coefficient = ((2* zeta3D * kT * dt)^0.5) * (1/zeta3D);
               random = coefficient .* rand_values;
               displacement=(random(1)^2+random(2)^2+random(3)^2)^0.5;
               new_position = enzyme_location + random;       
               check2D = 0;
               newposition = tubule_eject_2d(enzyme_location, new_position, tubule_size, enzyme_size, cilinder1(:, 1), cylinder_location([4,5,6])');
               counterE=0;
            else
               check2D = tubuleid;
            end;
        end
    end
end

