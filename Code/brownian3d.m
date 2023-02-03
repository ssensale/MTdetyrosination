function [new_position,check2D] = brownian3d(largos,num_tubules,tubule_size, i, enzyme_locations,zeta3D,kT,dt,center_rad,enzyme_size,cell_rad,verletC,cylinder_locations,binding_probability,discretization,vectors_cylinders,growing,verletCk,rand_values)
   
   coefficient = ((2* zeta3D * kT * dt)^0.5) * (1/zeta3D);
   random = coefficient .* rand_values;
   displacement=(random(1)^2+random(2)^2+random(3)^2)^0.5;
   new_position = enzyme_locations + random;

   check2D=0;
   if new_position(1)^2+new_position(2)^2+new_position(3)^2<(center_rad+enzyme_size)^2
      new_position=center_rebound(enzyme_locations, new_position, center_rad,enzyme_size);
   else
       if new_position(1)^2+new_position(2)^2+new_position(3)^2>(cell_rad-enzyme_size)^2
          new_position=cell_rebound(enzyme_locations, new_position, cell_rad,enzyme_size);
       end      
   end
   
  smallverlet=verletCk;
  
  
  if not(isempty(smallverlet))
       for j=1:size(smallverlet,1)
                if displacement>=smallverlet(j,3)%check if it moved enough to bind before checking if bind
                    POINTS = intersectLineCylinder([enzyme_locations,new_position-enzyme_locations], [cylinder_locations(smallverlet(j,1),[1,2,3]),cylinder_locations(j,[1,2,3])+vectors_cylinders(smallverlet(j,1),[1,2,3])*largos(smallverlet(j,1)),tubule_size], 'checkBounds', 1);
                    
                    if not(isempty(POINTS))
                        %%check if points are inside the segment
                        C1=POINTS(1,:);
                        A=enzyme_locations;
                        B=new_position;
                        AB=B-A;
                        AC1=C1-A;
                        distan1=(AC1(1)*AB(1)+AC1(2)*AB(2)+AC1(3)*AB(3));
                        distancia=AB(1)*AB(1)+AB(2)*AB(2)+AB(3)*AB(3);
                        if size(POINTS(:,1),1)>1
                            C2=POINTS(2,:);
                            AC2=C2-A;
                            distan2=(AC2(1)*AB(1)+AC2(2)*AB(2)+AC2(3)*AB(3));
                            if or(distan2>distancia,distan2*distancia<0)
                                POINTS(2,:)=[];
                                if or(distan1>distancia,distan1*distancia<0)
                                    POINTS(1,:)=[];
                                end
                            else
                                if or(distan1>distancia,distan1*distancia<0)
                                    POINTS(1,:)=[];
                                else
                                    if distan1>distan2
                                        X=POINTS(2,:);
                                        POINTS=X;
                                    else
                                        X=POINTS(1,:);
                                        POINTS=X;
                                    end
                                end
                            end
                            
                        else
                            if or(distan1>distancia,distan1*distancia<0)
                                POINTS=[];
                            end
                        end
                    end
                    if isempty(POINTS)
                        [distancias,closest,~]=dcylindersphere([cylinder_locations(smallverlet(j,1),[1,2,3])',(cylinder_locations(j,[1,2,3])+vectors_cylinders(smallverlet(j,1),[1,2,3])*largos(smallverlet(j,1)))'], tubule_size, new_position', enzyme_size);
                        differencebottom = closest - cylinder_locations(smallverlet(j,1),[1,2,3])'; 
                        AB=(cylinder_locations(j,[1,2,3])+vectors_cylinders(smallverlet(j,1),[1,2,3])*largos(smallverlet(j,1)))'-cylinder_locations(smallverlet(j,1),[1,2,3])';
                        AB_squared = AB(1)*AB(1)+AB(2)*AB(2)+AB(3)*AB(3);
                        NormalizedAB=AB/AB_squared^0.5;
                        t1 = (differencebottom(1)*NormalizedAB(1)+differencebottom(2)*NormalizedAB(2)+differencebottom(3)*NormalizedAB(3));
                        if distancias<=0
                            section=floor(discretization*t1/vectors_cylinders(smallverlet(j,1),4))+1; 
                            if section>=discretization
                                section=discretization;
                            end
                            if section==0
                                section=1;
                            end
                            %check if it sticks to a tubule
                            check2D=0;
                            chance = rand();
                            if and(chance <= binding_probability,section<=abs(growing(smallverlet(j,1))))
                                check2D = smallverlet(j,1);
                                %fprintf('Enzyme %d Bound to Tubule %d\n',i,j)  
                                else if section<=abs(growing(smallverlet(j,1)))
                                new_position = tubule_rebound_3d(enzyme_locations, new_position,...
                                tubule_size, enzyme_size, [cylinder_locations(smallverlet(j,1),[1,2,3])',(cylinder_locations(j,[1,2,3])+vectors_cylinders(smallverlet(j,1),[1,2,3])*largos(smallverlet(j,1)))']);
                                    end
                            end
                    %evaluate if it binds or rebounds
                        end
                    else
                        differencebottom = POINTS(1,:)' - cylinder_locations(smallverlet(j,1),[1,2,3])'; 
                        AB=(cylinder_locations(j,[1,2,3])+vectors_cylinders(smallverlet(j,1),[1,2,3])*largos(smallverlet(j,1)))'-cylinder_locations(smallverlet(j,1),[1,2,3])';
                        AB_squared = AB(1)*AB(1)+AB(2)*AB(2)+AB(3)*AB(3);
                        NormalizedAB=AB/AB_squared^0.5;
                        t1 = (differencebottom(1)*NormalizedAB(1)+differencebottom(2)*NormalizedAB(2)+differencebottom(3)*NormalizedAB(3));
                        section=floor(discretization*t1/vectors_cylinders(smallverlet(j,1),4))+1; 
                        if section>=discretization
                           section=discretization;
                        end
                        if section==0
                           section=1;
                        end
                        %check if it sticks to a tubule
                        check2D=0;
                        chance = rand();
                        if and(chance <= binding_probability,section<=abs(growing(smallverlet(j,1))))
                            check2D = smallverlet(j,1);
                            %fprintf('Enzyme %d Bound to Tubule %d\n',i,j)  
                            else if section<=abs(growing(smallverlet(j,1)))
                                    new_position = tubule_rebound_3d(POINTS(1,:), new_position,...
                                    tubule_size, enzyme_size, [cylinder_locations(smallverlet(j,1),[1,2,3])',(cylinder_locations(j,[1,2,3])+vectors_cylinders(smallverlet(j,1),[1,2,3])*largos(smallverlet(j,1)))']);
                                 end
                        end
                    end
                end
           end
       end
end

