function cylinder_locations = tubule_setup(num_tubules, tubule_size, center_rad, cell_rad, xcutoff, weight)
    areatotal=4*3.1416*cell_rad^2;
    arearight=2*3.1416*cell_rad*(cell_rad-xcutoff);
    arealeft=areatotal-arearight;
    tubesatright=floor(num_tubules*weight*arearight/areatotal);
    tubesatleft=num_tubules-tubesatright;
    areaofcylinderbase=3.1416*tubule_size^2;
    maximumnumleft=floor(arealeft/areaofcylinderbase);
    maximumnumright=floor(arearight/areaofcylinderbase);
    
    if or(tubesatleft>=maximumnumleft,tubesatright>maximumnumright)
        quit
    end
    
    cylinder_locations(1,:)=onetubule_setup(center_rad, cell_rad, xcutoff, 1, tubesatleft);
    for i = 2:num_tubules
        modify = 1;
        while modify
            cylinder_locations(i,:)=onetubule_setup(center_rad, cell_rad, xcutoff, i, tubesatleft);
            loop=1;j=1;
            while loop
                cilinder1(:,1)=cylinder_locations(i,[1,2,3]);
                cilinder1(:,2)=cylinder_locations(i,[4,5,6]);
                cilinder2(:,1)=cylinder_locations(j,[1,2,3]);
                cilinder2(:,2)=cylinder_locations(j,[4,5,6]);  
                
                if (cilinder1(1,1)-cilinder2(1,2))^2+(cilinder1(2,1)-cilinder2(2,1))^2+(cilinder1(3,1)-cilinder2(3,1))^2<=2*tubule_size
                    loop=0;
                else
                    [Distance,P0,P1] = dcylinders(cilinder1, tubule_size, cilinder2, tubule_size);
                    distancia(j)=Distance;
                    distancia(distancia==0)=NaN;
                    if distancia(j)<=0
                        loop=0;
                    else
                        j=j+1;
                        if j==i
                            loop=0;
                            modify=0;
                        end
                    end
                    minimo=min(distancia,[],'all','omitnan');                     
                    if minimo>0
                        modify = 0;
                        clearvars distancia minimo
                    end                    
                end
            end

        end
        disp(i);
    end
end