clc;clear;reset(RandStream.getGlobalStream, sum(100*clock));
InitialConditions
while not(accepted)
    if not(nodynamics)
        TubulinFreeNumb=TubulinTotalNumb-1624*sum(largos);
        Tub=TubulinFreeNumb/(602*Vcell); 
        growth=max(0,kon*Tub);
        growthmod=growth;        
        [nucleated,growing,largos] = nucleate(growing, pnuc, num_tubules, dt,largos);
        [modified,growing,largos]=signal(modified,growing,num_tubules,rescue,rescuemod,catas,catasmod,growth,growthmod,growth0,growthmod0,nucleated,shrink,shrinkmod,dt,discretization,dx,largos,CUTOFF);
    end
    
    if mod(counter,stepsverlet1)==0
        [verletC,warningsC] = verletcyl(largos,num_enzymes,enzyme_locations,enzyme_size,rc,num_tubules,cylinder_locations,tubule_size,growing,vectors_cylinders,dx);
        positions_temp=enzyme_locations;
        if mod(counter,frequency2)==0
			counter2=counter2+1;
            for i=1:num_tubules
				tracker(counter2)=counter;
				lengthMTs(counter2,i)=largos(i);
				propdet(counter2,i)=size(modified{i},1)*dx(i);
                sections{i,1}=i;
                sections{i,2}=modified{i};    
            end
            
            if counter2>1
                delete(name1);delete(name2);delete(name3);
            end
            
            name5=strcat('Locations',num2str(counter),'.dat');
            locations=cat(2,(1:num_enzymes)', enzyme_locations);
            name1=strcat('Detyrosinated',num2str(counter),'.dat');
            name2=strcat('Lengths',num2str(counter),'.dat');
            name3=strcat('Tracktime',num2str(counter),'.dat');
            name4=strcat('Chunks',num2str(counter),'.dat');
            writematrix(propdet,name1,'Delimiter','tab')
            writematrix(lengthMTs,name2,'Delimiter','tab')
            writecell(sections,name4,'Delimiter','tab')
			for i=counter2:size(propdet,1)
            	TotalLength(i,1)=(counter2-1)*frequency2*dt;
				TotalLength(i,2)=sum(lengthMTs(i,:));
				TotalLength(i,3)=sum(propdet(i,:));
                TotalLength(i,4)=toc;
                TotalLength(i,5)=countofbindsF;
                TotalLength(i,6)=countofbindsM;
            end
            writematrix(TotalLength,name3,'Delimiter','tab')   
        end								 
    else
        displacements=disps(positions_temp,enzyme_locations,num_enzymes);
        [warningsC,positions_temp,verletC]=check(largos,positions_temp,displacements,warningsC,num_enzymes,enzyme_locations,enzyme_size,rc,num_tubules,cylinder_locations,tubule_size,growing,vectors_cylinders,dx,verletC);
    end

    rand_values = randn(3,num_enzymes)';
    for i=1:num_enzymes 
        if not(check2D(i))
            idx=verletC(:,2)==i;
            verletCk=zeros(sum(idx),3);
            contador=0;
            if sum(idx)>0
                for peloduro=1:size(idx,1)
                    if idx(peloduro)
                        contador=contador+1;
                        verletCk(contador,:)=verletC(peloduro,:);
                    end
                end
            else
                verletCk=[];
            end  
            [enzyme_locations(i,:),check2D(i)]=brownian3d(largos,num_tubules,tubule_size,i, enzyme_locations(i,:),zeta3D,kT,dt,center_rad,enzyme_size,cell_rad,verletC,cylinder_locations,binding_probability,discretization,vectors_cylinders,growing,verletCk,rand_values(i,:));
        else
			if counterenzymes(i)==0 %enzyme not bound
				chance = rand();
				if chance<still_prob
                    countofbindsF=countofbindsF+1;
					counterenzymes(i)=2;					
					[enzyme_locations(i,:),check2D(i),section,modified{check2D(i)},counterenzymes(i)]=brownian2d_fixed(counterenzymes(i),vectors_cylinders(check2D(i),:),tubule_size, i, enzyme_locations(i,:),zetaANGLE,zetaHEIGHT,kT,dt,center_rad,enzyme_size,cell_rad,cylinder_locations(check2D(i),:),unbinding_probability,check2D(i),discretization,growing(check2D(i)),modified{check2D(i)}, probabilitymodifiestubule,zeta3D,rand_values(i,1));
                else
                    countofbindsM=countofbindsM+1;
                    counterenzymes(i)=1;
					[enzyme_locations(i,:),check2D(i),counterenzymes(i)]=brownian2d_mobile(counterenzymes(i),vectors_cylinders(check2D(i),:),tubule_size, i, enzyme_locations(i,:),zetaANGLE,zetaHEIGHT,kT,dt,center_rad,enzyme_size,cell_rad,cylinder_locations(check2D(i),:),unbinding_probabilityslide,check2D(i),discretization,growing(check2D(i)),dx(check2D(i)),zeta3D,rand_values(i,:)); 
                end
			else
				if counterenzymes(i)==1 %enzyme sliding
					[enzyme_locations(i,:),check2D(i),counterenzymes(i)]=brownian2d_mobile(counterenzymes(i),vectors_cylinders(check2D(i),:),tubule_size, i, enzyme_locations(i,:),zetaANGLE,zetaHEIGHT,kT,dt,center_rad,enzyme_size,cell_rad,cylinder_locations(check2D(i),:),unbinding_probabilityslide,check2D(i),discretization,growing(check2D(i)),dx(check2D(i)),zeta3D,rand_values(i,:)); 
                else
					[enzyme_locations(i,:),check2D(i),section,modified{check2D(i)},counterenzymes(i)]=brownian2d_fixed(counterenzymes(i),vectors_cylinders(check2D(i),:),tubule_size, i, enzyme_locations(i,:),zetaANGLE,zetaHEIGHT,kT,dt,center_rad,enzyme_size,cell_rad,cylinder_locations(check2D(i),:),unbinding_probability,check2D(i),discretization,growing(check2D(i)),modified{check2D(i)}, probabilitymodifiestubule,zeta3D,rand_values(i,:));
				end
			end
        end
    end 






    counter=counter+1;
    if counter>totalsteps
        accepted=1;
        fprintf('###Maximum Simulation Time Reached###\n')
    end
end