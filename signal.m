function [modified,growing,largos]=signal(modified,growing,num_tubules,rescue,rescuemod,catas,catasmod,growth,growthmod,growth0,growthmod0,nucleated,shrink,shrinkmod,dt,discretization,dx,largos,CUTOFF)
chance = rand(num_tubules,1);
rescue0=rescue;
catas0=catas;
growth00=growth;
growth000=growth0;
shrink0=shrink;

for i=1:num_tubules
    %make a movement of the end
    signo=sign(growing(i)); 
    if size(modified{i},1)>=CUTOFF
        rescue=rescuemod;
        catas=catasmod;
        growth=growthmod;
        growth0=growthmod0;
        shrink=shrinkmod;
    else
        rescue=rescue0;
        catas=catas0;
        growth=growth00;
        growth0=growth000;
        shrink=shrink0;
    end


    if signo>0
        if ismember(growing(i,1),modified{i})
            time1=((growing(i,1))*dx(i)-largos(i))/growthmod;
            largos(i)=largos(i)+growthmod*min(dt,time1)+max(0,dt-time1)*growth;
        else
            dummy=ismember(growing(i,1),nucleated);
            largos(i)=largos(i)+growth*dt*(1-dummy)+growth0*dt*dummy;
        end
        growing(i)=ceil(largos(i)/dx(i));
        if largos(i)>dx(i)*discretization
            growing(i,1)=-discretization;
            largos(i)=dx(i)*discretization;
        else
            dummy=ismember(abs(growing(i,1)),modified{i});
            pshrink=1-exp(-(catasmod*dummy+catas*(1-dummy))*dt);
            if chance(i)<pshrink
                growing(i,1)=-growing(i,1);
            end              
        end
    end
    if signo<0
        templist=[];
        maxdisplacement=-max(shrink,shrinkmod)*dt;
        maxreach=largos(i)+maxdisplacement;
        maxdivisionreach=max(ceil(maxreach/dx(i)),1);
        templist=[abs(growing(i,1)):-1:maxdivisionreach];
        if isempty(templist)
            regdisplacement=-shrink*dt;
            largos(i)=largos(i)+regdisplacement;            
            growing(i,1)=-ceil(largos(i)/dx(i));
        else
            if isempty(modified{i})
                regdisplacement=-shrink*dt;
                largos(i)=largos(i)+regdisplacement;            
                growing(i,1)=-ceil(largos(i)/dx(i));                
            else
                P = zeros(1, max(max(templist),max(modified{i})) ) ;
                %templist=templist(templist>0);
                P(templist) = 1;
                inter = modified{i}(logical(P(modified{i})));                
                if isempty(inter)
                    regdisplacement=-shrink*dt;
                    largos(i)=largos(i)+regdisplacement;            
                    growing(i,1)=-ceil(largos(i)/dx(i));
                else
                    %this is when we try to delete a modified location
                    d_t=dt;st=1;countweird=0;
                    while d_t>0
					    if st<=size(templist,2)
						    dummy=ismember(templist(st),inter);
                            time1=(largos(i)-(templist(st)-1)*dx(i))/(shrinkmod*dummy+shrink*(1-dummy));
                        else
                            countweird=countweird+1;
						    dummy=0;
                            time1=(largos(i)-(templist(st-countweird)-1-countweird)*dx(i))/(shrinkmod*dummy+shrink*(1-dummy));
                        end
                        largos(i)=largos(i)-(shrinkmod*dummy+shrink*(1-dummy))*min(d_t,time1);
                        d_t=d_t-time1;
                        st=st+1;
                    end
                    growing(i,1)=-ceil(largos(i)/dx(i));
                end               
            end
        end
        if largos(i)<=0
            growing(i,1)=0;
            largos(i)=0;
        else
            dummy=ismember(abs(growing(i,1)),modified{i});
            pgrow=1-exp(-(rescuemod*dummy+rescue*(1-dummy))*dt);
            if chance(i)<pgrow
                  growing(i,1)=abs(growing(i,1));
            end
        end
        modified{i}=modified{i}(modified{i}<=abs(growing(i,1)));     
    end
end
















end
