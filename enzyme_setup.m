function enzyme_locations = enzyme_setup(num_enzymes, enzyme_size,  center_rad, cell_size)
    enzyme_locations = zeros(num_enzymes, 3);
    enzyme_locations(1,:)=enzyme_pos1(enzyme_size,center_rad, cell_size); 
    number=2;
    if num_enzymes>1
        %Check there's no superposition of enzymes and enzymes
        counter=1;
        while counter==1
              enzyme_locations(number,:)=enzyme_pos1(enzyme_size,center_rad, cell_size);
              for number2=1:number-1
                  d(number2)=(enzyme_locations(number,1)-enzyme_locations(number2,1))^2+(enzyme_locations(number,2)-enzyme_locations(number2,2))^2+(enzyme_locations(number,3)-enzyme_locations(number2,3))^2;
              end;
              minimo=min(d);
              if min(d)>(2*enzyme_size)^2
                  number=number+1;
                  if number==num_enzymes+1
                      counter=0;
                  end;
              end;
        end;
    end;
end