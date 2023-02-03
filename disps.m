function displacements = disps(positions_temp,enzyme_locations,num_enzymes);
    for i=1:num_enzymes
        displacements(i,1)=((positions_temp(i,1)-enzyme_locations(i,1))^2+(positions_temp(i,2)-enzyme_locations(i,2))^2+(positions_temp(i,3)-enzyme_locations(i,3))^2)^0.5;
    end;
end

