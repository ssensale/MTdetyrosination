function new_position = tubule_rebound_3d(old_position, newposition, tubule_rad, enzyme_size, cylinder)
   vector = newposition - old_position;
   vector_mag = (vector(1) ^ 2 + vector(2) ^ 2 + vector(3) ^ 2) ^ 0.5;
   line = [old_position(1), old_position(2), old_position(3), vector(1), vector(2), vector(3)];
   cylinder = [cylinder(1) cylinder(2) cylinder(3) cylinder(4) cylinder(5) cylinder(6) tubule_rad+enzyme_size];
   intersection = intersectLineCylinder(line, cylinder);

   if ~isempty(intersection)
       if size(intersection, 1) > 1
           distance_1 = (newposition(1) - intersection(1, 1)) ^ 2 + ...
               (newposition(2) - intersection(1, 2)) ^ 2 + ...
               (newposition(3) - intersection(1, 3)) ^ 2;

           distance_2 = (newposition(1) - intersection(2, 1)) ^ 2 + ...
               (newposition(2) - intersection(2, 2)) ^ 2 + ...
               (newposition(3) - intersection(2, 3)) ^ 2;

           if distance_1 < distance_2
               intersection = intersection(1, :);
           else
               intersection = intersection(2, :);
           end 

       elseif size(intersection, 1) == 1
           %disp('Size issue');
           intersection = intersection(1, :);
       end

       distance_inside = ((newposition(1) - intersection(1)) ^ 2 + ...
           (newposition(2) - intersection(2)) ^ 2 + ...
           (newposition(3) - intersection(3)) ^ 2) ^ 0.5;

       change = -1 .* vector / vector_mag * distance_inside;
       new_position = old_position + change;
       
   else
       new_position = old_position;
   end
end