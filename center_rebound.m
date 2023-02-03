function new_position = center_rebound(old_position, newposition, center_rad, enzyme_size)
   vector = newposition - old_position;
   vector_mag = (vector(1) ^ 2 + vector(2) ^ 2 + vector(3) ^ 2) ^ 0.5;

   distance_inside = (center_rad + enzyme_size) - (newposition(1) ^ 2 + newposition(2) ^ 2 + newposition(3) ^ 2) ^ 0.5;
   
   change = -1 .* vector / vector_mag * distance_inside;
   new_position = old_position + change;
   
   new_position_mag = (new_position(1) ^ 2 + new_position(2) ^ 2 + new_position(3) ^ 2) ^ 0.5;
   
   if new_position_mag < center_rad + enzyme_size
       old_position_norm = (old_position(1) ^ 2 + old_position(2) ^ 2 + old_position(3) ^ 2) ^ 0.5;
       vector_to_outside = old_position / old_position_norm;
       change = vector_to_outside * distance_inside;
       new_position = old_position + change;
   end
end