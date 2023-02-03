function plot_tubule = plot_tubule(center_rad, cell_rad,cylinder_locations,tubule_size,discretization,modified,growing)

for i=1:size(cylinder_locations,1)
    plot_tubule_color(cylinder_locations(i,:),tubule_size,modified{i},growing(i),discretization);
    hold on
end;
view(3)
end