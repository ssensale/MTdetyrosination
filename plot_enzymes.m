function plot_enzymes = plot_enzymes(num_enzymes,enzyme_locations,enzyme_size)
       hold on
       for i=1:num_enzymes
           [X1,Y1,Z1]=sphere(5);
            XA1=X1*enzyme_size+enzyme_locations(i,1);
            YA1=Y1*enzyme_size+enzyme_locations(i,2);
            ZA1=Z1*enzyme_size+enzyme_locations(i,3);
            surface(XA1,YA1,ZA1,'FaceAlpha',1)
            hold on
       end;
end
