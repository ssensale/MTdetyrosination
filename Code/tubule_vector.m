function vector2=tubule_vector(cylinder_locations,num_tubules)
    for i=1:num_tubules
        cilinder1(i,:)=cylinder_locations(i,[1,2,3]);
        cilinder2(i,:)=cylinder_locations(i,[4,5,6]);
        %Get vector of cylinder
        vector(i,:) = cilinder2(i,:)-cilinder1(i,:);
        vector_mag(i) = (vector(i,1) ^ 2 + vector(i,2) ^ 2 + vector(i,3) ^ 2) ^ 0.5;
        vector(i,:)=vector(i,:)/vector_mag(i);    
    end
    vector2(:,:)=cat(2,vector(:,:),vector_mag(:));
    
end

