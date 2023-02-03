function [nucleate,growing,largos]=nucleate(growing,pnuc,num_tubules,dt,largos);
    nucleate = [];
    largos=largos;
    growing = growing;
    randoms=rand(num_tubules,1);
    counter = 0;
    for i=1:num_tubules
        if growing(i,1)==0
            if randoms(i,1) <= pnuc
                growing(i,1)=1;
                largos(i,1)=0;
                counter = counter + 1;
                nucleate(counter) = i;
            end
        end
    end;
end
