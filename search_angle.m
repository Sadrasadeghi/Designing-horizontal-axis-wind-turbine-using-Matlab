function [p]=search_angle(X,alpha)
    
    p=1;
    
    while (X(p)<alpha)
        p=p+1;
    end
    
    p=p-1;
end
        
