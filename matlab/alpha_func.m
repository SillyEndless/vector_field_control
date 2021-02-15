function a = alpha_func(x,y,path,W)

    np = length(path(1,:));
    S = 0;
    for k = 1:1:np
        S = S + W(k)*f_func([x;y],path(:,k));
%         fprintf('%d\t%f\n',k,S)
    end
    a = 1-S;
end