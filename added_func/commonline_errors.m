function [err] = commonline_errors(W,indeces,proj_num)
err = zeros(proj_num,1);
for i=1:size(indeces,1)
    
err(indeces(i,1)) = err(indeces(i,1)) + W(i);
err(indeces(i,2)) = err(indeces(i,2)) + W(i);
end
err = err / proj_num;
end

