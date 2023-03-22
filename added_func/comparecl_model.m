function [ thetadiff ] = comparecl_model( clstack, ...
    ref_clmatrix, n_theta )
% compare the estimated clstack to the true ref_clmatrix and compute the
% correctness.
%
% Input:
%   clstack: estimated commonline stack.
%   ref_clmatrix: true commonline stack.
%   n_theta: angular resolution.
%   max_angle: maximum tolerant error
%
% Output:
%   prob: the probability of good commonlines
%   clean_clmatrix: the commonline stack with good commonlines
%   n_cl: the number of good commonlines for each image.
%
% Lanhui Wang, July 3, 2013 (Based on cryo_clmatrix_v3.m by Yoel)


n_proj = size(clstack,1);
res=(360/n_theta);

thetadiff=zeros(size(clstack));

for k1=1:n_proj
    for k2=k1+1:n_proj
        tcl1=ref_clmatrix(k1,k2);
        tcl2=ref_clmatrix(k2,k1);
            l1=clstack(k1,k2)*res;
     l2=clstack(k2,k1)*res;
        d1s=((l1-1)-(tcl1-1));
        d2s=((l2-1)-(tcl2-1));
        
        
        % Estimation error in angles
        if d1s <= 90 && d1s>-90
            thetadiff(k1,k2)=d1s;
        elseif d1s >= 90 && d1s < 270
            thetadiff(k1,k2)=d1s-180;
        elseif d1s >= 270 && d1s <= 360
            thetadiff(k1,k2)=d1s-360;
            
        elseif d1s <= -90 && d1s>-270
            thetadiff(k1,k2)=d1s+180;
            
        elseif d1s <= -270 && d1s >= -360
            thetadiff(k1,k2)=d1s+360;
        else
            d1s
        end
        % Estimation error in angles
        if d2s <= 90 && d2s>-90
            thetadiff(k2,k1)=d2s;
        elseif d2s >= 90 && d2s < 270
            thetadiff(k2,k1)=d2s-180;
        elseif d2s >= 270 && d2s <= 360
            thetadiff(k2,k1)=d2s-360;
            
        elseif d2s <= -90 && d2s>-270
            thetadiff(k2,k1)=d2s+180;
            
        elseif d2s <= -270 && d2s >= -360
            thetadiff(k2,k1)=d2s+360;
        else
            d2s
        end

    end
end
end


