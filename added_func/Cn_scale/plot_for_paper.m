n_symm = 4;
%% create n_symm volume
im_size = 129;
for j = 10:20
gt_vol{j} = zeros(im_size,im_size,im_size);
x0 = round(im_size/2);
[X,Y] = meshgrid(linspace(-x0,x0,im_size),linspace(x0,-x0,im_size));
zmax = 15;
for z=0:zmax
plane=zeros(im_size,im_size);
m = double((X.^2+Y.^2).^0.5 < x0 * 0.8 * cos(z/zmax*pi/2)* cos(atan(Y./X)*n_symm+pi/10*j));
m = m .* (1+z);
for i=0:n_symm-1
plane = plane + imrotate(m,360/n_symm*i,'bilinear','crop');
end
gt_vol{j}(:,:,z+x0) = plane;
if z >0
    gt_vol{j}(:,:,-z+x0) = plane;
end
end
end
%%
gt_vol10 = gt_vol{10};
gt_vol11 = gt_vol{11};
gt_vol12 = gt_vol{12};
gt_vol13 = gt_vol{13};
gt_vol14 = gt_vol{14};
gt_vol15 = gt_vol{15};
gt_vol16 = gt_vol{16};
gt_vol17 = gt_vol{17};
gt_vol18 = gt_vol{18};
gt_vol19 = gt_vol{19};
gt_vol20 = gt_vol{20};
%%
[ rot ] = euler_to_rot( [-20,30,0] );
 

close all
figure;vol3d('Cdata',gt_vol9,'xdata', [-x0 x0], 'ydata', [-x0 x0], 'zdata', [-x0 x0]);view(3);
 hold on;
  quiver3([0;0;0],[0;0;0],[0;0;0],rot*[100;0;0],rot*[0;100;0],rot*[0;0;100])

%  quiver3([0;0;0],[0;0;0],[0;0;0],[100*cosd(30);0;-100*sind(30)],[0;100;0],[100*sind(30);0;100*cosd(30)])
 hold on;
 quiver3([0;0;0],[0;0;0],[0;0;0],[100;0;0],[0;100;0],[0;0;100])
view(70-20,35)
