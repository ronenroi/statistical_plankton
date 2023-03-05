function compareproj(gt_projs,est_projs,nrows,ncols,tit,MSE,showlabels)
%
% Display the first nrows*ncols projections from the stack projs.
%
% Input:
%   nrows       Number of rows to display.
%   ncols       Number of columns to display.
%   showlabels  Nonzero to show image numbering labels (default 1).
%
% Example:
%   [projs, noisy_projs, ~, ref_q] = cryo_gen_projections(65,10,1/8,0,1);
%   viewstack(projs,5,5)
%
% Yoel Shkolnisky, July 2013.
n1 = size(gt_projs,1) + size(est_projs,1);
n2 = size(gt_projs,2);
n3 = size(gt_projs,3);
projs = zeros([n1, n2, n3]);
projs(1:size(gt_projs,1),:,:) = gt_projs;
projs(size(gt_projs,1)+1:end,:,:) = est_projs;
if ~exist('showlabels','var')
    showlabels=1;
end

particles=draw(projs,1,nrows,ncols,tit,MSE,showlabels);
sz=size(projs);

nticks=ceil(sz(3)/(nrows*ncols))-1;
if nticks>1
    %disp([sz(3) nrows ncols nticks]);
    
    ax=get(particles,'Parent');
    pos=get(ax,'Position');
    units=get(ax,'Units');
    sld = uicontrol('Parent',gcf,'Style', 'slider',...
        'Min',0,'Max',nticks,'Value',0,...
        'SliderStep',[1/nticks 1/nticks],'Units',units,...
        'Position', [pos(1) 0.02 pos(3)*0.75 0.05],...
        'UserData',struct('projs',projs,'nrows',nrows,'ncols',ncols,'tit',tit,'MSE',MSE,'showlabels',showlabels),...
        'Callback', @sldCallback);
    align([ax,sld],'Center','None');
end


end


function particles=draw(projs,startidx,nrows,ncols,tit,MSE,showlabels)
sz=size(projs);
borderwidth=2;
img=zeros(nrows*sz(1)+borderwidth*(nrows+1),ncols*sz(2)+borderwidth*(ncols+1),3); % RGB matrix

% Draw red borders
loc=1;
for k=1:nrows+1
    img(loc:loc+borderwidth-1,:,1)=1;
    loc=loc+sz(1)+borderwidth;
end

loc=1;
for k=1:ncols+1
    img(:,loc:loc+borderwidth-1,1)=1;
    loc=loc+sz(2)+borderwidth;
end

% Put in the images
idx=startidx;
colloc=borderwidth+1;
for k=1:ncols
    rowloc=borderwidth+1;
    for j=1:nrows
        if idx>sz(3) % No more images left
            continue;
        end
        p=projs(:,:,idx);
        
        % Scale image to be between 0 and 1
        pmin=min(p(:));
        pmax=max(p(:));
        p=(p-pmin)./(pmax-pmin);
        
        % Put image in the right place in the grid.
        img(rowloc:rowloc+sz(1)-1,colloc:colloc+sz(2)-1,1)=p;
        img(rowloc:rowloc+sz(1)-1,colloc:colloc+sz(2)-1,2)=p;
        img(rowloc:rowloc+sz(1)-1,colloc:colloc+sz(2)-1,3)=p;
        rowloc=rowloc+sz(1)+borderwidth;
        idx=idx+1;
    end
    colloc=colloc+sz(2)+borderwidth;
end

particles=imagesc(im2double(img));
axis off
axis image
title(tit);

if showlabels
    
    % Add text labels
    idx=startidx;
    colloc=borderwidth+1;
    for k=1:ncols
        rowloc=borderwidth+1;
        for j=1:nrows
            if idx>sz(3) % No more images left
                continue;
            end
            ht=text(colloc+5*nrows,rowloc+5*ncols,[num2str(idx), ' MSE = ',num2str(round(MSE(idx),3))]);
            set(ht,'FontUnits','normalized');
            set(ht,'FontWeight','bold');
            %set(ht,'FontSize',40);
            set(ht,'Color','b');
            set(ht,'BackgroundColor','y');
            rowloc=rowloc+sz(1)+borderwidth;
            idx=idx+1;
        end
        colloc=colloc+sz(2)+borderwidth;
    end
end
end

function sldCallback(source,event)
figure(get(source,'parent'));
val = source.Value;
% For R2014a and earlier:
% val = get(source,'Value');
ud=source.UserData;
startidx=floor(val*ud.nrows*ud.ncols+1);
draw(ud.projs,startidx,ud.nrows,ud.ncols,ud.tit,ud.MSE,ud.showlabels);

end
