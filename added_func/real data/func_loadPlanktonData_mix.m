function stackedProjections = func_loadPlanktonData_mix(dataFolder, projectionType)

jpgImageFiles  = dir([dataFolder  '/main/*.jpeg']);   
pngImageFiles  = dir([dataFolder  '/main/*.png']);
if isempty(jpgImageFiles)
    imageFiles = pngImageFiles;
else 
    imageFiles = jpgImageFiles;
end
    
nproj = length(imageFiles);    % Number of files found

max_imsize = 0;
for ii=1:nproj
    currentfilename = imageFiles(ii).name;
    currentImage = imread([dataFolder '/main/' currentfilename]);

    if strcmp(projectionType, 'BrightField')
        image = -log(double(currentImage)/255.0);
    else
        image = double(rgb2gray(currentImage));
    end
    
    % crop out the Plankton in the image
    image = image.*(image>graythresh(image)/1.3);
    [y,x] = ind2sub(size(image), find(image));
    coord = [x, y];
    mc = min(coord)-0.5;
    Mc = max(coord)+0.5;
    bbox = [mc Mc-mc]; 
    images{ii} = imcrop(image, bbox);
    max_imsize = max(max_imsize, max(size(images{ii})));
end


jpgImageFiles  = dir([dataFolder  '/mix/*.jpeg']);   
pngImageFiles  = dir([dataFolder  '/mix/*.png']);
if isempty(jpgImageFiles)
    imageFiles = pngImageFiles;
else 
    imageFiles = jpgImageFiles;
end
    
nproj1 = length(imageFiles);    % Number of files found

for ii=1:nproj1
    currentfilename = imageFiles(ii).name;
    currentImage = imread([dataFolder '/mix/' currentfilename]);

    if strcmp(projectionType, 'BrightField')
        image = -log(double(currentImage)/255.0);
    else
        image = double(rgb2gray(currentImage));
    end
    
    % crop out the Plankton in the image
    image = image.*(image>graythresh(image)/1.3);
    [y,x] = ind2sub(size(image), find(image));
    coord = [x, y];
    mc = min(coord)-0.5;
    Mc = max(coord)+0.5;
    bbox = [mc Mc-mc]; 
    images1{ii} = imcrop(image, bbox);
    max_imsize = max(max_imsize, max(size(images1{ii})));
end




max_imsize = max_imsize + ~mod(max_imsize,2);
stackedProjections = zeros([max_imsize,max_imsize, nproj+nproj1]);
for ii=1:nproj
   img = images{ii};
   imsize = size(img);
   padded = padarray(img, ceil((max_imsize-imsize)/2), 'pre');
   padded = padarray(padded, floor((max_imsize-imsize)/2), 'post');
   stackedProjections(:,:,ii) = padded;
end

for ii=1:nproj1
   img = images1{ii};
   imsize = size(img);
   padded = padarray(img, ceil((max_imsize-imsize)/2), 'pre');
   padded = padarray(padded, floor((max_imsize-imsize)/2), 'post');
   stackedProjections(:,:,ii+nproj) = padded;
end

end

