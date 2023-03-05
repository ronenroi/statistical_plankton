function stackedProjections = func_loadPlanktonData(dataFolder, projectionType)

jpgImageFiles  = dir([dataFolder  '/*.jpeg']);   
pngImageFiles  = dir([dataFolder  '/*.png']);
if isempty(jpgImageFiles)
    imageFiles = pngImageFiles;
else 
    imageFiles = jpgImageFiles;
end
    
nproj = length(imageFiles);    % Number of files found

max_imsize = 0;
min_val = 999;
for ii=1:nproj
    currentfilename = imageFiles(ii).name;
    currentImage = imread([dataFolder '/' currentfilename]);

    if strcmp(projectionType, 'BrightField')
        image = -log(double(currentImage)/255.0);
    else
        image = double(rgb2gray(currentImage));
    end
    
    % crop out the Plankton in the image
    image = image.*(image>graythresh(image));
    [y,x] = ind2sub(size(image), find(image));
    coord = [x, y];
    mc = min(coord)-0.5;
    Mc = max(coord)+0.5;
    bbox = [mc Mc-mc]; 
    images{ii} = imcrop(image, bbox);
    max_imsize = max(max_imsize, max(size(images{ii})));
end

max_imsize = max_imsize + ~mod(max_imsize,2);
stackedProjections = zeros([max_imsize,max_imsize, nproj]);
for ii=1:nproj
   img = images{ii};
   imsize = size(img);
   padded = padarray(img, ceil((max_imsize-imsize)/2), 'pre');
   padded = padarray(padded, floor((max_imsize-imsize)/2), 'post');
   stackedProjections(:,:,ii) = padded;
end

end

