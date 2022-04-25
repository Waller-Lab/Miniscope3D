function stack = read_tiff_stack(path,varargin)
%stack = read_tiff_stack(path,downsample_ratio,list_of_images)
%list_of_images is a vector of indices to read out. if empty, read all.



info = imfinfo(path);
num_images_in = numel(info);

if nargin>1
    ds = varargin{1};
    if nargin > 2
        k_list = varargin{2};
    else
        k_list = 1:num_images_in;
    end
else
    ds = 1;
end
num_planes = length(k_list);
if info(1).SamplesPerPixel == 1
    
    stack = zeros(info(1).Height/ds, info(1).Width/ds,num_planes);
else
    stack = zeros(info(1).Height/ds, info(1).Width/ds,info(1).SamplesPerPixel,num_planes);
end
n = 0;
for k = k_list
    n = n+1;
    if info(1).SamplesPerPixel == 1
        stack(:,:,n) = imresize(imread(path, k, 'Info', info),1/ds,'box');
    else
        stack(:,:,:,n) = imresize(imread(path, k, 'Info', info),1/ds,'box');
    end
    
end