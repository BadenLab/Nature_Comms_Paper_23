function Colour_noise_temp = load_noise_from_hdf5 (hdf5_file, varargin)
%This function loads the noise frequency from the hdf5 file. The filename 
%and location of the noise is a required input argument. 
%By default, the whole sequence is load into the workspace, but
%based on user input only a subset of Noise can be loaded as well. Two
%options are available: Give number of boxes as input, or set askforboxnr =
%true in the input
%@MSeifert 2020
%
%
%Input: (Noise_file,askforboxnr,(true or
%false),nr_boxes,nr_frames,start_loc,[frame,box])

defaultAnswer = false;
defaultStartloc = [1 1];
defaultNrFrames = Inf;
defaultNrBoxes = Inf;


p = inputParser;
validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
validMatrixShape = @(x) isnumeric(x) && (size(x,1) ==1) && (size(x,2) == 2);
addRequired(p,'noise_path',@isstring);
addOptional(p,'askforboxnr',defaultAnswer,@islogical);
addOptional(p,'nr_boxes',defaultNrBoxes,validScalarPosNum);
addOptional(p,'nr_frames',defaultNrFrames,validScalarPosNum);
addOptional(p,'start_loc', defaultStartloc,validMatrixShape)



parse(p,hdf5_file,varargin{:});

%out = p.Results

%Check if hdf5 file exists at location
if exist(hdf5_file,'file') ~= 2
    error("Hdf5 file specified does not exist in the directory")
end


%Ask for the number of boxes if whished for
if p.Results.askforboxnr == 1
    
    prompt = {'Enter number of boxes in x:','Enter number of boxes in y:'};
    dlgtitle = 'Input';
    dims = [1 35];
    definput = {'40','40'};
    nr_boxes_temp = inputdlg(prompt,dlgtitle,dims,definput);
    nr_boxes = str2double(nr_boxes_temp{1});
    nr_boxes(2) = str2double(nr_boxes_temp{2});
    nr_boxes = nr_boxes(1)*nr_boxes(2); %This is the total number of boxes 
    %in both dimensions which we can use to load the hdf5 file
else
    nr_boxes = p.Results.nr_boxes;
end

%If no nr of frames are given, load all frames
if p.Results.nr_frames == Inf
try
file_info = h5read(hdf5_file, '/Info');
nr_frames = file_info(3);
catch
    disp('No info variable found in hdf5 file, add info variable, or define nr_frames')
end
else
    nr_frames = p.Results.nr_frames;
end

%If nr of boxes isnt given, load all frames
if p.Results.nr_boxes == Inf
try
file_info = h5read(hdf5_file, '/Info');
nr_boxes = file_info(1)*file_info(2);
catch
    disp('No info variable found in hdf5 file, add info variable, define nr boxes or set askforboxnr = true')
end
end

%Load the data

hdf5_start = p.Results.start_loc;
hdf5_count = ceil([nr_frames nr_boxes]);

Colour_noise_temp = zeros([uint16(hdf5_count),4]);
Colour_noise_temp(:,:,1) = h5read(hdf5_file, '/Red_Noise', hdf5_start, hdf5_count);
Colour_noise_temp(:,:,2) = h5read(hdf5_file, '/Green_Noise', hdf5_start, hdf5_count);
Colour_noise_temp(:,:,3) = h5read(hdf5_file, '/Blue_Noise', hdf5_start, hdf5_count);
Colour_noise_temp(:,:,4) = h5read(hdf5_file, '/UV_Noise', hdf5_start, hdf5_count);








    
    


end