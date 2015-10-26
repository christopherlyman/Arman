function Arman()
% Good reference file
% Two key things:
% (1) Convert Y to double (not single, or anything else) before multiplying
% (2) Map appropriately for SliceLocation
%

clear
clear global


% uiwait(msgbox('Please select directory you would like to convert.','SPA'));
main_dir = uigetdir('Z:\Hopkins-data\', 'Select directory you would like to convert...');

box = {'Enter Name of Output File:'};
box_title = 'SPA';
num_lines = 1;
InputFileName = inputdlg(box,box_title,num_lines);
if isempty(InputFileName{1}) == 1;
    OutFileName = ('NIfTI');
end
if isempty(InputFileName{1}) == 0;
    OutFileName = InputFileName{1};
end

% uiwait(msgbox('Please select Output directory.','SPA'));
out_dir = uigetdir(main_dir, 'Select directory for output...');


list = ls(main_dir);
struct1 = dir(main_dir);
listsize1 = size(list);
listlength1 = listsize1(1);
struct = struct1(3:listlength1,:);
listsize = size(struct);
listlength = listsize(1);
clear struct1;
clear listsize1;
clear listlength1;

% cell1 = struct2cell(struct);
% cell2 = cell1(1,:);
% cell = cell2.';
% clear cell1;
% clear cell2;

DICOM_names = {};

for kk = 1:1:listlength;
    isdir = struct(kk).isdir;
    if isdir == 0,
        name = [main_dir '\' struct(kk).name];
        filename = struct(kk).name;
        if isdicom(name) == 1,
            DICOM_size = size(DICOM_names,1)+1;
            DICOM_names{DICOM_size,1} = filename;
        end
    end
end
clear DICOM_size

DICOM_size = size(DICOM_names,1);

name = [main_dir '\' DICOM_names{1}];
info = dicominfo(name);
Y = dicomread(info);

if isfield(info,'NumberOfTimeSlices') == 1,
    frames = info.NumberOfTimeSlices;
else
    frames = 1;
end

xdim = info.Rows;
ydim = info.Columns;
zdim = info.NumberOfSlices;
z_width = info.SliceThickness;

x_width=info.PixelSpacing(1);
y_width=info.PixelSpacing(2);

image=zeros(xdim,ydim,zdim,frames);
image_final=zeros(xdim,ydim,zdim,frames);

baddata = 0;
Totalfiles = (frames*zdim);
check = Totalfiles/DICOM_size;
if check ~= 1,
    if check == 2,
        frames = 1;
    else
        uiwait(msgbox('The number of files does not match. Please check that you have all the DICOM data.','VWI'));
        baddata = 1;
    end
end

if baddata ~= 1,
    for frame=1:frames
        if frame == 1,
            for j=1:zdim   % Bottom to top
                slice_number = j+(frame-1)*zdim
                
                isdir = struct(j).isdir;
                if isdir == 0,
                    name = struct(j).name
                    fullname = [main_dir '\' struct(j).name];
                    if isdicom(fullname) == 1,
                        
                        info=dicominfo(sprintf('%s',fullname));
                        Y = dicomread(info);
                        
                        image(:,:,j,frame)=flipud(rot90((double(Y))*info.RescaleSlope+info.RescaleIntercept,1));
                        
                        SliceLocation(j,frame)=info.SliceLocation;
                        AcquisitionTime(j,frame)=str2num(info.AcquisitionTime);
                        FrameReferenceTime(j,frame)=info.FrameReferenceTime;
                    end
                end
                
            end
        else
            for j=zdim+1:1:listlength   % Bottom to top
                slice_number = j
                
                isdir = struct(j).isdir;
                if isdir == 0,
                    name = struct(j).name
                    fullname = [main_dir '\' struct(j).name];
                    if isdicom(fullname) == 1,
                        
                        info=dicominfo(sprintf('%s',fullname));
                        Y = dicomread(info);
                        
                        % We do rotation and flipping to correspond to Vinci DICOM images
                        image(:,:,j-zdim,frame)=flipud(rot90((double(Y))*info.RescaleSlope+info.RescaleIntercept,1));
                        
                        SliceLocation(j-zdim,frame)=info.SliceLocation;
                        AcquisitionTime(j-zdim,frame)=str2num(info.AcquisitionTime);
                        FrameReferenceTime(j-zdim,frame)=info.FrameReferenceTime;
                    end
                end
                
            end
        end
    end
    
    
    SliceThickness=double(info.SliceThickness);
    % SliceLocation;
    min_SliceLocation=min(min(SliceLocation));
    
    % Map to appropriate slice number
    for frame=1:frames;
        hold=round((SliceLocation(:,frame)-min_SliceLocation)/SliceThickness)+1;
        SliceLocationMap(:,frame)=double(zdim)-hold+1;
    end
    
    for j=1:zdim;
        [hold I]=sort(FrameReferenceTime(j,:));
        frame_map(j,:)=I;
    end
    
    
    % Sorts planes; also switches axis direction to correspond to Vinci
    % DICOM images
    for frame=1:frames;
        for j=1:zdim;   % Bottom to top
            image_final(:,:,SliceLocationMap(j,frame),frame)=image(:,:,j,frame_map(j,frame));
        end
    end
    
    
    for frame=1:frames;
        for j=1:zdim;   % Bottom to top
            image_final2(:,:,zdim-j+1,frame)=flipud(fliplr(image_final(:,:,j,frame)));
        end
    end
    
    % Saving all dynamic frames into a single NIFTI file
    %nifti_image=make_nii(image_final2,[x_width y_width z_width],16)
    %nifti_image.hdr.dime.datatype=16
    %nifti_image.hdr.dime.bitpix=16
    %save_nii(nifti_image,'nifti_image')
    
    
    % Puts origin at center of FOV by taking half the X,Y,Z dimensions.
    X = xdim/2;
    Y = ydim/2;
    Z = zdim/2;
    origin = [X Y Z];
    
    % Saving each dynamic frames into a separate NIFTI file
    for frame=1:frames;
        nifti_image=make_nii(image_final2(:,:,:,frame),[x_width y_width z_width],[origin],16);
        nifti_image.hdr.dime.datatype=16;
        nifti_image.hdr.dime.bitpix=16;
        if frames > 1,
            output_filename=sprintf('%s/%s%s%d%s',out_dir,OutFileName,'_fr-', frame, '.nii');
        else
            output_filename=sprintf('%s/%s%s',out_dir,OutFileName,'.nii');
        end
        save_nii(nifti_image,output_filename);
    end
    
    
    clear all
    clc
    
    disp('DONE!');
end

end