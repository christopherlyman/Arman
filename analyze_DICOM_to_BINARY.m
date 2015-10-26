%%% Good reference file
% Two key things:
% (1) Convert Y to double (not single, or anything else) before multiplying
% (2) Map appropriately for SliceLocation

clear
clear global

%shift=0;
%skip=0;

%main_dir='/Users/armanrahmim/linux/pet_programs/DICOM_read/GwennSmith/151802/';
%main_dir='/Users/armanrahmim/linux/pet_programs/DICOM_read/GwennSmith/152802/';
%main_dir='/Users/armanrahmim/linux/pet_programs/DICOM_read/GwennSmith/merged/';
main_dir='Z:\Hopkins-data\AD-DBS\!Sites\02_Hopkins\Raw-Data\FNMI-02-001\DICOM\PET_2\DICOM\20121228\PT_1509P\DCM_229E_AC\test';

cd(main_dir);

zero_string_list={'000','00','0'}';
ref_string='PT00';
info=dicominfo(sprintf('%s/%s%s1',main_dir,ref_string,zero_string_list{1}));
Y = dicomread(info);


xdim=info.Rows
ydim=info.Columns
zdim=info.NumberOfSlices
frames=info.NumberOfTimeSlices
z_width=info.SliceThickness
%x_width=double(info.ReconstructionDiameter)/double(info.Width)
x_width=info.PixelSpacing(1)
y_width=info.PixelSpacing(2)
%StudyTime=info.StudyTime %for entire study    
%InjectionTime=info.RadiopharmaceuticalInformationSequence.Item_1.RadiopharmaceuticalStartTime
%Dose=info.RadiopharmaceuticalInformationSequence.Item_1.RadionuclideTotalDose

image=zeros(xdim,ydim,zdim,frames);
image_final=zeros(xdim,ydim,zdim,frames);
%'pause'
%pause


for frame=1:frames
  for slice_number0=1:zdim   % Bottom to top 
    slice_number = slice_number0+(frame-1)*zdim;
    
    lookup_index=floor(log10(double(slice_number)))+1
    
    info=dicominfo(sprintf('%s/%s%s%d',main_dir,ref_string,zero_string_list{lookup_index},slice_number));
    Y = dicomread(info);

    % We do rotation and flipping to correspond to Vinci DICOM images
    image(:,:,slice_number0,frame)=flipud(rot90((double(Y))*info.RescaleSlope+info.RescaleIntercept,1));

    SliceLocation(slice_number0,frame)=info.SliceLocation;
    AcquisitionTime(slice_number0,frame)=str2num(info.AcquisitionTime);
    FrameReferenceTime(slice_number0,frame)=info.FrameReferenceTime;

  end
end

SliceThickness=double(info.SliceThickness);
SliceLocation
min_SliceLocation=min(min(SliceLocation));

%min_FrameReferenceTime=min(min(FrameReferenceTime));

%SliceNumbers=zeros(zdim,frames);

% Map to appropriate slice number
for frame=1:frames
  hold=round((SliceLocation(:,frame)-min_SliceLocation)/SliceThickness)+1;
  SliceLocationMap(:,frame)=double(zdim)-hold+1;
end

for slice_number0=1:zdim
  [hold I]=sort(FrameReferenceTime(slice_number0,:));
  frame_map(slice_number0,:)=I;
end


% Sorts planes; also switches axid direction to correspond to Vinci
% DICOM images
for frame=1:frames
  for slice_number0=1:zdim   % Bottom to top 
    image_final(:,:,SliceLocationMap(slice_number0,frame),frame)=image(:,:,slice_number0,frame_map(slice_number0,frame));
  end
end



% 
% figure(1)
% clf
% figure(2)
% clf

% for frame=1:frames
% 
%   figure(1)
%   subplot(1,double(frames),double(frame))
%   pcolor(image_final(:,:,ceil(zdim/2),frame)'); shading flat; colorbar
% 
%   figure(2)
%   subplot(1,double(frames),double(frame))
%   pcolor(squeeze(image_final(:,55,:,frame))'); shading flat; colorbar
%   
%   % Save as binary
%   fid = fopen(sprintf('image_final_frame%d.i',frame), 'wb');
%   fwrite(fid,image_final(:,:,:,frame),'float32');
% end

SliceLocation
SliceLocationMap
%AcquisitionTime
FrameReferenceTime
frame_map

for frame=1:frames
  for slice_number0=1:zdim   % Bottom to top 
    image_final2(:,:,zdim-slice_number0+1,frame)=flipud(fliplr(image_final(:,:,slice_number0,frame)));
  end
end

% Saving all dynamic frames into a single NIFTI file
%nifti_image=make_nii(image_final2,[x_width y_width z_width],16)
%nifti_image.hdr.dime.datatype=16
%nifti_image.hdr.dime.bitpix=16
%save_nii(nifti_image,'nifti_image')

% Saving each dynamic frames into a separate NIFTI file
for frame=1:frames
  nifti_image=make_nii(image_final2(:,:,:,frame),[x_width y_width z_width],16)
  nifti_image.hdr.dime.datatype=16
  nifti_image.hdr.dime.bitpix=16
  output_filename=sprintf('NIfTI_frame-%d.nii',frame);
  save_nii(nifti_image,output_filename)

end



