function ImStack = load_spectra_stack(ImSize)

[FileNames,path]=uigetfile(['*','.IMA'],'Select Spectra Files','multiselect','on');

num_files=length(FileNames);
all_slices = cell(num_files,1);
for i = 1:num_files
    
    obj = dicom_open(fullfile(path,FileNames{i}));
    tmp_spectra = dicom_get_spectrum_siemens(obj);
    tmp_spectra = reshape(tmp_spectra,[], ImSize(1), ImSize(2));
    if i == 1
        ImStack = zeros([size(tmp_spectra),length(FileNames)]);
    end
    ImStack(:,:,:,i) = tmp_spectra;
end
