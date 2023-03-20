function [disfid,gasfid,traj,gastraj,H1_im,Gas_Image,H1_Mask,Vent_Mask,RBC_Mask,RBC2Bar,TR,ImSize,scanDateStr,write_path] = debug_wiggles(file)

%Quick function to make my life easier when debugging wiggles

if nargin<1
    [myfile,mypath]=uigetfile();
    file = fullfile(mypath,myfile);
    my_write_path = mypath;
else
    idcs = strfind(file,filesep);%determine location of file separators
    my_write_path = file(1:idcs(end)-1);%remove file
end

load(file);


try
    disfid = Dis_Fid;
    gasfid = Gas_Fid;
    traj = Dis_Traj;
    gastraj = Gas_Traj;
    H1_im = H1_Image_Dis;
    Gas_Image = LoRes_Gas_Image;
    H1_Mask = Proton_Mask;
    Vent_Mask = VentBinMask;
    RBC_Mask = RBC_Mask;
    RBC2Bar = -RBC2Bar;
    TR = TR;
    ImSize = ImSize;
    scanDateStr = scanDateStr;
    write_path = my_write_path;
catch
    disfid = Xe_Raw_Dis;
    gasfid = Xe_Raw_Gas;
    traj = Xe_Traj;
    gastraj = Xe_Traj;
    H1_im = H1_Image;
    Gas_Image = LoRes_Gas_Image;
    H1_Mask = Proton_Mask;
    Vent_Mask = VentBinMask;
    RBC_Mask = RBC_Mask;
    RBC2Bar = RBC2Bar;
    TR = TR;
    ImSize = ImSize;
    scanDateStr = scanDateStr;
    write_path = my_write_path;
end



