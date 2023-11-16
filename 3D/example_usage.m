% Must have eirods installed "https://eidors3d.sourceforge.net/"

% Ensure the config file is up do date

% Create an inverse model
shape_str = 'pig_23kg';
imdl = mk_3d_imdl(shape_str);

% Load ventilation and perfusion data
ventilation_file = "test_ventilation_file.eit";
perfusion_file = "test_perfusion_file.eit";
pth = '~/data/2021/Araos-Cornell/Pig-3D-EIT/Aston Martin/';
ventilation_file = [pth,'AE-B1.eit'];
perfusion_file   = [pth,'AE-B1H.eit'];

% Ventilation image reconstruction
[dd, ~] = eidors_readdata(ventilation_file,'LQ4');
dd = real(dd); 
%fs = aux.frame_rate;
window = 10;
ventilation_imgr = calc_ventilation(imdl,dd,window);

% perfusion image reconstruction
[dd, aux] = eidors_readdata(perfusion_file,'LQ4');
dd = real(dd); 
fs = aux.frame_rate;
perfusion_imgr = calc_perfusion(imdl,dd,fs);


