function imgr_perf = calc_perfusion(imdl,dd,fs)
%calc_perfusion - Calculate the perfusion given EIT data and a bolus injection time
% 
% Requires: 
% imgr -> inverse model
% dd -> data
% fs -> sampling rate of the data
	% Get values from the config file
	fname = ("config.json");
	fid = fopen(fname); 
	raw = fread(fid,inf); 
	str = char(raw'); 
	fclose(fid); 
	data = jsondecode(str);
	bolus_inj = floor(data.bolus_injection*fs);
	apnoea_start = floor(data.apnoea_start*fs);
	card_duration = floor(data.card_time*fs);
	resp_duration = floor(data.resp_time*fs);
	p = [];
	p.heart_x = data.heart_pixel_x;
	p.heart_y = data.heart_pixel_y;
	p.lung_x = data.lung_pixel_x;
	p.lung_y = data.lung_pixel_y;
	if data.apnoea_start < 1
		apnoea_start = apnoea_start + fs; % Can't be the very start of file
	end
	% Reconstruct the images with a reference set as the time before bolus injection
	dd = detrend(dd')'; 
	dd = freq_filt(dd,@(f) (f<0.8),fs,2);
	% note I am dropping the first second from the ref to try and remove any filtering nonsnse 
	imgr = inv_solve(imdl,mean(dd(:,apnoea_start:bolus_inj),2),dd); % Use the mean of the pre-injection time as reference

	% Generate an image of a 2d plane to identify a pixel corresponding to the heart and lung bolus passing
	test_plane = calc_slices(imgr,[inf inf 1]);
	heart_pix = squeeze(test_plane(p.heart_y,p.heart_x,:));
	lung_pix = squeeze(test_plane(p.lung_y,p.lung_x,:));

	% Crop the data to get some regions of interest
	cropped_data = lung_pix(bolus_inj:resp_duration)'; % pulmonary data first
	seg_offset = min(cropped_data)*-1; % Get the distance from 0
	pulm_sig = cropped_data+seg_offset; % Set the data to an offset of 0
	x1 = 1:numel(cropped_data)+numel(lung_pix(resp_duration:end))-1; % Add a tail to keep things smooth

	% Cardiac pixel signal
%%% DISCUSSION WITH SYMON - 20231116
%%% card_duration must be larger than bolus_inj
	cropped_data = heart_pix(bolus_inj:card_duration)'; % pulmonary data first
	seg_offset = min(cropped_data)*-1; % Get the distance from 0
	card_sig = cropped_data+seg_offset; % Set the data to an offset of 0
	%x2 = 1:numel(cropped_data)+numel(heart_pix(c_t:end))-1; % Add a tail to keep things smooth

	% Now fit curves to the manually selected points to get the desired gama PDF functions!
	p_lung = nlinfit(1:numel(pulm_sig),pulm_sig,@h,[10 10 10]);
	lung_fit = h(p_lung,x1);
	p_heart = nlinfit(1:numel(card_sig),card_sig,@h,[10 10 10]);
	%heart_fit = h(p_heart,x2);
	[~,reconst_loc] = max(lung_fit);

	global g_params % This looks silly, but required (see function set_gamma_params())
	g_params.k1 = p_lung(1);
	g_params.th1 = p_lung(2);
	g_params.k2 = p_heart(1);
	g_params.th2 = p_heart(2);
	
	% now go voxel by voxel, fit the curve with predefined characteristics 
	% TODO test that this is working after coping from tested source
	imgr_mod = imgr; % duplicate the image structure
	for i=1:size(imgr.elem_data,1) % This is a real big image...
		% 1. Fit a gamma curve to the selected time window of the pixel
		% Select the window of data from the injection onwards
		pix_dat = imgr.elem_data(i,bolus_inj:end);
		x = 1:numel(pix_dat);
		p = nlinfit(1:numel(pix_dat),pix_dat,@g,[1 1]);
		%ex_fit = g(p,x );
		imgr_mod.elem_data(i,bolus_inj:end) = pix_dat - g([0 1],x)*p(1)*10000; 
	end

	imgr_mod.get_img_data.frame_select = reconst_loc + bolus_inj ;
	imgr_perf = imgr_mod;
	imgr_perf.elem_data = imgr.elem_data(:,reconst_loc + bolus_inj);
end

function y = h(abc,x)
	a = abc(1); b = abc(2); c = abc(3);
	y = c * x.^(a-1) .* exp(-x/b) / (b^a * gamma(a));
end

function yz = g(ab, x)
	[k1,th1,k2,th2] = set_gamma_params();
    a = ab(1); b = ab(2); 
	y = x.^(k1-1) .* exp(-x/th1) / (th1^k1 * gamma(k1)); % Heart component
	z = x.^(k2-1) .* exp(-x/th2) / (th2^k2 * gamma(k2)); % Lung component 
	yz = a*y + b*z;
end

function [k1,th1,k2,th2] = set_gamma_params()
    % This global variable thing is a bit weird but not many other
    % ways to do it with the anonamous functions???
	% If any suggestions let me know
    global g_params
    k1 = g_params.k1;
    th1 = g_params.th1;
    k2 = g_params.k2;
    th2 = g_params.th2;
end
