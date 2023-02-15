function imgr = calc_ventilation(imdl,dd,w)
%calc_ventilation - Description
%
% Input
% imdl -> the inverse model for image reconstruciton
% fs -> sampling rate of the EIT data
% w -> moving mean filter size for breath detection
% 
% Returns:
% imgr -> reconstructed EIT image of the EA breaths

	[breath] = detect_breaths(dd,w);
	% Reconstruct the averaged breath
	imgr = inv_solve(imdl,mean(dd(:,unique([breath(:).trgh])),2),mean(dd(:,[breath(:).pk]),2));

end