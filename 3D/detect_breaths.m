function [breath] = detect_breaths(dd,w)
	% Really really basic filtering and peak detection
	if size(dd,1) == 1
		breath_sig = dd;
		breath_sig = movmean(breath_sig,w);
	else
		breath_sig = sum(dd);
	end
	breath_sig = rescale(breath_sig,0,1);
	prominence = 0.25; 
	pks = islocalmax(breath_sig,'MinProminence',prominence,'FlatSelection','first');
	p = prominence;
	while isempty(pks)
		p = p-0.05;
		if p < 0
			breath = [];
			return
		end
		pks = islocalmax(breath_sig,'MinProminence',p,'FlatSelection','first');
	end

	pks = find(pks);

	% Find the troughs using the signal min (end expiration)
	trghs_sig = -breath_sig;
	trghs = islocalmax(trghs_sig,'MinProminence',prominence,'FlatSelection','first');
	p = prominence;
	while isempty(pks)
		p = p-0.05;
		if p < 0
			breath = [];
			return
		end
		trghs = islocalmax(trghs_sig,'MinProminence',p,'FlatSelection','first');
	end
	trghs = find(trghs);
	if isempty(trghs)
		breath = [];
		return
	end
	trghs(1) = [];
	pks = pks(pks<max(trghs) & pks>min(trghs));
	if isempty(pks)
		breath = [];
		return
	end
	breath = zeros(1,length(pks));
	for i=1:numel(pks)
		breath(i).pk = pks(i);
		trghs_diff = trghs-pks(i);
		low_diffs = trghs_diff;
		high_diffs = trghs_diff;
		low_diffs(low_diffs>0) = inf;
		high_diffs(high_diffs<0) = inf;
		[~,T1] = min(abs(low_diffs));
		[~,T2] = min(high_diffs);
		breath(i).trgh(:) = [trghs(T1), trghs(T2)];
		length_inhale(i) = breath(i).pk - breath(i).trgh(1);
		length_exhale(i) = breath(i).trgh(2) - breath(i).pk; 
	end
	% Remove high outliers in breath lengths - assume these are not regular breaths
	while max(length_exhale)>1.5*mean(length_exhale) 
		[~,reject] = max(length_exhale);
		length_exhale(reject) = [];
		length_inhale(reject) = [];
		breath(reject) = [];
		pks(reject) = [];
	end