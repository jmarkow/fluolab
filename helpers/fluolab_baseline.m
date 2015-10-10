function BASELINE=fluoflab_sliding_window(SOURCE_FILES,varargin)
% plot window centered on ttl pulse
%
%
%
%

thinning=.3;
led_cut=3;
channel=1;
per=5;

nparams=length(varargin);

if mod(nparams,2)>0
	error('ephysPipeline:argChk','Parameters must be specified as parameter/value pairs!');
end

for i=1:2:nparams
	switch lower(varargin{i})
		case 'thinning'
			thinning=varargin{i+1};
		case 'led_cut'
			led_cut=varargin{i+1};
		case 'channel'
			channel=varargin{i+1};
		case 'per'
			per=varargin{i+1};
	end
end

BASELINE=nan(1,length(SOURCE_FILES));

[uniq_names,~,uniq_idx]=unique(SOURCE_FILES);

for i=1:length(uniq_names)

	disp([uniq_names{i}]);

	% get the source data

	load(uniq_names{i},'adc');

	% get points where the LED is on

	tmp=adc.data>led_cut;
	nsamples=length(adc.data);

	% convert to thinned idxs

	if all(tmp)
		thinned_idxs=[ round(thinning*adc.fs) nsamples-round(thinning*adc.fs)];
	else
		thinned_idxs=markolab_collate_idxs(tmp,round(-thinning*adc.fs));	
	end

	on_points=zeros(size(tmp));

	for j=1:size(thinned_idxs,1)
		on_points(thinned_idxs(j,1):thinned_idxs(j,2))=1;
	end

	% case to log

	on_points=on_points==1;
	length(find(on_points))

	BASELINE(uniq_idx==i)=prctile(adc.data(on_points),11);

end

