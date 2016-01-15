function [MAT POS_CROSS]=fluoflab_sliding_window(TTL,varargin)
% plot window centered on ttl pulse
%
%
%
%
MAT=[];

colors='jet';
cbar_dist=.1;
cbar_width = .025;
clim_order=1e3;
datenums=[];
time_dist=.005;
time_width=.075;
escape_dist=.005;
escape_width=.075;
escape_idx=[];
padding=[.2 .2];
thresh=.5;
channel=1;
binsize=.001;

nparams=length(varargin);

if mod(nparams,2)>0
	error('ephysPipeline:argChk','Parameters must be specified as parameter/value pairs!');
end

for i=1:2:nparams
	switch lower(varargin{i})
		case 'cbar_dist'
			cbar_dist=varargin{i+1};
		case 'colors'
			colors=varargin{i+1};
		case 'cbar_width'
			cbar_width=varargin{i+1};
		case 'datenums'
			datenums=varargin{i+1};
		case 'label'
			label=varargin{i+1};
		case 'escape_idx'
			escape_idx=varargin{i+1};
		case 'padding'
			padding=varargin{i+1};
		case 'channel'
			channel=varargin{i+1};
		case 'binsize'
			binsize=varargin{i+1};
	end
end

[nsamples,ntrials,nchannels]=size(TTL.data);

% get onsets

idx=(1:nsamples-1);

% POS_CROSS

POS_CROSS=zeros(1,ntrials);

padding_smps=round(padding.*TTL.fs);

new_fs=1/binsize;
new_samples=round((nsamples/TTL.fs)*new_fs);
new_t=[0:new_samples-1]/new_fs;

MAT=zeros(new_samples,ntrials);

for i=1:ntrials

	cur_ttl=TTL.data(:,i,channel);

	hit=TTL.t(find(cur_ttl(idx)<thresh&cur_ttl(idx+1)>=thresh));
	hit(hit<padding(1)|hit>TTL.t(end)-padding(2))=[];

	% map the hit to the nearest timepoint in the singing data

	if isempty(hit)
		continue;
	end

	hit=hit(1)
	[~,loc]=min(abs(new_t-hit));

	POS_CROSS(i)=loc;
	MAT(loc,i)=1;

end
