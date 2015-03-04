function [TRIALS]=fluolab_fb_proc(AUDIO,TTL,varargin)
%
%
%
%
%
%

range=[]; % check for hits in a specific time range?
catch_pattern=[]; % check for this pattern in the TTL signal

nparams=length(varargin);

if mod(nparams,2)>0
	error('Parameters must be specified as parameter/value pairs');
end

blanking=[.05 .05];
channel=1;
daf_level=20;
trial_cut=2;
ttl_level=1;
normalize=1;

for i=1:2:nparams
	switch lower(varargin{i})
		case 'blanking'
			blanking=varargin{i+1};
		case 'daf_level'
			daf_level=varargin{i+1};
		case 'normalize'
			normalize=varargin{i+1};
	end
end

[nsamples,ntrials]=size(AUDIO.data);

daf_mat=TTL.data.^2>ttl_level;
[~,daf_trials]=find(daf_mat);
daf_trials=unique(daf_trials);
first_daf=min(min(daf_trials),ntrials+1);

% catch trials
%
TRIALS.catch=setdiff(first_daf:ntrials,daf_trials);
TRIALS.pre=1:max((first_daf-1),1);
TRIALS.daf=daf_trials';
