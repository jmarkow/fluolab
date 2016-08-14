function DATA=fluolab_datascrub(DATA,varargin)
%
%
%
%



channel=1;
trial_cut=3;
trial_cut_factor=.75;
nmads=0;

nparams=length(varargin);
if mod(nparams,2)>0
	error('ephysPipeline:argChk','Parameters must be specified as parameter/value pairs!');
end

for i=1:2:nparams
	switch lower(varargin{i})
		case 'channel'
			channel=varargin{i+1};
		case 'trial_cut'
			trial_cut=varargin{i+1};
		case 'trial_cut_factor'
			trial_cut_factor=varargin{i+1};
		case 'nmads'
			nmads=varargin{i+1};
  end
end

proc_data=double(DATA.data(:,:,channel));

[nsamples,ntrials,nchannels]=size(proc_data);

if strcmp(trial_cut,'auto')
	trial_cut=trial_cut_factor*median(max(proc_data));
	fprintf('Trial cut: %g\n',trial_cut);
end

trial_cut_idx=any(proc_data(:,:,1)<trial_cut);
[~,bad_trial]=find(trial_cut_idx);

if nmads>0

	dt=max(abs(diff(proc_data(:,~trial_cut_idx,1))));
	mu=median(dt);
	v=mad(dt,2);

	tmp=dt>(mu+nmads*v);

	dt_idx=zeros(size(trial_cut_idx));
	dt_idx(~trial_cut_idx)=tmp;

	bad_trial=find(trial_cut_idx|dt_idx);

end

include_trials=setdiff(1:ntrials,unique(bad_trial));
ntrials=length(include_trials);

% preserve included trials, preserve trial indices

DATA.data=DATA.data(:,include_trials,:);
DATA.trial_idx=include_trials;
