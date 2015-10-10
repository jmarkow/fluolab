function [TRIALS]=fluolab_fb_proc(DATA,AUDIO,TTL,varargin)
%  fiberdata
%
%
%
%
%

TRIALS=[];

nparams=length(varargin);

if mod(nparams,2)>0
	error('Parameters must be specified as parameter/value pairs');
end

blanking=[0 0];
channel=1;
daf_level=.3;
trial_cut=2;
normalize='m';
newfs=400;
dff=1;
tau=.1;
detrend_win=.3;
classify_trials='t';
detrend_method='p';
tau_regress=.05;
nmads=4;
detrend_sliding=0;
neg=0;
padding=[];

for i=1:2:nparams
	switch lower(varargin{i})
		case 'blanking'
			blanking=varargin{i+1};
		case 'daf_level'
			daf_level=varargin{i+1};
		case 'normalize'
			normalize=varargin{i+1};
		case 'classify_trials'
			classify_trials=varargin{i+1};
		case 'newfs'
			newfs=varargin{i+1};
		case 'dff'
			dff=varargin{i+1};
		case 'tau'
			tau=varargin{i+1};
		case 'detrend_win'
			detrend_win=varargin{i+1};
		case 'channel'
			channel=varargin{i+1};
		case 'detrend_method'
			detrend_method=varargin{i+1};
		case 'tau_regress'
			tau_regress=varargin{i+1};
		case 'nmads'
			nmads=varargin{i+1};
		case 'trial_cut'
			trial_cut=varargin{i+1};
		case 'neg'
			neg=varargin{i+1};
		case 'padding'
			padding=varargin{i+1};
	end
end

if isempty(TTL) & isempty(AUDIO)
	classify_trials='';
end

proc_data=double(DATA.data(:,:,channel));

if neg
	proc_data=-proc_data;
end

[nsamples,ntrials]=size(proc_data);

trial_cut_idx=any(proc_data<trial_cut);
[~,bad_trial]=find(trial_cut_idx);


if nmads>0
	%[bad_trial2]=fluolab_hampel_filt(proc_data(blanking_idx,:),'nmads',nmads);
	% detect jumps in dt

	dt=max(abs(diff(proc_data(:,~trial_cut_idx))));
	mu=median(dt)
	v=mad(dt,2)

	tmp=dt>(mu+nmads*v);
	
	dt_idx=zeros(size(trial_cut_idx));
	dt_idx(~trial_cut_idx)=tmp;

	bad_trial=find(trial_cut_idx|dt_idx);

end

include_trials=setdiff(1:ntrials,unique(bad_trial));
ntrials=length(include_trials);

%pause();

% where are the feedback trials?

if isempty(classify_trials)
	C=zeros(size(DATA.data,2),1);
else
	[TRIALS,C]=fluolab_classify_trials(TTL,AUDIO,'include_trials',include_trials,...
		'method',classify_trials,'blanking',[0 0],'daf_level',daf_level,'padding',padding);
end

TRIALS.fluo_include.catch=find(C(include_trials)==2);
TRIALS.fluo_include.daf=find(C(include_trials)==1);
TRIALS.fluo_include.other=find(C(include_trials)==0);
TRIALS.fluo_include.catch_other=find(C(include_trials)==0|C(include_trials)==2);

TRIALS.fluo_include.all=[1:length(include_trials)];
TRIALS.all.fluo_include=include_trials;


