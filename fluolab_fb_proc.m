function [RAW,REGRESS,TRIALS]=fluolab_fb_proc(DATA,AUDIO,varargin)
%  fiberdata
%
%
%
%
%
nparams=length(varargin);

if mod(nparams,2)>0
	error('Parameters must be specified as parameter/value pairs');
end

blanking=[.2 .2];
channel=1;
daf_level=20;
trial_cut=2;

for i=1:2:nparams
	switch lower(varargin{i})
		case 'blanking'
			blanking=varargin{i+1};
		case 'daf_level'
			daf_level=varargin{i+1};
	end
end

proc_data=double(DATA.data(:,:,channel));
[nsamples,ntrials]=size(proc_data);

blanking_idx=[round(blanking(1)*DATA.fs):nsamples-round(blanking(2)*DATA.fs)];

% include these trials

[~,bad_trial]=find(proc_data(blanking_idx,:)<trial_cut);
include_trials=setdiff(1:ntrials,unique(bad_trial));
ntrials=length(include_trials);

% where are the feedback trials?

daf_mat=(AUDIO.data(:,include_trials).^2>daf_level);
[~,daf_trials]=find(daf_mat);
daf_trials=unique(daf_trials);
first_daf=min(min(daf_trials),ntrials);

% catch trials

TRIALS.catch=setdiff(first_daf:ntrials,daf_trials);
TRIALS.pre=1:max((first_daf-1),1);
TRIALS.daf=daf_trials';
TRIALS.all=1:ntrials;

trial_types=fieldnames(TRIALS);
ntypes=length(trial_types);

[new_data,time]=fluolab_condition(proc_data(:,include_trials),DATA.fs,'tau',.05,'detrend_win',.15,'newfs',100,'blanking',blanking);
[nsamples,ntrials]=size(new_data);

for i=1:ntypes
	if isempty(TRIALS.(trial_types{i})), continue; end
	RAW.ci.(trial_types{i})=bootci(1e3,{@mean,new_data(:,TRIALS.(trial_types{i}))'},'type','cper');
	RAW.mu.(trial_types{i})=mean(new_data(:,TRIALS.(trial_types{i}))');
end

new_data_regress=markolab_deltacoef(new_data',4,2)'; % approx 13 ms regression

for i=1:ntypes
	if isempty(TRIALS.(trial_types{i})), continue; end
	REGRESS.ci.(trial_types{i})=bootci(1e3,{@mean,new_data_regress(:,TRIALS.(trial_types{i}))'},'type','cper');
	REGRESS.mu.(trial_types{i})=mean(new_data_regress(:,TRIALS.(trial_types{i}))');
end

REGRESS.mat=new_data_regress;
REGRESS.t=time;

RAW.mat=new_data;
RAW.t=time;
