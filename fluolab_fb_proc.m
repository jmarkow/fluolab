function [RAW,REGRESS,TRIALS]=fluolab_fb_proc(DATA,AUDIO,TTL,varargin)
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

blanking=[.05 .05];
channel=1;
daf_level=.1;
trial_cut=2;
normalize='m';
newfs=100;
dff=1;
tau=.1;
detrend_win=.3;
classify_trials='t';

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
	end
end

proc_data=double(DATA.data(:,:,channel));
[nsamples,ntrials]=size(proc_data);

% include these trials

blanking_idx=[round(blanking(1)*DATA.fs):nsamples-round(blanking(2)*DATA.fs)];

[~,bad_trial]=find(proc_data(blanking_idx,:)<trial_cut);
include_trials=setdiff(1:ntrials,unique(bad_trial));

%pause();

% where are the feedback trials?

[TRIALS,C]=fluolab_classify_trials(TTL,AUDIO,'include_trials',include_trials,'method',classify_trials,'blanking',blanking);
ntrials=length(include_trials);

TRIALS.fluo_include.catch=find(C(include_trials)==2);
TRIALS.fluo_include.daf=find(C(include_trials)==1);
TRIALS.fluo_include.other=find(C(include_trials)==0);
TRIALS.fluo_include.catch_other=find(C(include_trials)==0|C(include_trials)==2);

TRIALS.fluo_include.all=[1:length(include_trials)];
TRIALS.all.fluo_include=include_trials;

trial_types=fieldnames(TRIALS.fluo_include);
ntypes=length(trial_types);

[new_data,time]=fluolab_condition(proc_data(blanking_idx,include_trials),DATA.fs,blanking_idx/DATA.fs,'tau',tau,'detrend_win',detrend_win,...
	'newfs',newfs,'normalize',normalize,'dff',dff);
[nsamples,ntrials]=size(new_data);
new_data_regress=markolab_deltacoef(new_data',4,2)'; % approx 13 ms regression

for i=1:ntypes

	if isempty(TRIALS.fluo_include.(trial_types{i})), continue; end
	if strcmp(trial_types{i},'include'), continue; end
	if strcmp(trial_types{i},'idx'), continue; end

	RAW.ci.(trial_types{i})=bootci(1e3,{@mean,new_data(:,TRIALS.fluo_include.(trial_types{i}))'},'type','cper');
	RAW.mu.(trial_types{i})=mean(new_data(:,TRIALS.fluo_include.(trial_types{i}))');

	REGRESS.ci.(trial_types{i})=bootci(1e3,{@mean,new_data_regress(:,TRIALS.fluo_include.(trial_types{i}))'},'type','cper');
	REGRESS.mu.(trial_types{i})=mean(new_data_regress(:,TRIALS.fluo_include.(trial_types{i}))');

end

REGRESS.mat=new_data_regress;
REGRESS.t=time;
REGRESS.fs=newfs;

RAW.mat=new_data;
RAW.t=time;
RAW.original_t=blanking_idx/DATA.fs;
RAW.fs=newfs;
