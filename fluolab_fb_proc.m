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
normalize=1;
newfs=100;
dff=1;
tau=.1;
detrend_win=.3;
classify_trials=1;

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

%%%%

blanking_idx=[round(blanking(1)*DATA.fs):nsamples-round(blanking(2)*DATA.fs)];

% include these trials

[~,bad_trial]=find(proc_data(blanking_idx,:)<trial_cut);
include_trials=setdiff(1:ntrials,unique(bad_trial));
ntrials=length(include_trials);

%pause();

% where are the feedback trials?

if classify_trials==1

	C=fluolab_classify_trials(TTL.data(blanking_idx,include_trials),TTL.fs);
	TRIALS.catch=find(C==2|C==0);
	TRIALS.daf=find(C==1);

else

	[b,a]=ellip(4,.2,40,[10e3 12e3]/(AUDIO.fs/2),'bandpass');
	daf_mat=filtfilt(b,a,double(AUDIO.data(blanking_idx,include_trials))).^2>daf_level;
	[~,daf_trials]=find(daf_mat);
	daf_trials=unique(daf_trials);
	first_daf=1;

	% for now power detection 10e3-12e3 works fine

	% catch trials
	TRIALS.catch=setdiff(1:ntrials,daf_trials);
	TRIALS.daf=daf_trials;

end

TRIALS.all=1:ntrials;
TRIALS.include=include_trials;

trial_types=fieldnames(TRIALS);
ntypes=length(trial_types);

[new_data,time]=fluolab_condition(proc_data(:,include_trials),DATA.fs,'tau',tau,'detrend_win',detrend_win,...
	'newfs',newfs,'blanking',blanking,'normalize',normalize,'dff',dff);
[nsamples,ntrials]=size(new_data);
new_data_regress=markolab_deltacoef(new_data',4,2)'; % approx 13 ms regression

for i=1:ntypes
	if isempty(TRIALS.(trial_types{i})), continue; end
	if strcmp(trial_types{i},'include'), continue; end

	RAW.ci.(trial_types{i})=bootci(1e3,{@mean,new_data(:,TRIALS.(trial_types{i}))'},'type','cper');
	RAW.mu.(trial_types{i})=mean(new_data(:,TRIALS.(trial_types{i}))');
	
	REGRESS.ci.(trial_types{i})=bootci(1e3,{@mean,new_data_regress(:,TRIALS.(trial_types{i}))'},'type','cper');
	REGRESS.mu.(trial_types{i})=mean(new_data_regress(:,TRIALS.(trial_types{i}))');

end

REGRESS.mat=new_data_regress;
REGRESS.t=time;
REGRESS.fs=newfs;

RAW.mat=new_data;
RAW.t=time;
RAW.fs=newfs;
