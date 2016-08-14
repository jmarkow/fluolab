function WIN_DATA=fluoflab_sliding_window(DATA,TRIAL_TIMES,CHANGE_IDX,varargin)
% plot window centered on ttl pulse
%
%
%
%

thresh=.5;
time_win=[.005 .4];
fs=200;

nparams=length(varargin);

if mod(nparams,2)>0
	error('ephysPipeline:argChk','Parameters must be specified as parameter/value pairs!');
end

for i=1:2:nparams
	switch lower(varargin{i})
    case 'time_win'
    	time_win=varargin{i+1};
		case 'fs'
			fs=varargin{i+1};
	end
end

[nsamples,ntrials]=size(DATA);

% for each session, get the average time to fill in for trials for with no TTL

time_smps=round(time_win.*fs);
time_vec=[-time_smps(1):time_smps(2)];
WIN_DATA=nan(length(time_vec),ntrials);

changes=unique(CHANGE_IDX);
nchanges=length(changes);
mu_times=zeros(nchanges,1);

for i=1:length(changes)
	mu_times(i)=nanmedian(TRIAL_TIMES(CHANGE_IDX==changes(i)));
end

for i=1:ntrials

	% if we have a nan (no trigger), use last good catch trial time_win

	cur_point=round(TRIAL_TIMES(i)*fs);
	cur_idx=CHANGE_IDX(i);

    if isnan(cur_point)
		cur_point=round(mu_times(cur_idx==changes)*fs);
    end

	if cur_point-time_smps(1)>0 & cur_point+time_smps(2)<=nsamples
		WIN_DATA(:,i)=DATA(cur_point-time_smps(1):cur_point+time_smps(2),i);
	end
end
