%%% user parameters
%%
trial_cut=1;
nmads=5;
tau=0.05;
tau_regress=.015;
smooth_type='b';
channel=1;
newfs=200;

padding=[ 0 0 ];
detrend_win=.6;
detrend_method='p';
normalize='n';

classify_trials='t';
daf_level=0;

% clean up the data
%%
fluo_data=fluolab_datascrub(adc,'channel',channel,'trial_cut',trial_cut,'nmads',nmads);

% where are the feedback trials?
%%
[trials,trials_idx]=fluolab_classify_trials(ttl,audio,...
		'include_trials',fluo_data.trial_idx,'method',classify_trials,'daf_level',daf_level,...
		'padding',template.extract_options.padding);
trials_idx_fluo=trials_idx(trials.all.fluo_include);

%%
% convert to df/f_0

[fluo.mat,fluo.t]=fluolab_condition(fluo_data.data(:,:,channel),fluo_data.fs,fluo_data.t,'tau',tau,'detrend_win',detrend_win,...
	'newfs',newfs,'normalize',normalize,'detrend_method',detrend_method,'smooth_type',smooth_type);
fluo.fs=newfs
fluo_regress=fluo;

% regression coefficients (smooth deriv.)

% normalize trials by looking outside of singing???

fluo_regress.mat=markolab_deltacoef(fluo.mat',round(tau_regress*newfs),2)';

% now align to interesting events
%%
[trial_times change_points change_trials,change_idx]=fluolab_ttl_proc(ttl,trials,'padding',template.extract_options.padding);
trial_times_fluo=trial_times(trials.all.fluo_include);
change_idx_fluo=change_idx(trials.all.fluo_include);
% align to trial_times window, plot peak time? value? integrated signal?

% this should be very straightforward
%%

win_data=fluolab_window_data(zscore(fluo.mat),fluo.t,trial_times_fluo,change_idx_fluo,'fs',fluo.fs);