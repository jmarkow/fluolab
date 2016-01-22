%%% user parameters

trial_cut=1;
nmads=0;
tau=0.1;
tau_regress=.015;
smooth_type='b';
channel=1:2;
newfs=200;

padding=[ .4 .25 ];
detrend_win=.2;
detrend_method='p';
normalize='n';

classify_trials='t';
daf_level=0;

% clean up the data

fluo_data=fluolab_datascrub(adc,'channel',channel,'trial_cut',trial_cut,'nmads',nmads);

% where are the feedback trials?

trials=fluolab_classify_trials(ttl,audio,...
		'include_trials',fluo_data.trial_idx,'method',classify_trials,'daf_level',daf_level,...
		'padding',padding);

% convert to df/f_0

[fluo.mat,fluo.t]=fluolab_condition(fluo_data.data(:,:,channel),fluo_data.fs,fluo_data.t,'tau',tau,'detrend_win',detrend_win,...
	'newfs',newfs,'normalize',normalize,'detrend_method',detrend_method,'smooth_type',smooth_type);
fluo.fs=newfs
fluo_regress=fluo;

% regression coefficients (smooth deriv.)

fluo_regress.mat=markolab_deltacoef(fluo.mat',round(tau_regress*newfs),2)';

% now align to interesting events

[trial_times change_points change_trials]=fluolab_ttl_proc(ttl,trials,'padding',padding);

% plot results...
