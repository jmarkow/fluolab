function [TRIAL_TIMES,CHANGE_POINTS,CHANGE_TRIALS]=fluolab_ttl_proc(TTL,TRIALS,varargin)
%
%
%
%
%
%
%


% loop through trials, find noise onset and changepoints (when we changed
% timing of playback)

ttl_thresh=2;
padding=[];

nparams=length(varargin);
if mod(nparams,2)>0
	error('ephysPipeline:argChk','Parameters must be specified as parameter/value pairs!');
end

for i=1:2:nparams
	switch lower(varargin{i})
		case 'ttl_thresh'
			ttl_thresh=varargin{i+1};
    case 'padding'
      padding=varargin{i+1};
  end
end

[nsamples,ntrials]=size(TTL.data);



pad_smps=round(padding.*TTL.fs);

if isempty(pad_smps) | pad_smps(1)==0
  pad_smps(1)=1;
end

if length(pad_smps)<2
  pad_smps(2)=0;
end

use_data=TTL.data(pad_smps(1):end-pad_smps(2),:);
nsamples=size(use_data,1);
idx=1:(nsamples-1);

TRIAL_TIMES=nan(1,ntrials);

for i=1:ntrials
    cur_trial=use_data(:,i)>ttl_thresh;

    onsets=~cur_trial(idx)&cur_trial(idx+1);
    onsets=min(find(onsets));
    if ~isempty(onsets)
      TRIAL_TIMES(i)=(onsets(1)+pad_smps(1))/TTL.fs;
    end
end

trial_idx=1:ntrials;

to_del=isnan(TRIAL_TIMES);

TRIAL_TIMES(to_del)=[];
trial_idx(to_del)=[];
TRIAL_TIMES=medfilt1(TRIAL_TIMES,5);

figure();plot(TRIAL_TIMES)
pause();

change=abs(diff([TRIAL_TIMES(1) TRIAL_TIMES]));
CHANGE_POINTS=trial_idx(find(change>mad(change,1)*100))+1;

% arrange trial structure to reflect changepoints (array of structs, copy everything)?


CHANGE_POINTS(CHANGE_POINTS>ntrials)=[];
CHANGE_POINTS=unique([1 CHANGE_POINTS ntrials]);
epochs=length(CHANGE_POINTS)-1;

conditions=fieldnames(TRIALS.all);
conditions(strcmp(conditions,'fluo_include'))=[];

for i=1:length(conditions)
  cur_trials=TRIALS.all.(conditions{i});
  fluo_trials=TRIALS.fluo_include.(conditions{i});
  map_trials=TRIALS.all.fluo_include(fluo_trials);

  for j=1:epochs
    CHANGE_TRIALS(j).all.(conditions{i})=cur_trials(cur_trials>CHANGE_POINTS(j)&cur_trials<CHANGE_POINTS(j+1));
    CHANGE_TRIALS(j).fluo_include.(conditions{i})=fluo_trials(map_trials>CHANGE_POINTS(j)&map_trials<CHANGE_POINTS(j+1));
  end
end
