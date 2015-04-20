function FLUOFIGS=fluolab_fb_plots(AUDIO,RAW,TTL,TRIALS,varargin)

blanking=[.2 .2];
facecolor=[0 0 1];
edgecolor=[0 0 1];
facecolor_fb=[1 0 0];
edgecolor_fb=[1 0 0];
ylimits=[.2 .7];
facealpha=.5;
win_size=30;
win_overlap=29;
colors='jet';
cbar_width = .025; 
cbar_dist=.01;
visible='on';
datenums=[];

nparams=length(varargin);

if mod(nparams,2)>0
	error('ephysPipeline:argChk','Parameters must be specified as parameter/value pairs!');
end

for i=1:2:nparams
	switch lower(varargin{i})
		case 'blanking'
			blanking=varargin{i+1};
		case 'facecolor'
			facecolor=varargin{i+1};
		case 'edgecolor'
			edgecolor=varargin{i+1};
		case 'edgecolor_fb'
			edgecolor_fb=varargin{i+1};
		case 'facecolor_fb'
			facecolor_fb=varargin{i+1};
		case 'ylimits'
			ylimits=varargin{i+1};
		case 'facealpha'
			facealpha=varargin{i+1};
		case 'win_size'
			win_size=varargin{i+1};
		case 'win_overlap'
			win_overlap=varargin{i+1};
		case 'visible'
			visible=varargin{i+1};
		case 'datenums'
			datenums=varargin{i+1};
	end
end

isdaf=length(TRIALS.fluo_include.daf)>1;
iscatch=length(TRIALS.fluo_include.catch)>1;
isother=length(TRIALS.fluo_include.other)>1;
iscatchandother=length(TRIALS.fluo_include.catch_other)>1;

% take first catch or all trial for AUDIO sample

sample=min(TRIALS.all.catch);
if isempty(sample)
	sample=min(TRIALS.all.other);
	if isempty(sample)
		sample=1;
	end
end

[nsamples,ntrials]=size(AUDIO.data);

blanking_idx=round(blanking(1)*AUDIO.fs):nsamples-round(blanking(2)*AUDIO.fs);
blanking_idx_ttl=round(blanking(1)*TTL.fs):nsamples-round(blanking(2)*TTL.fs);

[s,f,t]=zftftb_pretty_sonogram(AUDIO.data(:,sample),...
	AUDIO.fs,'filtering',300,'clipping',[-2 2],'len',80,'overlap',79,'zeropad',0,'norm_amp',1);
f=f/1e3;

if isdaf
	if isother
		TTL.plot=TTL.data(:,sort([TRIALS.all.daf(:);TRIALS.all.other(:)]));
	else
		TTL.plot=TTL.data(:,TRIALS.all.daf(:));
	end
else
	TTL=[];
end

if isdaf & iscatch
	FLUOFIGS.songalign_daf_catch=figure('paperpositionmode','auto','visible',visible);
	ax=fluolab_dualshade_plot(t,f,s,RAW.t,RAW.ci.daf,RAW.mu.daf,RAW.ci.catch,RAW.mu.catch,TTL,'labels',{'DAF','Catch'});
	axes(ax(1));
	title(['ntrials catch: ' num2str(length(TRIALS.fluo_include.catch)) ' ntrials daf:  ' num2str(length(TRIALS.fluo_include.daf))]);
end


if isdaf & isother
	FLUOFIGS.songalign_daf_other=figure('paperpositionmode','auto','visible',visible);
	ax=fluolab_dualshade_plot(t,f,s,RAW.t,RAW.ci.daf,RAW.mu.daf,RAW.ci.other,RAW.mu.other,TTL,'labels',{'DAF','Other'});
	axes(ax(1));
	title(['ntrials other ' num2str(length(TRIALS.fluo_include.other)) ' ntrials daf:  ' num2str(length(TRIALS.fluo_include.daf))]);
end

if iscatch & isother
	FLUOFIGS.songalign_other_catch=figure('paperpositionmode','auto','visible',visible);
	ax=fluolab_dualshade_plot(t,f,s,RAW.t,RAW.ci.other,RAW.mu.other,RAW.ci.catch,RAW.mu.catch,TTL,'labels',{'Other','Catch'});
	axes(ax(1));
	title(['ntrials other ' num2str(length(TRIALS.fluo_include.other)) ' ntrials catch:  ' num2str(length(TRIALS.fluo_include.catch))]);
end

if isdaf & iscatchandother
	FLUOFIGS.songalign_daf_catchandother=figure('paperpositionmode','auto','visible',visible);
	ax=fluolab_dualshade_plot(t,f,s,RAW.t,RAW.ci.daf,RAW.mu.daf,RAW.ci.catch_other,RAW.mu.catch_other,TTL,'labels',{'DAF','Catch & Other'});
	axes(ax(1));
	title(['ntrials catch and other ' num2str(length(TRIALS.fluo_include.catch_other)) ' ntrials daf:  ' num2str(length(TRIALS.fluo_include.daf))]);
end

%  two sliding window plots, then mic amplitude plots


if iscatch
	% sliding window plots

	FLUOFIGS.smoothtrials_catch=figure('paperpositionmode','auto','visible',visible);

	if ~isempty(datenums)
		tmp=datenums(TRIALS.all.fluo_include(TRIALS.fluo_include.catch));
	else
		tmp=[];
	end

	fluolab_slidingwindow_plot(t,f,s,RAW.mat(:,TRIALS.fluo_include.catch),RAW.t,win_size,win_overlap,'datenums',tmp,'label','Catch');
end


if isdaf
	FLUOFIGS.smoothtrials_daf=figure('paperpositionmode','auto','visible',visible);

	if ~isempty(datenums)
		tmp=datenums(TRIALS.all.fluo_include(TRIALS.fluo_include.daf));
	else
		tmp=[];
	end


	fluolab_slidingwindow_plot(t,f,s,RAW.mat(:,TRIALS.fluo_include.daf),RAW.t,win_size,win_overlap,'datenums',tmp,'label','DAF');
end

if isother
	FLUOFIGS.smoothtrials_other=figure('paperpositionmode','auto','visible',visible);

	if ~isempty(datenums)
		tmp=datenums(TRIALS.all.fluo_include(TRIALS.fluo_include.other));
	else
		tmp=[];
	end


	fluolab_slidingwindow_plot(t,f,s,RAW.mat(:,TRIALS.fluo_include.other),RAW.t,win_size,win_overlap,'datenums',tmp,'label','Other');
end

if iscatchandother
	FLUOFIGS.smoothtrials_catchandother=figure('paperpositionmode','auto','visible',visible);

	if ~isempty(datenums)
		tmp=datenums(TRIALS.all.fluo_include(TRIALS.fluo_include.catch_other));
	else
		tmp=[];
	end


	fluolab_slidingwindow_plot(t,f,s,RAW.mat(:,TRIALS.fluo_include.catch_other),RAW.t,win_size,win_overlap,'datenums',tmp,'label','Catch & other');
end

FLUOFIGS.songalign_all=figure('paperpositionmode','auto','visible',visible);
ax=fluolab_dualshade_plot(t,f,s,RAW.t,RAW.ci.all,RAW.mu.all);

