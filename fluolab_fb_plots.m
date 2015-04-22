function FLUOFIGS=fluolab_fb_plots(AUDIO,RAW,TTL,TRIALS,varargin)

blanking=[.2 .2];
facecolor=[0 0 1];
edgecolor=[0 0 1];
facecolor_fb=[1 0 0];
edgecolor_fb=[1 0 0];
ylimits=[.2 .7];
facealpha=.5;
win_size=20;
win_overlap=19;
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

conditions.names=fieldnames(TRIALS.fluo_include);
nconditions=length(conditions.names);
conditions.ntrials=zeros(1,nconditions);
conditions.proc=zeros(1,nconditions);

for i=1:nconditions
	conditions.ntrials(i)=length(TRIALS.fluo_include.(conditions.names{i}));
	conditions.proc(i)=conditions.ntrials(i)>1;
end

% take first catch or all trial for AUDIO sample

sample=min(TRIALS.all.catch);
if isempty(sample)
	sample=min(TRIALS.all.other);
	if isempty(sample)
		sample=1;
	end
end

[nsamples,ntrials]=size(AUDIO.data);
[s,f,t]=zftftb_pretty_sonogram(AUDIO.data(:,sample),...
	AUDIO.fs,'filtering',300,'clipping',[-2 2],'len',80,'overlap',79,'zeropad',0,'norm_amp',1);
f=f/1e3;

% do we have any feedback trials?

isdaf=conditions.proc(strcmp(conditions.names,'daf'));
iscatch=conditions.proc(strcmp(conditions.names,'catch'));
isother=conditions.proc(strcmp(conditions.names,'other'));

if isdaf & ~isempty(TTL)
	if isother
		TTL.plot=TTL.data(:,sort([TRIALS.all.daf(:);TRIALS.all.other(:)]));
	else
		TTL.plot=TTL.data(:,TRIALS.all.daf(:));
	end
else
	TTL=[];
end

idx=~strcmp(conditions.names,'daf');
other_names=conditions.names(idx)
other_proc=conditions.proc(idx);
other_names(other_proc==0)=[];
other_names(strcmp(other_names,'all'))=[];

if isdaf
	for i=1:length(other_names)
		fig_name=['songalign_daf_' other_names{i}];
		FLUOFIGS.(fig_name)=figure('paperpositionmode','auto','visible',visible);

		label=other_names{i};
		label(label=='_')='&';

		ax=fluolab_dualshade_plot(t,f,s,RAW.t,RAW.ci.daf,RAW.mu.daf,RAW.ci.(other_names{i}),RAW.mu.(other_names{i}),TTL,'labels',{'DAF',label});

		if ~isempty(ax)
			axes(ax(1));
			title(['ntrials catch: ' num2str(length(TRIALS.fluo_include.(other_names{i}))) ...
			   	' ntrials daf:  ' num2str(length(TRIALS.fluo_include.daf))]);
		end
	end
end

if iscatch & isother
	FLUOFIGS.songalign_other_catch=figure('paperpositionmode','auto','visible',visible);
	ax=fluolab_dualshade_plot(t,f,s,RAW.t,RAW.ci.other,RAW.mu.other,RAW.ci.catch,RAW.mu.catch,TTL,'labels',{'Other','Catch'});
	if ~isempty(ax)
		axes(ax(1));
		title(['ntrials other ' num2str(length(TRIALS.fluo_include.other)) ' ntrials catch:  ' num2str(length(TRIALS.fluo_include.catch))]);
	end
end

%  two sliding window plots, then mic amplitude plots


ntrials=size(AUDIO.data,2);
escape_idx=zeros(ntrials,1);
escape_idx(TRIALS.all.other)=1;
escape_idx=markolab_smooth(escape_idx,min(ntrials,100),'z');

idx=~strcmp(conditions.names,'all');
other_names=conditions.names(idx);
other_proc=conditions.proc(idx);
other_names(other_proc==0)=[]

FLUOFIGS

for i=1:length(other_names)

	fig_name=['smoothtrials_' other_names{i}];
	FLUOFIGS.(fig_name)=figure('paperpositionmode','auto','visible',visible);

	trials=TRIALS.all.fluo_include(TRIALS.fluo_include.(other_names{i}));

	tmp=[];
	tmp2=[];

	if ~isempty(datenums)
		tmp=datenums(trials);
	end

	if ~isempty(escape_idx)
		tmp2=escape_idx(trials);
	end

	trials=TRIALS.fluo_include.(other_names{i});

	label=other_names{i};
	label(label=='_')='&';

	fluolab_slidingwindow_plot(t,f,s,RAW.mat(:,trials),RAW.t,win_size,win_overlap,'datenums',tmp,'label',label,'escape_idx',tmp2);

end

FLUOFIGS.songalign_all=figure('paperpositionmode','auto','visible',visible);
ax=fluolab_dualshade_plot(t,f,s,RAW.t,RAW.ci.all,RAW.mu.all);

