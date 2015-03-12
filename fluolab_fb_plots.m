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
	end
end

if length(TRIALS.daf)>1
	isdaf=1;
end

% take first catch or all trial for AUDIO sample

sample=min(TRIALS.include(TRIALS.catch));
if isempty(sample)
	sample=min(TRIALS.all);
end

[nsamples,ntrials]=size(AUDIO.data);

blanking_idx=round(blanking(1)*AUDIO.fs):nsamples-round(blanking(2)*AUDIO.fs);
blanking_idx_ttl=round(blanking(1)*TTL.fs):nsamples-round(blanking(2)*TTL.fs);

[s,f,t]=zftftb_pretty_sonogram(AUDIO.data(blanking_idx,sample),...
	AUDIO.fs,'filtering',300,'clipping',-3,'len',80,'overlap',79,'zeropad',0);
f=f/1e3;

if isdaf
	nplots=4;
	ttl_plot=TTL.data(blanking_idx_ttl,TRIALS.include);

	ntrials=length(TRIALS.daf);
	step_size=win_size-win_overlap;
	steps=1:step_size:ntrials-win_size-1;

	RAW.sliding_mu.win_size=win_size;
	RAW.sliding_mu.win_overlap=win_overlap;
	RAW.sliding_mu.mat.daf=zeros(size(RAW.mat,1),length(steps));

	counter=1;

	for i=steps
		RAW.sliding_mu.mat.daf(:,counter)=mean(RAW.mat(:,TRIALS.daf(i:i+win_size))');
		counter=counter+1;
	end

else
	nplots=3;
end

idx=1;

ax=[];
FLUOFIGS.songalign=figure('paperpositionmode','auto','visible',visible);
ax(idx)=subplot(nplots,1,idx);
imagesc(t,f,s);axis xy;box off;
set(gca,'xtick',[]);
colormap(colors);

idx=idx+1;

if isdaf

	title(['ntrials catch: ' num2str(length(TRIALS.catch)) ' ntrials daf:  ' num2str(length(TRIALS.daf))]);
	
	ax(idx)=subplot(nplots,1,idx);

	% plot probability of feedback

	[nsamples,ntrials]=size(ttl_plot);
	plot([1:nsamples]/TTL.fs,mean(ttl_plot'>1),'r-','linewidth',2);
	set(gca,'xtick',[],'ytick',[0 1]);ylim([0 1]);
	box off;

	idx=idx+1;

	ax(idx)=subplot(nplots,1,idx:idx+1);
	markolab_shadeplot(RAW.t,RAW.ci.daf,facecolor_fb,edgecolor_fb);
	hold on;
	plot(RAW.t,RAW.mu.daf,'k-');

	markolab_shadeplot(RAW.t,RAW.ci.catch,facecolor,edgecolor);
	hold on;
	plot(RAW.t,RAW.mu.catch,'k-');

	obs_min=min([RAW.ci.catch(:);RAW.ci.daf(:)]);
	obs_max=max([RAW.ci.catch(:);RAW.ci.daf(:)]);

	alpha([facealpha]);

else

	ax(idx)=subplot(nplots,1,idx);

	title(['ntrials: ' num2str(length(TRIALS.all))]);
	
	markolab_shadeplot(RAW.t,RAW.ci.all,facecolor,edgecolor,1);
	hold on;
	plot(RAW.t,RAW.mu.all,'k-');
	
	obs_min=min([RAW.ci.all(:)]);
	obs_max=max([RAW.ci.all(:)]);

end

%

ylabel('Fluo. (normalized)');
xlabel('Time (s)');
box off;
linkaxes(ax,'x');

if obs_min<ylimits(1)
	ylimits(1)=obs_min;
end

if obs_max>ylimits(2)
	ylimits(2)=obs_max;
end

ylim([ylimits]);
set(gca,'YTick',ylimits);

% sliding window plots

ntrials=length(TRIALS.catch);
step_size=win_size-win_overlap;
steps=1:step_size:ntrials-win_size-1;
RAW.sliding_mu.mat.catch=zeros(size(RAW.mat,1),length(steps));

counter=1;

for i=steps
	RAW.sliding_mu.mat.catch(:,counter)=mean(RAW.mat(:,TRIALS.catch(i:i+win_size))');
	counter=counter+1;
end

%  two sliding window plots, then mic amplitude plots

ax=[];

FLUOFIGS.smoothtrials_catch=figure('paperpositionmode','auto','visible',visible);

ax(1)=subplot(3,1,1);
imagesc(t,f,s);axis xy;box off;
set(gca,'xtick',[]);
colormap(colors);

ax(2)=subplot(3,1,2:3);
imagesc(RAW.t,[],RAW.sliding_mu.mat.catch');
ylabel('Trials');
xlabel('Time(s)');
box off;
linkaxes(ax,'x');

pos=get(gca,'position');
clims=caxis();
clims=round(clims*100)/100;
caxis(clims);

hc = colorbar('location','eastoutside','position',...
	[pos(1)+pos(3)+cbar_dist pos(2) cbar_width pos(4)],'fontsize',10,'ytick',clims);


ax=[];

if isdaf

	FLUOFIGS.smoothtrials_daf=figure('paperpositionmode','auto','visible',visible);

	ax(1)=subplot(4,1,1);
	imagesc(t,f,s);axis xy;box off;
	set(gca,'xtick',[]);
	colormap(colors);


	[nsamples,~]=size(ttl_plot);

	ax(2)=subplot(4,1,2);
	plot([1:nsamples]/TTL.fs,mean(ttl_plot'>1),'r-','linewidth',2);
	set(gca,'xtick',[],'ytick',[0 1]);

	ax(3)=subplot(4,1,3:4);
	imagesc(RAW.t,[],RAW.sliding_mu.mat.daf');
	ylabel('Trials');
	xlabel('Time(s)');
	box off;
	linkaxes(ax,'x');

	clims=caxis();
	clims=round(clims*100)/100;
	caxis(clims);

	pos=get(gca,'position');
	hc = colorbar('location','eastoutside','position',...		
		[pos(1)+pos(3)+cbar_dist pos(2) cbar_width pos(4)],'fontsize',10,'ytick',clims);

end


