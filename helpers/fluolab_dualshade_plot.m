function AX=fluolab_dualshade_plot(T,F,S,CI_T,CI1,MU1,CI2,MU2,TTL,varargin)
%
%
%

if nargin<9, TTL=[]; end
if nargin<8, MU2=[]; end
if nargin<7, CI2=[]; end

facecolor=[0 0 1];
edgecolor=[0 0 1];
facecolor_fb=[1 0 0];
edgecolor_fb=[1 0 0];
datenums=[];
%ylimits=[.2 .7];
ylimits=[];
ylim_order=100;
facealpha=.5;
colors='jet';
labels={};

AX=[];

nparams=length(varargin);

if mod(nparams,2)>0
	error('ephysPipeline:argChk','Parameters must be specified as parameter/value pairs!');
end

for i=1:2:nparams
	switch lower(varargin{i})
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
		case 'cbar_dist'
			cbar_dist=varargin{i+1};
		case 'labels'
			labels=varargin{i+1};
	end
end

nsamples=length(CI_T);
nsamples2=size(CI1,2);
nsamples3=size(CI2,2);

if nsamples~=nsamples2 | (~isempty(CI2)&nsamples~=nsamples3)
	warning('Time vector and confidence interval vectors do not match, skipping...');
	return;
end


nplots=3;

if ~isempty(TTL)
	nplots=nplots+1;
end

idx=1;
AX(end+1)=subplot(nplots,1,1);
imagesc(T,F,S);axis xy;box off;
set(gca,'xtick',[]);
colormap(colors);

idx=idx+1;

if ~isempty(TTL)
	% plot probability of feedback
	
	AX(end+1)=subplot(nplots,1,idx);
	[nsamples,ntrials]=size(TTL.plot);
	plot([1:nsamples]/TTL.fs,mean(TTL.plot'>1),'r-','linewidth',2);
	ylabel('P(feedback)');
	set(gca,'xtick',[],'ytick',[0 1],'TickLength',[0 0]);ylim([0 1]);
	box off;

	idx=idx+1;

end

h=[];

AX(end+1)=subplot(nplots,1,idx:idx+1);
h(1)=markolab_shadeplot(CI_T,CI1,facecolor_fb,edgecolor_fb);
hold on;
plot(CI_T,MU1,'k-');

if ~isempty(CI2)
	h(2)=markolab_shadeplot(CI_T,CI2,facecolor,edgecolor);
	hold on;
	plot(CI_T,MU2,'k-');

	obs_min=min([CI1(:);CI2(:)]);
	obs_max=max([CI1(:);CI2(:)]);
	alpha([facealpha]);

	if ~isempty(labels)
		L=legend(h,labels);
		legend boxoff;
	end
else
	obs_min=min(CI1(:));
	obs_max=max(CI1(:));
end

idx=idx+2;

ylabel('Fluo. (normalized)');
xlabel('Time (s)');
box off;
set(gca,'TickLength',[0 0]);
linkaxes(AX,'x');
xlim([max(CI_T(1),T(1)) min(CI_T(end),T(end))]);

new_ylimits=ylimits;

if isempty(ylimits) | obs_min<new_ylimits(1)
	new_ylimits(1)=floor(obs_min*ylim_order)/ylim_order;
end

if isempty(ylimits) | obs_max>new_ylimits(2)
	new_ylimits(2)=ceil(obs_max*ylim_order)/ylim_order;
end

ylim([new_ylimits]);
set(gca,'YTick',new_ylimits);


