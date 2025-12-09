%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This code uses freely available software chunkie (installation required):
% https://chunkie.readthedocs.io/en/latest/getchunkie.html
% https://github.com/fastalgorithms/chunkie

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%          Create a geometry with a stadium cavity and a few ellipses inside
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
global_start_time=tic;
omega=2*pi;
k=omega;

%%%%%%%%%%%%%%%%%%      Define cavity, leads and ellipses parameters
%%%     Stadium curvature radius and cavity dimensions
geom_params.cava=2.0;
geom_params.cavb=10.0;
geom_params.stadr=geom_params.cavb/2;
%%%     Lengths, widths and y-positions of centers of all the leads
geom_params.numleads=3;
geom_params.leadsl=2.0*ones(1,geom_params.numleads);
geom_params.leadsd=0.6*ones(1,geom_params.numleads);
geom_params.leadscens=[3.0,0.1,-2.4];
%%%     Ellipses centers, semiaxes and angles
geom_params.numell=3;
geom_params.ellcens=[-0.5,-0.2,-0.75;
          1.0,0.5,-0.5];
geom_params.ellsemi=[0.1*ones(1,geom_params.numell);
          0.2*ones(1,geom_params.numell)];
geom_params.elltheta=[pi/4,pi/3,-pi/5];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                  Define all scattering properties of all leads and all modes in each lead
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       These are mostly internal objects - you don't need to change them
lead_params.numchannels=3;

%linspace(a/2+0.1*d1,a/2+0.9*d1,2E1)
lead_params.lead_x_terminations=geom_params.cava/2+geom_params.leadsl;
lead_params.lead_x_calcs=geom_params.cava/2+0.8*geom_params.leadsl;
lead_params.lead_y_calcs=zeros(geom_params.numleads,2E1);
for i2=1:geom_params.numleads
lead_params.lead_y_calcs(i2,:)=geom_params.leadscens(i2)+linspace(-0.4*geom_params.leadsd(i2),0.4*geom_params.leadsd(i2),size(lead_params.lead_y_calcs,2));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%          Get S-matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[S,overall_geom,kern,kern_der,all_channels_info,a_in_inv,sysmat,leads]=get_S_matrix_v15(omega,geom_params,lead_params);
S'*S
S*S'
S

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%          Example of geometry plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig=figure('units','normalized','position',[0,0,1,1]);
hold on
plot(overall_geom,'b-x')
quiver(overall_geom,'r')
hold off
axis equal tight

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%          Examples of field plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%   !WARNING!  %%%%%%%%%
%%%%%%%%%   Make sure you only evaluate the fields some distance away
%%%%%%%%%   from the boundary. If you try to evaluate too close to the
%%%%%%%%%   boundary, chunkie will use its nasty internal quadratures and
%%%%%%%%%   the evaluation will be very very slow, without telling you
%%%%%%%%%   anything
%%%%%%%%%   The boundary includes ellipses and lead ends

%%%   In these examples, we will set up an incoming excitation through one of
%%%   the channels and evaluate the fields in one of the other leads, and
%%%   then inside the cavity due to this excitation

%%  Field in one lead
% Evaluation points grid
x=linspace(geom_params.cava/2,geom_params.cava/2+0.9*geom_params.leadsl(1),5E1);
y=geom_params.leadscens(1)+linspace(-0.45*geom_params.leadsd(1),0.45*geom_params.leadsd(1),3E1);
[X,Y]=meshgrid(x,y);
XY=[X(:).';Y(:).'];

% We find the correct linear combination of solved sigma charge densities
% to make exactly one of the channels (inc_channel_ind) purely incoming
inc_channel_ind=3;
tic
c1=a_in_inv(1,inc_channel_ind);
c2=a_in_inv(2,inc_channel_ind);
c3=a_in_inv(3,inc_channel_ind);
special_sigma=c1*all_channels_info{1}.sigma+c2*all_channels_info{2}.sigma+c3*all_channels_info{3}.sigma;
% Call chunkie routine to evaluate field given geometry, kernel, charge
% density sigma and 2-by-N matrix of [x;y] evaluation points coordinates
psitot_chunkie=chunkerkerneval(overall_geom,kern,special_sigma,XY).';
psitot_eval_time=toc;
% Notice that this evaluates the total field

% Plot
%thingtoplot=abs( reshape(psiscat_chunkie,size(X)) );
%thingtoplot=reshape(thingtoplot/max(thingtoplot,[],'all'),[],1);
%colors=[reshape(thingtoplot/max(thingtoplot,[],'all'),[],1),zeros(size(X,1)*size(X,2),2)];
tic
fig=figure('units','normalized','position',[0,0,1,1]);
hold on
colorplot=imagesc(   x,y,abs(  reshape(psitot_chunkie,size(X))  )   );
plot(overall_geom,'k')
%scatter(X,Y,40,'black','filled')
%scatter(reshape(X,[],1),reshape(Y,[],1),30,colors,'filled')
hold off
colorbar
set(gca,'ydir','normal')
colormap(redblue);
daspect([1,1,1])
plot_time=toc;

%%  Field inside the cavity
% Evaluation points grid
%x=linspace(geom_params.cava/2,geom_params.cava/2+0.9*geom_params.leadsl(1),5E1);
%y=geom_params.leadscens(1)+linspace(-0.45*geom_params.leadsd(1),0.45*geom_params.leadsd(1),3E1);
bdr_delta=3E-1;
x=linspace(-geom_params.stadr-geom_params.cava/2+bdr_delta,geom_params.cava/2-bdr_delta,5E1);
y=linspace(-geom_params.cavb/2+bdr_delta,geom_params.cavb/2-bdr_delta,5E1);
[X,Y]=meshgrid(x,y);
XX=X(:).';
YY=Y(:).';

within_cavity=( ((XX+geom_params.cava/2).^2+YY.^2)<=(geom_params.stadr-bdr_delta)^2 );
XX(~within_cavity)=[];
YY(~within_cavity)=[];

for i2=1:geom_params.numell
x0=geom_params.ellcens(1,i2);
y0=geom_params.ellcens(2,i2);
lx=geom_params.ellsemi(1,i2);
ly=geom_params.ellsemi(2,i2);
bdr_delta_now=bdr_delta/max([lx,ly]);
within_ellips=( ((XX-x0).^2/lx^2+(YY-y0).^2/ly^2)<=(1+bdr_delta_now)^2 );
XX(within_ellips)=[];
YY(within_ellips)=[];
end
XY=[XX(:).';YY(:).'];

% We find the correct linear combination of solved sigma charge densities
% to make exactly one of the channels (inc_channel_ind) purely incoming
inc_channel_ind=3;
tic
c1=a_in_inv(1,inc_channel_ind);
c2=a_in_inv(2,inc_channel_ind);
c3=a_in_inv(3,inc_channel_ind);
special_sigma=c1*all_channels_info{1}.sigma+c2*all_channels_info{2}.sigma+c3*all_channels_info{3}.sigma;
% Call chunkie routine to evaluate field given geometry, kernel, charge
% density sigma and 2-by-N matrix of [x;y] evaluation points coordinates
psitot_chunkie=chunkerkerneval(overall_geom,kern,special_sigma,XY).';
psitot_eval_time=toc;
% Notice that this evaluates the total field

% Plot real part
thingtoplot=real(psitot_chunkie);
redd=interp1([-max(abs(thingtoplot)),0,max(abs(thingtoplot))],[0,1,1],thingtoplot);
bluee=interp1([-max(abs(thingtoplot)),0,max(abs(thingtoplot))],[0,1,0],thingtoplot);
greenn=interp1([-max(abs(thingtoplot)),0,max(abs(thingtoplot))],[1,1,0],thingtoplot);
colors=[redd(:),bluee(:),greenn(:)];
colorbar_colors=[interp1([-max(abs(thingtoplot)),0,max(abs(thingtoplot))],[0,1,1],linspace(-max(abs(thingtoplot)),max(abs(thingtoplot)),10)).',...
    interp1([-max(abs(thingtoplot)),0,max(abs(thingtoplot))],[0,1,0],linspace(-max(abs(thingtoplot)),max(abs(thingtoplot)),10)).',...
    interp1([-max(abs(thingtoplot)),0,max(abs(thingtoplot))],[1,1,0],linspace(-max(abs(thingtoplot)),max(abs(thingtoplot)),10)).'];
pointsize=50;
tic
fig=figure('units','normalized','position',[0,0,1,1]);
hold on
plot(overall_geom,'k')
scatter(reshape(XX,[],1),reshape(YY,[],1),pointsize,colors,'filled')
hold off
colormap(colorbar_colors)
cb=colorbar;
cb.Ticks=[0,0.5,1];
cb.TickLabels=[-max(abs(thingtoplot)),0,max(abs(thingtoplot))];
set(gca,'ydir','normal')
daspect([1,1,1])
plot_time=toc;

%%
% Plot abs value
thingtoplot=abs(psitot_chunkie);
redd=interp1([0,max(abs(thingtoplot))],[1,1],thingtoplot);
bluee=interp1([0,max(abs(thingtoplot))],[1,0],thingtoplot);
greenn=interp1([0,max(abs(thingtoplot))],[1,0],thingtoplot);
colors=[redd(:),bluee(:),greenn(:)];
colorbar_colors=[interp1([0,max(abs(thingtoplot))],[1,1],linspace(0,max(abs(thingtoplot)),10)).',...
    interp1([0,max(abs(thingtoplot))],[1,0],linspace(0,max(abs(thingtoplot)),10)).',...
    interp1([0,max(abs(thingtoplot))],[1,0],linspace(0,max(abs(thingtoplot)),10)).'];
pointsize=50;
tic
fig=figure('units','normalized','position',[0,0,1,1]);
hold on
plot(overall_geom,'k')
%scatter(XX,YY,40,'black','filled')
scatter(reshape(XX,[],1),reshape(YY,[],1),pointsize,colors,'filled')
hold off
colormap(colorbar_colors)
cb=colorbar;
cb.Ticks=[0,1];
cb.TickLabels=[0,max(abs(thingtoplot))];
set(gca,'ydir','normal')
daspect([1,1,1])
plot_time=toc;
%}
%}