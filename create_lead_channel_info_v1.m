function [leads,all_channels_info] = create_lead_channel_info_v1(omega,boundary_xy,numleads,numchannels,...
                                                                 leadsl,leadsd,leadscens,lead_x_terminations,lead_x_calcs,lead_y_calcs)

leads=cell(1,numleads);
all_channels_info=cell(1,numchannels);
for i2=1:numleads
all_channels_info{i2}.x_termination=lead_x_terminations(i2);
all_channels_info{i2}.leadend_inds=( abs( boundary_xy(1,:)-all_channels_info{i2}.x_termination )<1E-5 & ...
                                     abs( boundary_xy(2,:)-leadscens(i2) )<=leadsd(i2)/2 );
all_channels_info{i2}.leadend_xy=boundary_xy(2,all_channels_info{i2}.leadend_inds);
all_channels_info{i2}.mode_fun=@(y) sin(pi/2+pi*(y-leadscens(i2))/leadsd(i2));
all_channels_info{i2}.beta=sqrt(omega^2-pi^2/leadsd(i2)^2);
all_channels_info{i2}.save_rhs_plot=false;
end
for i2=1:numleads
leads{i2}.x=lead_x_calcs(i2);
leads{i2}.y=lead_y_calcs(i2,:);
[X,Y]=meshgrid(leads{i2}.x,leads{i2}.y);
leads{i2}.XY=[X(:).';Y(:).'];
leads{i2}.XY_ptinfo.r=leads{i2}.XY;
leads{i2}.XY_ptinfo.d=repmat([0,1],length(leads{i2}.y)).';
leads{i2}.XY_ptinfo.d2=zeros(length(leads{i2}.y),2).';
leads{i2}.XY_ptinfo.n=repmat([-1,0],length(leads{i2}.y)).';
beta=sqrt(omega^2-pi^2/leadsd(i2)^2);
leads{i2}.norm_factor=sqrt( beta*leadsd(i2)/2/omega );
leads{i2}.out_mode_matrix=1/leads{i2}.norm_factor*sin(pi/2+pi*(leads{i2}.y-leadscens(i2))/leadsd(i2)).';
%leads{i2}.save_lsqr_plots=true;
leads{i2}.save_lsqr_plots=false;
end

end
