function [S,overall_geom,kern,kern_der,all_channels_info,a_in_inv,sysmat,leads] = get_S_matrix_v15(omega,geom_params,lead_params)
% This function creates chunkie geometry, kernels and forms system matrix based on user input
% The function then makes N system solves for N right hand sides to get
% S-matrix.
% This code uses freely available software chunkie (installation required):
% https://chunkie.readthedocs.io/en/latest/getchunkie.html
% https://github.com/fastalgorithms/chunkie
% 
% Syntax: [S,optional_outputs] = get_S_matrix_v15(omega,geom_params,lead_params)
% 
% Input:
%   omega - frequency
%
%   geom_params - structure containing geometric parameters. The following
%       is the current geometric structure this is set up for
%       geom_params.cava - length of the straight x-side of the cavity
%       geom_params.cavb - length of the straight y-side of the cavity
%       geom_params.stadr - radius of the stadium arc
%       geom_params.numleads - number of leads
%       geom_params.leadsl - lengths of leads (1-by-numleads vector)
%       geom_params.leadsd - widths of leads (1-by-numleads vector)
%       geom_params.leadscens - y-coordinates of centers of leads (1-by-numleads vector)
%       geom_params.numell - number of ellipses inside the cavity
%       geom_params.ellcens - coordinates of ellipses' centers (2-by-ellsemi matrix [x;y])
%       geom_params.ellsemi - ellipses' semiaxes (1-by-ellsemi vector)
%       geom_params.elltheta - ellipses' angles of rotation (1-by-ellsemi vector)
% 
%   lead_params - structure containing lead parameters used in S-matrix
%   calculation
%       lead_params.numchannels - total number of channels across all
%       leads; right now, only lead_params.numchannels=geom_params.numleads
%       is supported
%       lead_params.lead_x_terminations - x-coordinates of ends of all
%       leads (1-by-geom_params.numleads vector)
%       lead_params.lead_x_calcs - x-coordinates of line cuts where field
%       overlaps will be done for S-matrix evaluation (1-by-geom_params.numleads vector)
%       lead_params.lead_y_calcs - y-coordinates of points where field
%       overlaps will be done for S-matrix evaluation (geom_params.numleads-by-user_specified_number matrix)
% 
% Output:
%   S - S-matrix (lead_params.numchannels-by-lead_params.numchannels)
%
% Optional outputs
%   overall_geom - chunkie geometry object
%   kern - chunkie kernel object (required for field plots)
%   kern_der - chunkie kernel derivative object
%   all_channels_info - channel info structure (contains solved charge densities from S-matrix calculation; required for field plots)
%   a_in_inv - matrix arising in the S-matrix calculation (required for field plots)
%   sysmat - system matrix
%   leads - leads info structure


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%          Create geometry
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Calls a routine which returns a chunkie geometry object, parametrizing
%   all the boundaries in the problem
overall_geom=create_stadiumgeom_v1(geom_params.cava,geom_params.cavb,geom_params.stadr,geom_params.numleads,...
    geom_params.leadsl,geom_params.leadsd,geom_params.leadscens,geom_params.numell,geom_params.ellcens,geom_params.ellsemi,geom_params.elltheta);
boundary_xy=reshape(overall_geom.r,2,[]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%          Define all lead and channel info
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Calls a routine which creates leads and channels structures used in
%   S-matrix calculation (these have nothing to do with chunkie - this is my personal book keeping)
numchannels=lead_params.numchannels;
numleads=geom_params.numleads;
[leads,all_channels_info] = create_lead_channel_info_v1(omega,boundary_xy,geom_params.numleads,lead_params.numchannels,...
                                                                 geom_params.leadsl,geom_params.leadsd,geom_params.leadscens,...
                                                                 lead_params.lead_x_terminations,lead_params.lead_x_calcs,lead_params.lead_y_calcs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%          Set integral equation kernel and create discretized system matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Calls two chunkie routines which initialize a chunkie kernel object and
%   kernel derivative object
%tic
kern=kernel('h','c',omega,[1,-1i*omega]);
kern_der=kernel('h','cp',omega,[1,-1i*omega]);
%   Calls a chunkie routine which creates the system matrix based on
%   provided geometry and kernels. This is the most time-consuming part at
%   the moment
tic
sysmat=chunkermat(overall_geom,kern)+0.5*eye(overall_geom.npt);
sysmat_time=toc;
disp(['System matrix initialization took ',num2str(sysmat_time)]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This is the start of S-matrix evaluation
field_coefs=zeros(numchannels,numchannels);
field_der_coefs=zeros(numchannels,numchannels);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%          Excite each channel (leads+modes in each lead) and solve the scat problem
%%%         This is the computationally intense part
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i2=1:numchannels
    iter_start=tic;
%    i2=2;
    field_coefs_now=[];
    field_der_coefs_now=[];
    inc_channel_info=all_channels_info{i2};
%%          For the given channel (given lead and given mode), set the boundary condition
    rhs=zeros(1,size(overall_geom.r,2)*size(overall_geom.r,3));
    rhs(inc_channel_info.leadend_inds)=inc_channel_info.mode_fun(inc_channel_info.leadend_xy);
    rhs=rhs.';
%%          For the given channel excitation, solve the scat problem
    tic
    all_channels_info{i2}.sigma=sysmat\rhs;
    solve_time=toc;
%%          Find field and field derivative coefficients into every other channel from the given channel excitation
    for i3=1:numleads
    %i3=1;
    tic
    out_lead=leads{i3};
    psiscat_chunkie=chunkerkerneval(overall_geom,kern,all_channels_info{i2}.sigma,out_lead.XY).';
    dpsiscat_dx_chunkie=chunkerkerneval(overall_geom,kern_der,all_channels_info{i2}.sigma,out_lead.XY_ptinfo).';
    [new_psicoefs,~]=lsqr(out_lead.out_mode_matrix,psiscat_chunkie.');
    [new_psidercoefs,~]=lsqr(out_lead.out_mode_matrix,dpsiscat_dx_chunkie.');
    %disp(['i2=',num2str(i2),'; i3=',num2str(i3),'; new_coefs=',num2str(new_psicoefs),'; new_psip_coefs=',num2str(new_psidercoefs)])
    field_coefs_now=[field_coefs_now;new_psicoefs.'];
    field_der_coefs_now=[field_der_coefs_now;new_psidercoefs.'];
    out_lead_eval_time=toc;
%%
%%%%%%%%    Save a plot of field decomposition in a given lead, if desired
        if out_lead.save_lsqr_plots
        fig=figure('units','normalized','position',[0,0,1,1],'visible','off');
        subplot(1,2,1)
        hold on
        plot(out_lead.XY(2,:),real(psiscat_chunkie),'color','blue','linewidth',4)
        plot(out_lead.XY(2,:),real(new_psicoefs*out_lead.out_mode_matrix),'blue--','linewidth',4)
        plot(out_lead.XY(2,:),imag(psiscat_chunkie),'color','red','linewidth',4)
        plot(out_lead.XY(2,:),imag(new_psicoefs*out_lead.out_mode_matrix),'red--','linewidth',4)
        hold off
        legend('Re[Chunkie field]','Re[Scaled mode field]','Im[Chunkie field]','Im[Scaled mode field]','location','eastoutside')
        subplot(1,2,2)
        hold on
        plot(out_lead.XY(2,:),real(dpsiscat_dx_chunkie),'color','blue','linewidth',4)
        plot(out_lead.XY(2,:),real(new_psidercoefs*out_lead.out_mode_matrix),'blue--','linewidth',4)
        plot(out_lead.XY(2,:),imag(dpsiscat_dx_chunkie),'red','linewidth',4)
        plot(out_lead.XY(2,:),imag(new_psidercoefs*out_lead.out_mode_matrix),'red--','linewidth',4)
        hold off
        legend('Re[Chunkie derfield]','Re[Scaled mode derfield]','Im[Chunkie derfield]','Im[Scaled mode derfield]','location','eastoutside')
        exportgraphics(fig,['inc_channel_',num2str(i2),'_lead_',num2str(i3),'.png'])
        end
    end
    field_coefs(:,i2)=field_coefs_now;
    field_der_coefs(:,i2)=field_der_coefs_now;
    iter_time=toc(iter_start);
%    disp(['Channel ',num2str(i2),' total_time=',num2str(iter_time),...
%        '; solve_time=',num2str(solve_time),'; last_lead_eval_time=',num2str(out_lead_eval_time)])
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%          Find S-matrix from all the field and field derivative coefficients
all_channels_beta=zeros(1,numchannels);
for i2=1:numchannels
all_channels_beta(i2)=all_channels_info{i2}.beta;
end

a_in=(field_coefs+1./(1i*all_channels_beta).*field_der_coefs)/2;
a_in_inv=inv(a_in);
a_out=(field_coefs-1./(1i*all_channels_beta).*field_der_coefs)/2;
S=a_out*(a_in_inv);

end