function [overall_geom] = create_stadiumgeom_v1(cava,cavb,stadr,numleads,leadsl,leadsd,leadscens,numell,ellcens,ellsemi,elltheta)
%{
%%%     Get all params
cava=cav_params.cava;
cavb=cav_params.cavb;
stadr=cav_params.stadr;

numleads=leads_params.numleads;
leadsl=leads_params.leadsl;
leadsd=leads_params.leadsd;
leadscens=leads_params.leadscens;

numell=ell_params.numell;
ellcens=ell_params.ellcens;
ellsemi=ell_params.ellsemi;
elltheta=ell_params.elltheta;
%}
%%%%%%%%%%%%%%%%%%%     Create the chunkgraph
%nverts=4;
nverts=4+4*numleads;
%%%     Edge incidence matrix: in our case, just
%%%    M=[1,2,3,...,nverts;
%%%       2,3,4,...,1]
%%%     To be specific, M(:,i) gives the numbers of starting and ending
%%%     vertices for edge i
edgeinc=[1:nverts;
         circshift(1:nverts,-1)];
%%%%%%%%%     Define all the vertices of the outer cavity
verts=[-cava/2,cava/2;
       cavb/2,cavb/2];
for i2=1:numleads
current_lead_verts=[cava/2,cava/2+leadsl(i2),cava/2+leadsl(i2),cava/2;
                    leadscens(i2)+leadsd(i2)/2,leadscens(i2)+leadsd(i2)/2,leadscens(i2)-leadsd(i2)/2,leadscens(i2)-leadsd(i2)/2];
verts=[verts,current_lead_verts];
end
verts_end=[cava/2,-cava/2;
       -cavb/2,-cavb/2];
verts=[verts,verts_end];
%%%%%%%%%     Define straight edges between all the vertices
edgefuns=cell(1,nverts);
for i2=1:(nverts-1)
edgefuns{i2}=@(t) (1-t(:).').*verts(:,i2)+t(:).'.*verts(:,i2+1);
end
%edgefuns{nverts}=@(t) (1-t(:).').*verts(:,nverts)+t(:).'.*verts(:,1);
%%%%%%%%%     Define a curved stadium-like edge between the last pair of
%%%%%%%%%     vertices
phi0=pi+asin(cavb/2/stadr);
deltax=-cava/2-stadr*cos(phi0);
alpha=(phi0/pi-1)*2*pi;
edgestadium=@(t) [deltax+stadr*cos( -alpha*t(:).'+phi0 );stadr*sin( -alpha*t(:).'+phi0 )];
edgefuns{nverts}=edgestadium;
%%%%%%%%%     Define all the ellipses
for i2=1:numell
edgeinc=[edgeinc,nan(2,1)];
ellfun=@(t) ellcens(:,i2)+ellsemi(:,i2).*[cos( 2*pi*(t(:).'-elltheta(i2)) );sin( 2*pi*(t(:).'-elltheta(i2)) )];
edgefuns=[edgefuns,{ellfun}];
end

overall_geom=chunkgraph(verts,edgeinc,edgefuns);

end
