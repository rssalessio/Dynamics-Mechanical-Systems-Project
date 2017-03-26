function [str_title,nfig,nfig_aperte]=static(caso,handles,nfig,nfig_aperte)

nnod=handles.nnod;
ngdl=handles.ngdl;
M=handles.M;
K=handles.K;
size=handles.size;
idb=handles.idb;
incid_masse=handles.incid_masse;
nm_conc=handles.nm_conc;
xy=handles.xy;
lingua_it=handles.lingua_it;
beams=handles.beams;
funi=handles.funi;
Tbeams=handles.Tbeams;

% Computes the static deflection of the assigned structure, under its own
% weight or under a specified concentrated force

% keyboard

%ico = menu('Type of static load','Structure weight','Concentrated nodal loads') ;
if caso == 'g'
    % The vector of nodal loads corresponding to structure weight is computed
    gravity = zeros(ngdl,1) ;
    for j=1:nnod
        kk=idb(j,2) ;
        if kk<=ngdl
            gravity(kk)=-9.81;
        end
    end
    force = M*gravity ;
    totweight = sum(force) ;
%    print_stringa(['Peso complessivo della struttura : ' num2str(totweight)],handles);
%    print_stringa(['Total weight of the structure : ' num2str(totweight)],handles);
    if lingua_it
        str_title='Deformata statica (gravit�)';
    else
        str_title='Static Deflection (Structure weight)';
    end
else
    % The vector of nodal forces is received as input from the user    
    [nf,iexit]=assign_nf(lingua_it);
    if iexit;return;end;

    [force,node,idir,vf,idof,iexit]=assign_forces_data(idb,ngdl,nnod,nf,lingua_it);
    if iexit;return;end;
    if lingua_it
        str_title='Deformata statica (forze concentrate)';    
    else
        str_title='Static Deflection (Concentrated loads) ';
    end
end     %end of if construct own weight / concentrated loads

% Calculation of static displacements under the assigned static loads
defo=K\force;

% display of results
% 1. The static deformation vector is displayed (text) on the Matlab Command window

if lingua_it
    print_stringa('Deformata statica ',handles);
else
    print_stringa('Static Deflection ',handles);
end
vise(idb,ngdl,defo,handles) ;

% The dynamic deformation vector is plotted with an appropriate scale factor
nfig=nfig+1;
nfig_aperte=nfig_aperte+1;
fscala=0.2*size/max(abs(defo));
set_fscala(fscala,defo,str_title,'sta',nm_conc,incid_masse,idb,xy,nfig,lingua_it,beams,Tbeams,funi);


