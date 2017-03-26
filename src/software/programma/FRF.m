function [fre,x,xpp,vett_str,nfig_end]=FRF(handles,nfig,type);

nnod=handles.nnod;
ngdl=handles.ngdl;
M=handles.M;
K=handles.K;
R=handles.R;
size=handles.size;
idb=handles.idb;
incid_masse=handles.incid_masse;
nm_conc=handles.nm_conc;
xy=handles.xy;
lingua_it=handles.lingua_it;
beams=handles.beams;
funi=handles.funi;
Tbeams=handles.Tbeams;

fre=[];
x=[];
xpp=[];
vett_str=[];
nfig_end=nfig;


% Computes the Frequency Response Function for the assigned structure,
% given the structural matrices M, K, R, and displays the results in terms
% of bode diagrams and structure dynamic deflection


% Step 1 input of FRF analysis parameters
[nf,iexit]=assign_nf(lingua_it);
if iexit;return;end;

[force,node,idir,vf,idof,iexit]=assign_forces_data(idb,ngdl,nnod,nf,lingua_it);
if iexit;return;end;

if lingua_it
    lb1='Estremo inferiore del campo di frequenza da analizzare [Hz]: ';
    lb2='Estremo superiore del campo di frequenza da analizzare [Hz]: ';
    lb3='Ampiezza del passo in frequenza [Hz]: ';
    dlgTitle='Impostazione dei dati in frequenza';
else
    lb1='Lower boundary of frequency range for the analysis [Hz]: ';
    lb2='Upper boundary of frequency range for the analysis [Hz]: ';
    lb3='Amplitude of frequency step to span the frequency range [Hz]: ';
    dlgTitle='Frequency range data';
end

def_1=''; 
def_2=''; 
def_3=''; 
prompt={lb1,lb2,lb3};
default={def_1,def_2,def_3};   
lineNo=1;
answer=inputdlg(prompt,dlgTitle,lineNo,default,'on');
if isempty(answer)
    return;
else
    fmin=str2double(answer(1)); 
    fmax=str2double(answer(2)); 
    df=str2double(answer(3)); 
end 

% keyboard

np=round((fmax-fmin)/df)+1;

% Step 2: FRF analysis - results are stored in a complex matrix x with ngdl rows and np columns,
% each row representing the FRF for a single degree of freedom of the
% structure and each column representing a specified value of frequency

if lingua_it
    print_stringa('Attendere prego ...',handles);
else
    print_stringa('Please wait, this calculation may take a while ...',handles);
end

delta_stampa=round(np/10);
i_stampa=delta_stampa;
for i=1:np
    if ngdl>500 & i==i_stampa; print_stringa(['... ' num2str(i) '/' num2str(np)],handles); i_stampa=i_stampa+delta_stampa; end
    fre(i)=fmin+(i-1)*df;
    ome=2*pi*fre(i);
    a=-ome^2*M+j*ome*R+K;
    x(:,i)=a\force;
    xpp(:,i)=-a\force*ome^2;
end


if type=='s'
    mod=abs(x);
    fas=angle(x);
else
    mod=abs(xpp);
    fas=angle(xpp);
end

if lingua_it
    if type=='s'
        str_type=' - comp. spost.: ';
        str_ylabel_mod='Amp. [m/N]';
    else
        str_type=' - comp. acc.: ';
        str_ylabel_mod='Amp. [ms^2/N]';
    end
    str_ylabel_fase='Fase [deg]';
    str_nodo='Nodo: ';
    str_gdl=' - g.d.l.: ';
    str_menu1='FRF output';
    str_menu2='Diagrammi di Bode - scale lineari';
    str_menu3='Diagrammi di Bode - scale logaritmiche';
    str_menu4='Deformata dinamica per una specificata frequenza';
    str_menu5='Esci';
    str_msg_gdl='Il gdl considerato � vincolato, i grafici non sono disponibili ';
else
    if type=='s'
        str_type=' - displ. comp.: ';
        str_ylabel_mod='Amp. [m/N]';
    else
        str_type=' - acc. comp.: ';
        str_ylabel_mod='Amp. [ms^2/N]';
    end
    str_ylabel_fase='Phase [deg]';
    str_nodo='Node: ';
    str_gdl=' - d.o.f.: ';
    str_menu1='FRF output';
    str_menu2='Bode diagrams - linear scales';
    str_menu3='Bode diagrams - log scales';
    str_menu4='Dynamic deflection for a specified frequency';
    str_menu5='Exit';
    str_msg_gdl='The considered d.o.f. is subjected to a constraint, no output can be produced ';
end

% outputs
iout=0;
ra2deg=180/pi;

while iout ~= 4
    if type=='s'
        iout=menu(str_menu1,str_menu2,str_menu3,str_menu4,str_menu5);
    else
        iout=menu(str_menu1,str_menu2,str_menu3,str_menu5);
        if iout==3;iout=4;end;
    end

    if or((iout == 1),(iout == 2)) == 1
        % Bode diagrams (linear or log)
        [node,idir,iexit]=assign_FRFout_data(nnod,lingua_it);
        if iexit;return;end;
        
        idof=idb(node,idir);
        if idof > ngdl
            msgbox(str_msg_gdl)
        else
            nfig_end=nfig_end+1;
            figure(nfig_end);
            subplot(211)
            if iout == 1
                plot(fre,mod(idof,:));grid
            else
                % first point is omitted to avoid the possibility of having freq(1)=0, not suitable for logscale representation        
                loglog(fre(2:np),mod(idof,2:np));grid
            end
            xlabel('Freq. [Hz]')
            ylabel('Amp.')
            str_title=[str_nodo,int2str(node),str_type,int2str(idir),str_gdl,int2str(idof)];
            vett_str=strvcat(vett_str,str_title);
            title(str_title)
            subplot(212)
            if iout == 1
                plot(fre,fas(idof,:)*ra2deg);grid
                a=axis;a(3)=-200;a(4)=200;axis(a)
            else
                % first point is omitted to avoid the possibility of having freq(1)=0, not suitable for logscale representation        
                semilogx(fre(2:np),fas(idof,2:np)*ra2deg);grid
                a=axis;a(3)=-200;a(4)=200;axis(a)
            end
            xlabel('Freq. [Hz]')
            ylabel(str_ylabel_fase)
        end     %end of check that the displacement to be plotted is not constrained
        
    elseif iout == 3

        % display of dynamic deflection
        [fout,iexit,iout]=set_frequency(fmin,fmax,np,df,fre,lingua_it);
        
        % the complex vector of motion x is projected along the phase of the d.o.f. with maximum amplitude
        [y,imax]=max(abs(x(:,iout)));
        defo=real(x(:,iout)/x(imax,iout))*abs(x(imax,iout));

        if lingua_it
            str_print=['Deformata dinamica a f =',num2str(fout),' [Hz]'];
            str_title=['Deformata dinamica a f =',num2str(fout),' [Hz] - fattore di scala '];
        else
            str_print=['Dynamic Deflection at f =',num2str(fout),' [Hz]'];
            str_title=['Dynamic Deflection at f =',num2str(fout),' [Hz] - Scale factor '];
        end

        % The dynamic deformation vector is displayed (text) on the Matlab Command window
        print_stringa(str_print,handles);
        vise(idb,ngdl,defo,handles) ;
        
        % The dynamic deformation vector is plotted with an paapropriate scale factor
        fscala=0.2*size/y;
        nfig_end=nfig_end+1;
        vett_str=strvcat(vett_str,str_title);
        set_fscala(fscala,defo,str_title,'frf',nm_conc,incid_masse,idb,xy,nfig_end,lingua_it,beams,Tbeams,funi);
    end   %end of if structure
end     % end of while structure




