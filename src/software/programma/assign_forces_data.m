function [force,node,idir,vf,idof,iexit]=assign_forces_data(idb,ngdl,nnod,nf,lingua_it);

force=zeros(ngdl,1);
for i=1:nf
    leggi_dati=true;
    while leggi_dati
        if lingua_it
            lb1='Nodo forzato: ';
            lb2='Direzione della componente della forza 1=x, 2=y, 3=theta: ';
            lb3='Ampiezza della forza ([N]): ';
            dlgTitle=['Componente di forza n. ' int2str(i)];
        else
            lb1='Forced node: ';
            lb2='Forced displacement component 1=x, 2=y, 3=theta: ';
            lb3='Force amplitude ([N]): ';
            dlgTitle=['Force component n. ' int2str(i)];
        end
        def_1=''; 
        def_2=''; 
        def_3=''; 
        prompt={lb1,lb2,lb3};
        default={def_1,def_2,def_3};   
        lineNo=1;
        answer=inputdlg(prompt,dlgTitle,lineNo,default,'on');
        if isempty(answer)
            iexit=true;
            node=[]; 
            idir=[]; 
            vf=[]; 
            idof=[];
            return
        else
            iexit=false;
            node=str2double(answer(1)); 
            idir=str2double(answer(2)); 
            vf=str2double(answer(3)); 
            if node>=1 & node<=nnod & idir>=1 & idir<=3 
                leggi_dati=false;
            end
        end 
    end    
    
    idof=idb(node,idir);
    if idof > ngdl
        if lingua_it
            str_warn='ATTENZIONE: il gdl considerato � vincolato; questa componente di forza non produrr� alcuna deformazione';
        else
            str_warn='WARNING: the considered d.o.f. is subjected to a constraint, this force component will not produce any displacement';
        end
        h=msgbox(str_warn);
        waitfor(h);
        iexit=true;
    else
        force(idof)=force(idof)+vf;
    end
end 




