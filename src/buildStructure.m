function [] = buildStructure(fileName, dampingCoefficients, frequency, approximationType, approximationParam)
% Usage: buildstructure('file.txt',[a,b], ApproxType.DerivativeRule,0.5);
%
%
%
    if nargin < 3 || nargin > 5
        error('buildStructure:TooManyInputs', ...
            'requires at least 3 inputs');
    end
    
    if (size(dampingCoefficients) ~= 2)
        error('buildStructure:WrongDampingCoefficients', ...
            'You need to provide atleast 2 parameters for the damping coefficients');
    end
    
    dampingCoefficients = abs(dampingCoefficients);
    frequency = abs(frequency);
    
    if ~exist('approximationType','var')
        approximationType = ApproxType.HalfPower;
    else
        if (approximationType == ApproxType.HalfPower)
            approximationParam=0;
        elseif (approximationType == ApproxType.FreqRange)
            if ~exist('approximationParam','var') || ~isreal(approximationParam)
                disp('Not supplied an approximation param for FreqRange or it is not a number. A default one is provided (10)');
                approximationParam = 10;
            end
           % approximationParam = abs(approximationParam);
        elseif (approximationType == ApproxType.DerivativeRule)
           if ~exist('approximationParam','var') || ~isreal(approximationParam)
                disp('Not supplied an approximation param for DerivativeRule or it is not a number. A default one is provided (10)');
                approximationParam = 0.5;
            end
            approximationParam = abs(approximationParam);
        else
            error('buildStructure:WrongApproximationType', ...
                'Wrong approximation type');
        end
    end

    
    disp(['Damping Coefficients: ' num2str(dampingCoefficients(1)) ' - ' num2str(dampingCoefficients(2))]);
    disp(['Approximation Type: ' char(approximationType)]);
    disp(['Approximation Parameter: ' num2str(approximationParam)]);

    beams = java.util.LinkedList;
    nodes = java.util.LinkedList;
   
    readStructure(fileName,beams.listIterator,nodes.listIterator);
    
    beams = unique(listToMatrix(beams),'rows');
    nodes = unique(listToMatrix(nodes),'rows');
    
    disp(['Read ' num2str(size(beams,1)) ' beams']);
    disp(['Read ' num2str(size(nodes,1)) ' nodes']);
    
    
    [nodesTree,beamsTree] = buildNodesStructure(beams,nodes,frequency,approximationType, approximationParam);
    writeStructure(fileName, nodesTree,beamsTree, dampingCoefficients);
end


function [nodesTree,beamsTree] = buildNodesStructure(beams,nodes,frequency,approximationType, approximationParam)
    
    [nodesTree,beamsTree] = (placeNodes(beams,nodes,frequency, approximationType,approximationParam));  
    nodesTree=round(nodesTree*100)/100.0;
    nodesTreeSize = size(nodesTree,1);
    
    for i=1:nodesTreeSize
       for j=i+1:nodesTreeSize
          if (nodesTree(i,2:3)==nodesTree(j,2:3))
              temp= nodesTree(j,1);
              nodesTree(j,1) = nodesTree(i,1);
              
              for w=1:size(beamsTree,1)
                for z=1:2
                    if (beamsTree(w,z+1)==temp)
                        beamsTree(w,z+1)=nodesTree(i,1);
                    end
                end
              end
              
          end
       end
    end
    
    nodesTree = unique(nodesTree,'rows');
    
    nodesTreeSize = size(nodesTree,1);
    for i=1:nodesTreeSize-1
       if (nodesTree(i+1,1) > nodesTree(i,1)+1)
          temp = nodesTree(i+1,1);
          n = nodesTree(i,1)+1;
          nodesTree(i+1,1) = n;
          for j=1:size(beamsTree,1)
              for z=1:2
                  if (beamsTree(j,z+1)== temp)
                      beamsTree(j,z+1) = n;
                  end
              end
          end
          
       end
    end

end

function [nodesTree,beamsTree] = placeNodes(beams,nodes,frequency,approximationType,approximationParam)
    j=1;
    k=1;
    w=1;
    q=0;
    beamsTree = zeros(1,6);
    for z=1:size(beams,1)
        beam = beams(z,:);
        nlength = beamLength(beam(5:end), frequency, approximationType,approximationParam);
        nnodes = ceil(beam(4)/nlength);
        nlength = beam(4)/nnodes;
        nnodes=nnodes+1;
        bn = zeros(nnodes, 6);
        beam(3) = beam(3)*pi/180;
        q=0;
        for i=1:nnodes
          %  beam(1:3)
            bn(i,1) = w;
            bn(i,2) = beam(1)+nlength*(i-1)*cos(beam(3));
            bn(i,3) = beam(2)+nlength*(i-1)*sin(beam(3));
            bn(i,4:end) = [0,0,0];
            %bn
            
            
            beamsTree=placeNodesBeams(beamsTree,beam,j,i+q,nnodes,w,k);
            %beamsTree(:,1:3)
            
            for p=1:size(nodes,1)
                if(nodes(p,1:2)==bn(i,2:3))
                    bn(i,4:end) = nodes(p,3:end);
                elseif ( i < nnodes)
                    angle = atan2((nodes(p,2)-beam(2)),(nodes(p,1)-beam(1)));
                    r1 = beam(3)-angle;
                    if (r1 < 1e-10)
                        x = nodes(p,1);
                        y = nodes(p,2);
                        x0 = bn(i,2)+nlength*cos(beam(3));
                        y0 = bn(i,3)+nlength*sin(beam(3));
                        if ( ((x < x0 && y <= y0 ) || (x <= x0 && y < y0 ) )&& ( (bn(i,2) < x && bn(i,3) <=y) ||  (bn(i,2) <= x && bn(i,3) <y)))
                            w=w+1;
                            q=q+1;
                            nnodes=nnodes+1;
                            bn = [bn; w, nodes(p,:)];
                            beamsTree=placeNodesBeams(beamsTree,beam,j,i+q,nnodes,w,k);
                        end
                    end
                end
            end
            
            w=w+1;
        end
        if j==1
            nodesTree = [ bn];
        else
            nodesTree = [nodesTree;bn];
        end
        j = j+nnodes-1;
        k=k+1;
    end
end

function [beamsTree] = placeNodesBeams(beamsTree,beam,j,i,nnodes,w,k)

    if i==1
        kn=nnodes-1;
        beamsTree= [beamsTree; k*ones(kn,1) zeros(kn,2) beam(5)*beam(6)*ones(kn,1) beam(8)*beam(6)*ones(kn,1) beam(8)*beam(7)*ones(kn,1)];
        if (j==1)
            beamsTree = beamsTree(2:end,:);
        end
        beamsTree(j+i-1, 2) =w;
        %5 = p, 6=A, 7=J, 8=E
    elseif i==nnodes
        beamsTree(j+i-2, 3) = w; 
    else
        if (size(beamsTree,1) == j+i-2)
            kn=1;
            beamsTree= [beamsTree; k*ones(kn,1) zeros(kn,2) beam(5)*beam(6)*ones(kn,1) beam(8)*beam(6)*ones(kn,1) beam(8)*beam(7)*ones(kn,1)];
        end
            beamsTree(j+i-1, 2) =  w;
            beamsTree(j+i-2, 3) =  w;
    end
end



function [A] = listToMatrix(list)
    if(list.size()>0)
        A=zeros(list.size(), size(list.get(0),1));
        for i=1:list.size()
           A(i,:) = str2num(list.get(i-1))'; 
        end
    end
end
%--------------

function [] = readStructure(fileName, beams, nodes)
    fID = fopen(fileName);
    
    if fID ~= 1
         tline = fgetl(fID);
         while ischar(tline)
            strs = strsplit(tline,{'(',',',')'},'CollapseDelimiters',true);
            strs = strrep(strs, ' ', '');
            strs = strs(~cellfun('isempty',strs));
            parseLine(strs,beams,nodes);
            tline = fgetl(fID);
         end
    else
       error(['Cannot open the file: ' fileName]); 
    end
    
    fclose(fID);
end

function [] = parseLine(line, beams,nodes)
    if size(line,2) > 0
        c = line(1,1);
        switch(c{:})
            case 'b'
                beams.add(line(2:end));
            case 'n'
                nodes.add(line(2:end));
            otherwise
                error(['Unidentified command found during the parsing of the structure file']);
        end
    end
end

function [] = writeStructure(fileName,nodesTree, beamsTree, dampingCoefficients)
    fileName = strcat(fileName, '.inp');
    fID = fopen(fileName,'w');
    if fID ~= 1
        fprintf(fID, '*NODES\n');
        for i=1:size(nodesTree,1)
            fprintf(fID, '%d \t  %d %d %d \t %f \t %f\n', nodesTree(i,1), nodesTree(i,4),nodesTree(i,5),nodesTree(i,6),nodesTree(i,2),nodesTree(i,3));
        end
        fprintf(fID, '*ENDNODES\n');
        fprintf(fID, '*BEAMS\n');
        for i=1:size(beamsTree,1)
            fprintf(fID, '%d \t \t %d %d \t %f \t %f \t %f\n', i, beamsTree(i,2), beamsTree(i,3),beamsTree(i,4),beamsTree(i,5),beamsTree(i,6));
        end
        fprintf(fID, '*ENDBEAMS\n');
        fprintf(fID, '*DAMPING\n');
        fprintf(fID, '%f %f\n', dampingCoefficients(1), dampingCoefficients(2));
    else
       error(['Cannot open the file: ' fileName]);
    end
    fclose(fID);
end

