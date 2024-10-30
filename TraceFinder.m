function [alltraces,Ielectorder] = TraceFinder(dataref, electrode_pts,nelectdisplay,Fiducials)
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
% Authors: Sandhya Tiku, Luize Vasconcelos                                   %
% Nanshu Lu Research Group, UT austin                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Rotate mesh and electrode coordinates to facilate cut traces
    % Chose shift the ref arrays for interconnect
    xshiftarray=0; %mm
    %number of electrodes
    nelect=size(electrode_pts,1);
    meshconnectivity=dataref.ConnectivityList;
    meshcoordinates(:,1)=dataref.Points(:,1);%x
    meshcoordinates(:,2)=dataref.Points(:,2);%y
    meshcoordinates(:,3)=dataref.Points(:,3);%z
    %% determine the order to match with terminal positions (left to right)
    if nelect==19
        % Manually ordered to avoid overlapping order for 10-20 system
        if nelectdisplay==19
            Ielectorder= [12;17;7;11;16;6;19;2;5;15;10;1;18;4;14;9;3;13;8]; % OPTION FOR ALL 19 ELECTRODES
            sorted_pts = electrode_pts(Ielectorder,:);
        else
            if nelectdisplay==15
                Ielectorder= [3;2;17;12;11;16;19;6;5;15;10;4;18;14;9;8;13;1;7];  % OPTION FOR ALL  W/O forehead electrodes (15 electrodes)
                sorted_pts = electrode_pts(Ielectorder,:);
                sorted_pts=sorted_pts(3:17,:);
                nelect=nelectdisplay;
            end
        end
    else
        %sorts electrodes to match terminal array by x
        [sorted_pts,Ielectorder] = sortrows(electrode_pts,1);
    end

    Vertices = dataref.Points; % call vertices coordinates
    %% plots the mesh
    % figure
    % trimesh(dataref,'FaceColor','none','EdgeColor',[0,0,0]+0.8) %plot mesh
    % hold on
    %% Terminals location
    if nelectdisplay==19
        gap=40;%30
    else
        if nelectdisplay==17
            gap=45;%30
        else
            if nelectdisplay==15
                gap=0.7*abs(Fiducials(3,1)-Fiducials(4,1))/2-3;
            else
                gap=0.75*abs(Fiducials(3,1)-Fiducials(4,1))/2-3;
            end
        end
    end
    xt=(-gap:gap*2/(nelect-1):gap)';
    xt=xt+Fiducials(2,1)+xshiftarray; % if the head is shifted in x (x=0 not centered, shift array)
    yt=mean([Fiducials(1,2);Fiducials(2,2)])*ones(nelect,1)+30;
    zt=-40*ones(nelect,1);
    if nelectdisplay==15 %terminal corrections 
        xt(9)=xt(9)-7; %Cz 
        xt(7)=xt(7)+7; %Fz
        %xt(13)=xt(13)-5; %C4 
        xt(12)=xt(12)+5; 
        xt(11)=xt(11)+4; 
        xt(10)=xt(10)+3; 
        %xt(3)=xt(3)+5; %C3 
        xt(4)=xt(4)-5; 
        xt(15)=xt(15)-3;% P8 
        xt(1)=xt(1)+3;% P7
        xt(2)=xt(2)-5; %T7 
        xt(5)=xt(5)-2; %O2
        xt(14)=xt(14)+5; %T8
    end
    terminals = [xt yt zt];
    % hold on
    % plot3(terminals(:,1),terminals(:,2),terminals(:,3),'sg','MarkerSize',10,'LineWidth',2)
    % axis equal

    clear xt
    %% electrode-terminal intermediate ref points
    gap=0.8*abs(Fiducials(3,1)-Fiducials(4,1))/2;
    xt=(-gap:gap*2/(nelect-1):gap)';
    xt=xt+Fiducials(2,1)+xshiftarray; % if the head is shifted in x (x=0 not centered, shift array)
    yt=Fiducials(2,2)*ones(nelect,1)+100; %100
    zt=mean([max(Vertices(:,3));Fiducials(2,3)])*ones(nelect,1);
    
    if nelectdisplay==15 
        xt(15)=xt(15)+30;% P8 
        xt(1)=xt(1)-30;% P7 
        xt(9)=xt(9)+15; %Cz 
        xt(7)=xt(7)-15; %Fz 
        xt(13)=xt(13)+15; %C4
        xt(3)=xt(3)-10; %C3 
        %xt(4)=xt(4)-5; %C3 
        xt(6)=xt(6)+0; %F4
        xt(10)=xt(10)-0; %F3
        xt(2)=xt(2)+5; %T7
        xt(14)=xt(14)-5; %T8
    end

    hold on
    mid_pts = [xt yt zt];

    % runs build_trace function for all the electrodes and adds the trace
    % to alltraces
    alltraces = cell(nelect,1);

    % chosen_trace = 8; % CHOOSE TRACE, IF 0, PLOT ALL
    for i = 1:nelect   
      
        % plot3(mid_pts(i, 1), mid_pts(i, 2), mid_pts(i, 3),'or','MarkerSize',15,'LineWidth',2);
        cut_pts = build_trace(meshconnectivity, meshcoordinates, sorted_pts(i, 1:3), mid_pts(i, 1:3), terminals(i,1:3));
        alltraces{i} = cut_pts;

    end

end

function [a,b,c,d] = equation_plane(p1, p2, p3)
    normal = cross(p1 - p2, p1 - p3);
    a = normal(1,1);
    b = normal(1,2);
    c = normal(1,3);
    d = -1*(normal(1,1)*p2(1,1) + normal(1,2)*p2(1,2) + normal(1,3)*p2(1,3));
end

function bool = elem_intersect(A, B, C, D, p1, p2, p3)
    % All points that have the same sign are on the same side of the plane
    % therefore all elements that have 2 points with one sign and 2 points
    % with the other sign intersect the plane
    val1 = A*p1(1,1) + B*p1(1,2) + C*p1(1,3) + D;
    val2 = A*p2(1,1) + B*p2(1,2) + C*p2(1,3) + D;
    val3 = A*p3(1,1) + B*p3(1,2) + C*p3(1,3) + D;

    % Goes through all the combinations to see when there are 2 points on
    % the same side of the plane in which case it returns true
    if val1 > 0 && val2 > 0 && val3 > 0
        bool = [false val1 val2 val3];
    elseif val1 < 0 && val2 < 0 && val3 < 0
        bool = [false val1 val2 val3];
    else

        if (val1 > 0 && val2 < 0 && val3 < 0)
            bool = [true p1 p2 p3];
        elseif (val1 < 0 && val2 > 0 && val3 < 0)
            bool = [true p2 p1 p3];
        elseif (val1 < 0 && val2 < 0 && val3 > 0)
            bool = [true p3 p1 p2];
        elseif (val1 > 0 && val2 > 0 && val3 < 0)
            bool = [true p3 p1 p2];
        elseif (val1 > 0 && val2 < 0 && val3 > 0)
            bool = [true p2 p1 p3];
        elseif (val1 < 0 && val2 > 0 && val3 > 0)
            bool = [true p1 p2 p3];
        end
    end

end

function cut_pts = find_intersection_pts(connectivity_list, points, A,B,C,D, p1, p3) 
    countint=0;%count intercepts
    intersection_pts=zeros(1,6); %reserve memory

    for i = 1:length(connectivity_list) - 1
        % Define variables for the 3 points in the connectivitylist element
        cl1 = connectivity_list(i,1);
        cl2 = connectivity_list(i,2);
        cl3 = connectivity_list(i,3);
        % Create an array of the x,y,z values of the 3 points in the
        % element
        
        %Check if the element intersects the plane
        % The first index in the array represents whether it intersects
        % with the plane or not
        bool = elem_intersect(A, B, C, D, points(cl1,:), points(cl2,:), points(cl3,:));
        if (bool(1,1) == true)

            opp = bool(2:4);
            same1 = bool(5:7);
            same2 = bool(8:10);

            u1 = (A*opp(1,1) + B*opp(1,2) + C*opp(1,3) + D)/(A*(opp(1,1) - same1(1,1)) + (B*(opp(1,2) - same1(1,2))) + C*(opp(1,3) - same1(1,3)));%not used - WHY? -luize
            i1 = opp + u1 * (same1 - opp);
            u2 = (A*opp(1,1) + B*opp(1,2) + C*opp(1,3) + D)/(A*(opp(1,1) - same2(1,1)) + (B*(opp(1,2) - same2(1,2))) + C*(opp(1,3) - same2(1,3)));
            i2 = opp + u2 * (same2 - opp);

            if (ismember(round(i2,0), round(intersection_pts(:,1:3),0),'rows') == false)
                countint = 1 + countint;
                intersection_pts(countint,1:3) = i2;
                intersection_pts(countint,4:6) = cross(points(cl1,:) - points(cl2,:), points(cl1,:) - points(cl3,:));
            end
            if (ismember(round(i1,0), round(intersection_pts(:,1:3),0),'rows') == false)
                countint = 1 + countint;
                intersection_pts(countint,1:3) = i1;
                intersection_pts(countint,4:6) = cross(points(cl1,:) - points(cl2,:), points(cl1,:) - points(cl3,:));
            end        

        end
    end
    intersection_pts = sortrows(intersection_pts,1);
    cut_pts = cut_sort(p1, p3, intersection_pts, countint);

end

function cut_pts = build_trace(connectivity_list, points, p1, p2, p3)

    [A,B,C,D] = equation_plane(p1, p2, p3) ; %Build plan containing electrode, midpoint, and terminal point
    cut_pts = find_intersection_pts(connectivity_list, points, A,B,C,D, p1, p3);

end

function cut_pts = cut_sort(p1, p3, intersection_pts, list_length)
%p1 : electrode
%p3 : terminal
%p1_adj : electrode in front

%Keep only points in between terminal and electrode by creating a plane
%parallel to the x-axis containing both points and deleting everything
%below and all the points below the Z coordinate of the terminal.
    cut_pts = zeros(length(list_length),6);
    j = 0;
    for i = 1:list_length
        % point that will create a plane that is parallel to x so that all
        % the points above the plane define the interconnect
        pterminalref=[p3(1)+10 p3(2) p3(3)]; 
         
        [A, B, C, D]=equation_plane(p1, p3, pterminalref); %find the plane terminal and electrode

        pt = intersection_pts(i,1:3); %take coordinates of the current point
        val1 = A*pt(1) + B*pt(2) + C*pt(3) + D;

        if val1 < 0 && pt(3) > p3(3) % keep if current point in between terminal and electrode
            j = j + 1;
            cut_pts(j, :) = intersection_pts(i,:);
        end
    end
    ordered_trace = zeros(length(cut_pts),6);
    dist  = zeros(length(cut_pts),1);
    %Sort trace points from terminal to electrode
    cPt = p1; % current point
    candidates = cut_pts; %initialize candidate list
    for j = 1: size(cut_pts,1)
        % find the closest point to current point
        for i=1:length(cut_pts)
            dist(i,1) = sqrt(((cPt(1,1) - candidates(i,1))^ 2) + ((cPt(1,2) - candidates(i,2))^ 2) + ((cPt(1,3) - candidates(i,3))^ 2));
        end
        [~, Imindist]=min(dist);
        ordered_trace(j,1:6)=cut_pts(Imindist,:); % store closest point in ordered_trace
        candidates(Imindist,1:6)=NaN; %erase point from candidate list
        cPt=ordered_trace(j,1:3); %make the closest point the current point
    end

    cut_pts=ordered_trace;

    %% inspect if trace pts are in correct order
%     plot3(cut_pts(:,1), cut_pts(:,2), cut_pts(:,3),'-*m','LineWidth',2)
%     hold on

end
