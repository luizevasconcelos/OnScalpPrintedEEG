


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% This postprocessor loads a toolpath and transforms it in 5-axis G-code  %
% Toolpath input file format (xlsx):
% Columns 1 to 3 contain x,y,z coordinates of the path                    %
% Columns 4 to 6 contain the coordinates of the normal vector for each pt %                                                         %
%                                                                         %
% Authors: Eric Li, Luize Vasconcelos (luizescalco@gmail.com)             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Workspace Setup
clear
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                 USER INPUTS                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Choose file
subject='hair_mannequin';
electplot=7; %chose electrode to plot paths
chooseprint=0;% # from Names % 0 to print all
choosetrace=1;%1 for interconnect 2 for electrode

%% CHOOSE STANDARD ROTATION
rot0z=90;%first rotation
rot0y=65;%second rotation

%% MATCH PHYSICAL AND DIGITAL HEAD POSITIONS BASED ON LANDMARKS
pm=[0,0 0; -81.56 -39.4 -4.25; -120.21 -17.33 -27.56]; % input physical landmarks


%% Jetting frequency (jet/s)
jetfreq=12; %Hz

%% Version name/ rotate mesh
ExportNameVersion=append(num2str(rot0y),'deg'); % Name to append to export file name: 'filename''spec'.txt
%% tool feed rate
speed=200; %%mm/min for linear movement and deg/min for rotatory
maxspeed=speed*3;
transitionspeed=1000;
%% angle during (0,0,0)
c0=90; % either: C= 0 or 90 deg; B = 0
b0=0;
%% Choose the jetting distance from the surface
gapSize = 8;
calgapZ = 37-5.88;

%% Input the machine information
% Vertical distance from the B-axis of rotation and the tip of the jetting device
d =57.59;
% Horizontal distance from the C-axis of rotation and the tip of the jetting device
a =145.5+7.7+27.5; % longitudinal config


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                 PROCESSING                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LOAD SCANNED DATA
% Load raw toolpaths ------------------------------------------------------
alldata=load(append(subject,'/','InterconnectElectrodePaths.mat'));
alldata=alldata.InterconnectElectrodePaths; % load all date
Interconnect=alldata(:,1); % interconnect data
Electrode=alldata(:,2);% electrode data
Names=alldata(:,3);% electrode names
% Load scanned landmarks --------------------------------------------------
Landmarks=load(append(subject,'/','Landmarks.mat'));
ps=Landmarks.Landmarks;
nelect=size(alldata,1); %number of electrodes
% Load scanned head mesh --------------------------------------------------
HeadMesh=load(append(subject,'/','HeadMesh.mat'));
HeadMesh=HeadMesh.dataref;
% Load landmark names -----------------------------------------------------
LandmarkNames=load(append(subject,'/','LandmarkNames.mat'));
LandmarkNames=LandmarkNames.LandmarkNames;

%% MAKE LANDMARK POINT 1 THE DIGITAL ORIGIN
HeadMeshPoints=HeadMesh.Points-ps(1,:); % head mesh
HeadMesh = triangulation(HeadMesh.ConnectivityList, HeadMeshPoints);
for i=1:nelect
    Interconnect{i}(:,1:3)=Interconnect{i}(:,1:3)-ps(1,:);% interconnect
    Electrode{i}(:,1:3)=Electrode{i}(:,1:3)-ps(1,:);% electrode
end
ps=ps-ps(1,:); % landmarks

%% ROUGH ROTATION
HeadMeshPoints=(roty(rot0y)*rotz(rot0z)*HeadMesh.Points')'; %headmesh
HeadMesh = triangulation(HeadMesh.ConnectivityList, HeadMeshPoints);
for i=1:nelect
    Interconnect{i}(:,1:3)=(roty(rot0y)*rotz(rot0z)*Interconnect{i}(:,1:3)')'; %interconnect pts
    Interconnect{i}(:,4:6)=(roty(rot0y)*rotz(rot0z)*Interconnect{i}(:,4:6)')'; %interconnect normals
    Electrode{i}(:,1:3)=(roty(rot0y)*rotz(rot0z)*Electrode{i}(:,1:3)')'; %electrode pts
    Electrode{i}(:,4:6)=(roty(rot0y)*rotz(rot0z)*Electrode{i}(:,4:6)')'; %electrode normals
end
ps=(roty(rot0y)*rotz(rot0z)*ps')'; %landmarks
ps=ps(1:3,:);
ps=ps([1 3 2],:);
pm=pm([1 3 2],:);


%% Pose matching
%Plot vector normal to plane formed by the three digital and physical landmarks
psR=scan2phys(ps,ps,pm); %1st rotation
psb=ps([1 3 2],:);
pmb=pm([1 3 2],:);
psR2=scan2phys(psb,psb,pmb); %2nd rotation
psRm=(psR+psR2([1 3 2],:))/2; %average

% Plot pose matching landmarks
figure(1)
hold on
title('Landmark pose matching')
plot3(pm(:,1),pm(:,2),pm(:,3),'k*','MarkerSize', 15,'LineWidth',2)
plot3(ps(:,1),ps(:,2),ps(:,3),'b*','MarkerSize', 15,'LineWidth',2)
plot3(psR2(:,1),psR2(:,2),psR2(:,3),'r*','MarkerSize', 15,'LineWidth',2)
trimesh(HeadMesh,'FaceColor','none','EdgeColor',[0,0,0]+0.8) %plot mesh
axis equal
legend('Physical landmarks','Scanned landmarks before pose matching','Scanned landmarks after pose matching')

HeadMeshPointsR=scan2phys(HeadMesh.Points,ps,pm); %headmesh
HeadMeshR = triangulation(HeadMesh.ConnectivityList, HeadMeshPointsR);

for i=1:nelect
    InterconnectR{i}(:,1:3) = scan2phys(Interconnect{i}(:,1:3),ps,pm);
    InterconnectR{i}(:,4:6) = scan2phys(Interconnect{i}(:,4:6),ps,pm);
    ElectrodeR{i}(:,1:3) = scan2phys(Electrode{i}(:,1:3),ps,pm);
    ElectrodeR{i}(:,4:6) = scan2phys(Electrode{i}(:,4:6),ps,pm);
end

Interconnect=InterconnectR;
Electrode=ElectrodeR;
HeadMesh=HeadMeshR;
ps=psR;

%% Plot rotated traces
figure(3)
plotfig(Interconnect,Electrode,HeadMesh,ps,pm,nelect,Names)

%% Determine if chosen toolpath is electrode or interconnect and if print all or single one
stopnext=0;
for m=1:nelect
    if choosetrace==1 %chose to print interconnect
        if chooseprint==0 %chose to print all interconnect at once
            clear data
            data=Interconnect{m};
        else
            if stopnext==0
                data=Interconnect{chooseprint}; %chose to a single interconnect
                stopnext=1;
            else
                break % stop loop if printing a single interconnect
            end
        end
    else
        if chooseprint==0 %chose to print electrode
            clear data
            data=Electrode{m}; %chose to print all electrodes at once
        else
            if stopnext==0
                data=Electrode{chooseprint}; %chose to a single interconnect
                stopnext=1;
            else
                break % stop loop if printing a single electrode
            end
        end
    end
    npts= size(data,1);
    %% Calculate Vertical Offset
    t = d + gapSize; % adapter length + gap size offset
    %% Store Trace Coordinates
    g = data(:,1:3); % trace position vector
    %% Store Trace Normal Vectors
    en = data(:,4:6); % trace normal vector
    %normalize normal vectors
    for i = 1:npts
        en(i, :) = en(i, :) / norm(en(i, :));
    end

    %% Visualize original paths
    if m==electplot || chooseprint~=0
        figure(4)
        subplot(2,4,1)
        plot3(g(:,1),g(:,2),g(:,3),'-*') %plot path in 3D
        hold on  %plot normals in 3D
        quiver3(g(:,1),g(:,2),g(:,3),en(:,1),en(:,2),en(:,3),'AutoScaleFactor',0.5) %plot normals
        axis equal
        xlabel('x (mm)')
        ylabel('y (mm)')
        zlabel('z (mm)')
        title('Original trace path')
    end
    %% Standalone plot
    figure(3)
    plot3(g(:,1),g(:,2),g(:,3),'-*') %plot path
    hold on
    quiver3(g(:,1),g(:,2),g(:,3),en(:,1),en(:,2),en(:,3),1) %plot normals
    axis equal
    xlabel('x (mm)')
    ylabel('y (mm)')
    zlabel('z (mm)')
    title('Electrode-interconnect traces & normal vectors')
    %% Visualize rotated path
    if m==electplot || chooseprint~=0
        figure(4)
        subplot(2,4,2)
        hold on
        plot3(g(:,1),g(:,2),g(:,3),'-*') %plot path
        quiver3(g(:,1),g(:,2),g(:,3),en(:,1),en(:,2),en(:,3),1) %plot normals
        axis equal
        xlabel('x (mm)')
        ylabel('y (mm)')
        zlabel('z (mm)')
        title('Rotated  trace path')
    end
    %% Choose tool angle position for setting machine x,y,z zeros---------------------------
    if c0==0 && b0==0
        g(:,1)=g(:,1)-a;
        g(:,3)=g(:,3)+d+calgapZ;
    else
        if c0==90 && b0==0
            g(:,2)=g(:,2)-a;
            g(:,3)=g(:,3)-(d+calgapZ);
        else
            if c0==0 && b0==90
                g(:,2)=g(:,2)-a-d-calgapZ;
            else
                msg = 'Choose the B,C angles that the tool will be in when setting x,y,z machine zeros';
                error(msg)
            end
        end
    end
    %% Visualize trace path on machine frame
    if m==electplot || chooseprint~=0
        figure(4)
        subplot(2,4,2)
        plot3(g(:,1),g(:,2),g(:,3),'-*') %plot path
        hold on
        quiver3(g(:,1),g(:,2),g(:,3),en(:,1),en(:,2),en(:,3),'AutoScaleFactor',0.5) %plot normals
        title('Trace path on trace frame vs machine frame')
    end
    %% Compute C angle
    cAxisAngles = zeros(npts, 1);
    for i = 1:npts
        cAxisAngles(i) = findCaxisAngle(en, i); % C angle
    end
    bAxisAngles = zeros(npts, 1);
    for i = 1:npts
        bAxisAngles(i) = findBaxisAngle(en, i); % B angle
    end


    %% Probihit sudden C-angle sign flip (larger than 40deg in one increment)
    indexflip = [];
    nflips=0;
    countflips=zeros(npts,1);
    for i=1:npts-1
        k=npts-i;
        if ((cAxisAngles(k+1)>0 && cAxisAngles(k)<0) || (cAxisAngles(k+1)<0 && cAxisAngles(k)>0)) && abs(cAxisAngles(k))>20 && abs(cAxisAngles(k+1))>20
            countflips(k)=1;
            indexflip=[indexflip k]; %track number of corrections per trace
            nflips=1+nflips;
        end
    end

    if nflips==1
        cAxisAngles(1:indexflip(1))=-cAxisAngles(1:indexflip(1));
        bAxisAngles(1:indexflip(1))=-bAxisAngles(1:indexflip(1));
    else
        if nflips>1
            for f=1:nflips
                cAxisAngles(1:indexflip(f))=-cAxisAngles(1:indexflip(f));
                bAxisAngles(1:indexflip(f))=-bAxisAngles(1:indexflip(f));
            end
        end
    end

    %% Correcting for Tool Length and B & C Axis Offsets
    % NOTE: x & y correction depends on the C angle pick rule (C=180 vs C=0)
    % x coordinates
    for i = 1:npts
        xOffset = en(i, 1) * t - abs(cos(deg2rad(cAxisAngles(i)))) * a;
        g(i, 1) = g(i, 1) + xOffset;
    end

    % y coordinates
    for i = 1:npts
        yOffset = en(i, 2) * t + sin(deg2rad(cAxisAngles(i))) * a;
        g(i, 2) = g(i, 2) + yOffset;
    end

    % z coordinates
    for i = 1:npts

        crossP = cross(en(i, :), [0 0 1]);

        if abs(crossP) < 0.01 %avoid NaN for normal vector = z coordinate
            totalOffset = en(i, :) * t;
        else
            totalOffset = en(i, :) * t - crossP ./ norm(crossP) * a;
        end
        zOffset = totalOffset(3);
        g(i, 3) = g(i, 3) + zOffset;

    end


    %% round to machine resolution g-code
    g = round(g, 2);
    cAxisAngles = round(cAxisAngles, 2);
    bAxisAngles = round(bAxisAngles, 2);


    %% Set the speed of the nozzle to be constant based on variable translational speed of the C axis

    % Feed Per Minute (G94)
    % G94 G-code is a modal G-code. G94 instructs the control to interpret feed commands as
    %  *inches/minute or mm/minute for linear moves.
    %  *degrees/minute for rotary moves.
    %  *inches/minute or mm/minute for a combination of linear and rotary moves.
    % When a combination of linear and rotary moves is programmed, the rotary moves match the time it takes to make the linear moves.
    %
    % The G94 function selects feed F in mm/min or inches/minute. When this function is active the feed values will be programmed as follows: F50, F150, F500, F2000 and so forth.
    % G94 (feed per minute) G-code is used to perform movements with work feed when the spindle is stationary, or when it is necessary to release the axis feed from the spindle revolutions (e.g.: when milling with motor driven tools or live tooling).

    dispCaxis=zeros(npts,1);
    dispNozzle=zeros(npts,1);
    speedCaxis=zeros(npts,1);
    for i=2:npts
        dispCaxis(i,1)= sqrt((g(i,1)-g(i-1,1))^2+(g(i,2)-g(i-1,2))^2+(g(i,3)-g(i-1,3))^2); % translation of the C axis
        dispNozzle(i,1)= sqrt((data(i,1)-data(i-1,1))^2+(data(i,2)-data(i-1,2))^2+(data(i,3)-data(i-1,3))^2); % displacement of the nozzle
        speedCaxis(i,1)=speed*dispCaxis(i)/dispNozzle(i); % translation speed of the C axis to maintain nozzle speed constant
    end
    speedCaxisbounded=speedCaxis;
    for i=2:npts % Prevent C-axis translation to exceed max safe speed
        if speedCaxis(i,1)>maxspeed
            speedCaxisbounded(i,1)=maxspeed;
        end
    end
    speedCaxisbounded(1)= transitionspeed;
    timeperstep=(dispCaxis./speedCaxisbounded); %min
    timeperelectrode(m,1)=sum(timeperstep); %min
    timeperelectrode(m,2)=sum(dispNozzle./speed); %min
    njets(m,1)=timeperelectrode(m,1)*60*jetfreq; %totaltime (s)* jerfreq(#/s)
    if m==electplot || chooseprint~=0
        fig4plot(data,g,cAxisAngles,bAxisAngles) % plot print traces for electrode selected for display
    end
    speedCaxisbounded = round(speedCaxisbounded, 0); %round speed
    nozzlepath{i}=data;
    %% Build Matrix
    if chooseprint==0
        gcode0{m,1} = [g bAxisAngles cAxisAngles speedCaxisbounded]; %print all electrodes/interconnects at once
    else
        gcode0{1} = [g bAxisAngles cAxisAngles speedCaxisbounded]; %print a single electrode/interconnect
    end
end

%% alternate reversing paths to save print time
for i=1:size(gcode0)
    if mod(i,2)==0 %i is odd
        gcode0{i}=flip(gcode0{i},1); %flip
        %else %i is even
        % do nothing
    end
end

Zsafe=round(max(HeadMeshPointsR(:,3))+25);
%% Merge all traces and add transition
for i=1:size(gcode0)-1
    %% machine coord

    gcode0{i}=[gcode0{i}; gcode0{i}(end,:)]; %repeat last row
    gcode0{i}(end,7)=11; %add marker to start jetting: 11=ends, 10=starts

    gcode0{i}=[gcode0{i}; gcode0{i}(end,1:6) 0]; %repeat last row
    gcode0{i}(end,6)=transitionspeed;
    gcode0{i}(end,3)=Zsafe; %modify Z of the last row to Zsafe

    gcode0{i}=[gcode0{i}; gcode0{i}(end,1:6) 0]; %repeat last row
    gcode0{i}(end,4)=gcode0{i+1}(1,4); %modify B of the last row to next B
    gcode0{i}(end,5)=gcode0{i+1}(1,5); %modify C of the last row to next C

    gcode0{i}=[gcode0{i}; gcode0{i}(end,1:6) 0]; %repeat last row
    gcode0{i}(end,1)=gcode0{i+1}(1,1); %modify X of the last row to next X
    gcode0{i}(end,2)=gcode0{i+1}(1,2); %modify Y of the last row to next Y

    gcode0{i}=[gcode0{i}; gcode0{i+1}(1,:) 0]; %repeat next row
    gcode0{i}(end,6)=0.5 * transitionspeed;
    gcode0{i}=[gcode0{i}; gcode0{i+1}(1,:) 10]; %repeat last row and add marker to start jetting: 11=ends, 10=starts

end

%edit last cell
i = length(gcode0);
gcode0{i}=[gcode0{i}; gcode0{i}(end,:)]; %repeat last row
gcode0{i}(end,7)=11; %add marker to start jetting: 11=ends, 10=starts

gcode0{i}=[gcode0{i}; gcode0{i}(end,1:6) 0]; %repeat last row
gcode0{i}(end,3)=Zsafe; %modify Z of the last row to Zsafe

%add marker on second line to start jetting: 11=ends, 10=starts
gcode0{1}=[gcode0{1}(1,:); gcode0{1}(2,:); gcode0{1}(2,1:6) 10; gcode0{1}(2:end,:)];
gcode0{1}(1,3)=Zsafe; %modify Z of the first row to Zsafe


%% Format Gcode
gcodet = [];
if chooseprint==0
    for i = 1:size(gcode0)
        gcodet=[gcodet; gcode0{i}];
    end
else
    gcodet=gcode0{1}; %if printing single electrode
end

gcode0 = string(gcodet);
gcode = "";

for i = 1:size(gcode0) % <----------- adding letters to every line
    if i==1 % <----------- adding G94 and G1 to first line
        gcode(i, 1) = 'G94'; % add tool linear movement in feed rate mode
        gcode(i, 2) = 'G1'; % add tool linear movement in feed rate mode
        gcode(i, 3) = append('X', num2str(gcode0(i, 1)));
        gcode(i, 4) = append('Y', num2str(gcode0(i, 2)));
        gcode(i, 5) = append('Z', num2str(gcode0(i, 3)));
        gcode(i, 6) = append('B', num2str(gcode0(i, 4)));
        gcode(i, 7) = append('C', num2str(gcode0(i, 5)));
        gcode(i, 8) = append('F', num2str(gcode0(i, 6))); % chosen speed
    else
        gcode(i, 1) = append('X', num2str(gcode0(i, 1)));
        gcode(i, 2) = append('Y', num2str(gcode0(i, 2)));
        gcode(i, 3) = append('Z', num2str(gcode0(i, 3)));
        gcode(i, 4) = append('B', num2str(gcode0(i, 4)));
        gcode(i, 5) = append('C', num2str(gcode0(i, 5)));
        gcode(i, 6) = append('F', num2str(gcode0(i, 6))); % chosen speed
        if gcode0(i, 7) == '10' || gcode0(i, 7) == '11' %add marker to start jetting: M11=ends, M10=starts
            gcode(i, 7) = append('M', num2str(gcode0(i, 7)));
        else
            gcode(i, 7) = '';
        end
        gcode(i, 8) = '';

    end
end


%% Add command to turn on/off microjet
% m10 to turn on the relay and m11 to turn it off
g1 = strcat(gcode,{'   '});  %// add whitespacest
%// Convert to char array
outstr_char = char(g1{:});
%// Get size parameters
[m,n] = size(g1);
p = size(outstr_char,2);
%// Reshape + Permute Magic to create a
%// char array "replica" of input cell array
out = reshape(permute(reshape(outstr_char.',p,m,[]),[1 3 2]),n*p,m).';


%% Print g-code
if choosetrace==1%1 for interconnect 2 for electrode
    if chooseprint==0
        writematrix(gcode, append('exportgcode/',subject,'_post/','allinterconnects.txt'), 'Delimiter', 'space')
    else
        writematrix(gcode, append('exportgcode/',subject,'_post/',Names{chooseprint},'interconnect.txt'), 'Delimiter', 'space')
    end
else
    if chooseprint==0
        writematrix(gcode, append('exportgcode/',subject,'_post/','allelectrode.txt'), 'Delimiter', 'space')

    else
        writematrix(gcode, append('exportgcode/',subject,'_post/',Names{chooseprint},'electrode.txt'), 'Delimiter', 'space')
    end
end

%% clean up variables
clear angletemp crossP gode0 p g m n i x y z xOffset yOffset zOffset outstr_char chosenrotx chosenroty npts
%% display electrode names
for i=1:nelect
    Nnames{i,1}=i;
    Nnames{i,2}=Names{i};
end
disp(Nnames)
disp(['chosen electrode:' num2str(chooseprint)])
disp(['number of jets=' num2str(sum(njets))])
disp(['total print time=' num2str(round(sum(timeperelectrode(:,1)),1)) 'min'])
disp(['min total print time (unbounded speed)=' num2str(round(sum(timeperelectrode(:,2)),1)) 'min'])

%% Functions

function c = findCaxisAngle(normV, i)
% x-axis (comparison vector)
cVC = [1 0 0];
% project normal vector to x-y plane
pVC = [normV(i, 1) normV(i, 2) 0];%sign(normV(i,1))*

SignAngle=cross(pVC,cVC);
if SignAngle(3)==0 %if the normal vector is contained in the x-z plane
    c = NaN; %c is undefined, choice of c angle will depend on the c-history
else
    %c = sign(SignAngle(3))*(-90+findAngle(cVC, pVC));%sign(SignAngle(3))*(-90+findAngle(cVC, pVC));
    c = sign(-pVC(2))*(-90+findAngle(cVC, pVC));%sign(SignAngle(3))*(-90+findAngle(cVC, pVC));
end
end
function b = findBaxisAngle(normV, i)
%comparison vector
cVB = normV(i, :);
%projection of normV onto x-y
pVB = [normV(i, 1) normV(i, 2) 0];

%angle correction of -90 degrees
b =  sign(normV(i, 2))*(90-(findAngle(cVB, pVB)));

end
function angle = findAngle(compV, projV)
angle = rad2deg(acos((dot(projV, compV))/(norm(compV)*norm(projV))));
end

function points_global = moveElectrode(data,normalscan,normalm)

%% load points defined in a local coordinate system
points_local = data; % points in local coordinates
%% define plane using normal and local origin
nhat = normalscan; % hat{n} = hat{x}n_x+hat{y}n_y+hat{z}n_z
p0 = [0; 0; 0]; % p_0 = hat{x}x_0+hat{y}y_0+hat{z}z_0 is the origin of the local coordinate system
% the scalar equation of this plane is n_x*(x-x_0)+n_y*(y-y_0)+n_z*(z-z_0)=0
% when the normal is not perpendicular to the yz plane,
% let an arbitrary point on this plane p_1 be
% p_1 = hat{x}x_1+hat{y}y_1+hat{z}z_1
%     = hat{x}x_1+hat{y}(1.1*y0)+hat{z}(1.1*z0)
% then x_1 = x_0-(0.1*n_y*y_0+0.1*n_z*z_0)/n_x
if nhat(1)~=0
    p1 = [p0(1)-(0.1*nhat(2)*p0(2)+0.1*nhat(3)*p0(3))/nhat(1),1.1*p0(2),1.1*p0(3)];
elseif nhat(2)~=0
    p1 = [1.1*p0(1),p0(2)-(0.1*nhat(1)*p0(1)+0.1*nhat(3)*p0(3))/nhat(2),1.1*p0(3)];
elseif nhat(3)~=0
    p1 = [1.1*p0(1),1.1*p0(2),p0(3)-(0.1*nhat(1)*p0(1)+0.1*nhat(2)*p0(2))/nhat(3)];
else
    error('Invalid normal vector nhat = [0,0,0].\n');
end
% define the plane local coordinate system
zhat_prime = nhat/norm(nhat); % let the normal be the local z coordinate
xhat_prime = p1-p0; % the vector pointing to the point we found be \hat{x}^{\prime}
xhat_prime = xhat_prime/norm(xhat_prime); % we normalize because we are civilized
yhat_prime = cross(zhat_prime,xhat_prime); % \hat{y}^{\prime} is orthonormal to both \hat{x}^{\prime} and \hat{y}^{\prime}
yhat_prime = yhat_prime/norm(yhat_prime);
Ttilde = [xhat_prime;yhat_prime;zhat_prime]; % linear map to rotate from local to global
%% transform points coordinates from local to global
points_global = points_local*Ttilde+p0;

end

function plotfig(Interconnect,Electrode,HeadMesh,ps,pm,nelect,Names)

trimesh(HeadMesh,'FaceColor','none','EdgeColor',[0,0,0]+0.8) %plot mesh
axis equal; hold on

for i=1:nelect
    hold on
    plot3(Electrode{i}(:,1),Electrode{i}(:,2),Electrode{i}(:,3),'-b','LineWidth',1)%;quiver3(Electrode{i}(:,1),Electrode{i}(:,2),Electrode{i}(:,3),Electrode{i}(:,4),Electrode{i}(:,5),Electrode{i}(:,6),'AutoScaleFactor',2) %plot normals
end
xlabel('x (mm)'); ylabel('y (mm)'); zlabel('z (mm)'); axis equal

% Scanned Landmarks
plot3(ps(:,1),ps(:,2),ps(:,3),'ok','MarkerSize',15,'LineWidth',2)
text(ps(1,1),ps(1,2),ps(1,3)+15,'p1','FontSize',15)
text(ps(2,1),ps(2,2),ps(2,3)+15,'p2','FontSize',15)
text(ps(3,1),ps(3,2),ps(3,3)+15,'p3','FontSize',15)

% Measured Landmarks
plot3(pm(:,1),pm(:,2),pm(:,3),'*m','MarkerSize',15,'LineWidth',3)
end

%% Plot machine toolpath
function fig4plot(data,g,cAxisAngles,bAxisAngles)
figure(4)
hold on
subplot(2,4,3)
plot3(g(:,1),g(:,2),g(:,3),'-o')
hold on
title('Machine path')
xlabel('x (mm)')
ylabel('y (mm)')
zlabel('z (mm)')
axis equal

subplot(2,4,4)
plot(data(:,1),bAxisAngles,'-o')
hold on
ylabel('B angle (deg)')
xlabel('trace x coordinate (mm)')

subplot(2,4,5)
plot(data(:,1),cAxisAngles,'-o')
hold on
ylabel('C angle (deg)')
xlabel('trace x coordinate (mm)')

subplot(2,4,6)
plot(data(:,1),data(:,1),'-o')
hold on
plot(data(:,1),g(:,1),'-o')
hold on
ylabel('x (mm)')
xlabel('trace x coordinate (mm)')

subplot(2,4,7)
plot(data(:,1),data(:,2),'-o')
hold on
plot(data(:,1),g(:,2),'-o')
ylabel('y (mm)')
xlabel('trace x coordinate (mm)')

subplot(2,4,8)
plot(data(:,1),data(:,3),'-o')
hold on
plot(data(:,1),g(:,3),'-o')
legend('trace','machine')
ylabel('z (mm)')
xlabel('trace x coordinate (mm)')
end