
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                                                           %
% This code generates a personalized EEG electrode and interconnect layout   %                                                                            % 
% using a 3D scan of the head (slt file) for a selected EEG locations within % 
% the 10-20 system and its deviratives                                       % 
% Authors: Sandhya Tiku, Luize Vasconcelos                                   %
% Nanshu Lu Research Group, UT austin                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; clear; clc

%% Chose # electrodes to display
nelectdisplay=15;

%% Chose to replace 4 of 10-20 positions by 4 5-5 positions
disp329 = false; % display 5-5 positions with text
choose_replace = false; % whether or not to replace the 4 10-20 positions
EEGLab329_chosen = ["P7", "P8", "O2","FCz"];
EEGLab_replace = ["P7", "P8", "O2", "Fz"];

%% Choose # smoothing steps
nNormalSmoothrounds=40;

%% Chose electrode diameter (mm)
areaelectrode=1.1; %cm^2
diameter= 2*sqrt(areaelectrode/pi)*10; %mm

%% Select Fiducials [x,y,z] [Nz;Iz;M1;M2] coordinates
% nasion, the inion, and the left and right preauricular points
filename= 'hair_mannequin.stl';
ExportFolder='hair_mannequin';
if ismac
    data = stlread(append('Head mesh/',filename));
elseif ispc
    data = stlread(append('Head mesh\',filename));
end

Faces = data.ConnectivityList; %faces
System = 3; % 3:10-20, 2:10-10 1: 10-5
CoordFlag = 0; % 1: align coordinate system to fiducials 0: do nothing

%% Calculate the Landmarks and rotate the head the needed degrees
[pt1, pt2, pt3, LandmarkNames] = calc_landmarks(filename);
[Rx, Ry, Rz] = rotate_head(filename);

%% Rotate the whole head including the landmarks positions
Vertices_transposed = data.Points'; % call vertices coordinates
Vertices = Ry*Rz*Rx*Vertices_transposed; % rotate vertice coordinates
Vertices = Vertices'; clear Vertices_transposed 
%rotate landmarks
pt1 = Ry*Rz*Rx*pt1; 
pt2 = Ry*Rz*Rx*pt2; 
pt3 = Ry*Rz*Rx*pt3; 

%% Creates the mesh with the transposed and rotated vertices and the Connectivity List
dataref = triangulation(Faces, Vertices); %create mesh

Fiducials = def_fiducials(filename); % Input fiducials through function
Nasion = Fiducials(1,:); Inion = Fiducials(2,:); par = Fiducials(3,:); pal = Fiducials(4,:);

%% plots the mesh
plot_mesh(dataref, Fiducials)

%% Runs the function that computes the EEG pts or the 10-20 positions
[EEGPts,EEGLab,VerticesT] = ComputeEEGPos(Fiducials,Faces,Vertices,System,CoordFlag);
nelect = length(EEGPts);

%% Set the new origin as the Nasion
datarefpts(:,1) = dataref.Points(:,1) - Nasion(1);%-(-1.36);
datarefpts(:,2) = dataref.Points(:,2) - Nasion(2);%-13.89;
datarefpts(:,3) = dataref.Points(:,3) - Nasion(3);%-107.4;

%% Choose the close by positions to replace the 10-20 positions
if (choose_replace == true)
    [EEGPts329,EEGLab329,~] = ComputeEEGPos(Fiducials,Faces,Vertices,1,CoordFlag);
    [EEGPts, EEGLab] = choose_329traces(EEGPts, EEGLab, EEGLab329_chosen, EEGLab_replace, EEGPts329, EEGLab329);
end
%% Create Interconnect path
[alltraces, Ielectorder] = TraceFinder(dataref, EEGPts,nelectdisplay,Fiducials);

%% Create electrode path
[allelectrodes]=ElectrodeFinder(alltraces,nelectdisplay,diameter);

%% plots all the interconnects and electrodes in Figure 3
figure(3)
trimesh(dataref,'FaceColor','w','EdgeColor',[0,0,0]+0.8) %plot mesh
axis equal; hold on; xlabel('x(mm)'); ylabel('y(mm)'); zlabel('z(mm)'); grid off
for i = 1:nelectdisplay
    plot3(alltraces{i}(:,1),alltraces{i}(:,2),alltraces{i}(:,3),'-m','LineWidth',5) %plot interconnects
    plot3(allelectrodes{i}(:,1),allelectrodes{i}(:,2),allelectrodes{i}(:,3),'-b','LineWidth',5) %plot electrodes
end
title('custom electrode (blue) and interconnect (magenta) layout')
%% If disp329 is true then it displays 5-5 positions with text labels
if (disp329 == true)
    figure(3);
    text(EEGPts329(:,1),EEGPts329(:,2),EEGPts329(:,3),EEGLab329,'Color','black','FontSize',12);
end

%% Prepare for postprocessor 
lengthtrace = zeros(nelect,1);
PathNames = EEGLab(Ielectorder); 
if nelectdisplay == 15
    PathNames=PathNames(3:17,1);
end
for i = 1:nelectdisplay
    %% Define new origin at the electrode 
    trace = alltraces{i};
    tracenormals = alltraces{i}(:,4:6);

    %% Make all normals point outwards
    % rotate the vectors normal to the mesh
    norm_rot = rotx(45)*tracenormals'; %rotate normals by head angle on the chair
    norm_rot = norm_rot';

    % the orientation of the normal vectors is random (comes from the order of the cross product)
    % For the head position in the chair, we expect all normals to be point
    % upwards (i.e., jet happens at most at 90deg angle, not more. Hence, we
    % transform the orientation of the normal to be always up:
    for j = 1:length(norm_rot) % correct signs according to the z coord is +
        if norm_rot(j,3) < 0
            norm_rot(j,:) = -1* norm_rot(j,:); % if z is negative, multiply entire vector by -1
        end
    end

    %rotate normals back
    norm_rot = rotx(-45)*norm_rot'; %rotate normals by head angle on the chair
    tracenormalsup = norm_rot';

    %% make all normals point outwards
    %related to side normals
    % T3
    if i == nelectdisplay || i == nelectdisplay-1 || i == nelectdisplay-2
        for j = 1:length(tracenormalsup)
            if sign(tracenormalsup(j,1)) == 1
            else
                tracenormalsup(j,:) =- tracenormalsup(j,:);
            end

        end
    end
    %T4
    if i==1 || i == 2 || i == 3
        for j = 1:length(tracenormalsup)
            if sign(tracenormalsup(j,1)) == 1
                tracenormalsup(j,:) =- tracenormalsup(j,:);
            end

        end
    end
    
    for j=1:length(tracenormalsup)
        tracenormalsup(j,:) = tracenormalsup(j,:)/norm(tracenormalsup(j,:));%normalize trace vectors
    end
    tracenormalup_smooth = tracenormalsup;
    for k = 1:nNormalSmoothrounds % repeat smoothing step as many times as needed
        %smoothout trace normals through moving average
        tracenormalup_smooth(:,1) = movmean(tracenormalup_smooth(:,1),3); %average x coords
        tracenormalup_smooth(:,2) = movmean(tracenormalup_smooth(:,2),3); %average y coords
        tracenormalup_smooth(:,3) = movmean(tracenormalup_smooth(:,3),3); %average z coords
    end
    %Plot and export interconnects
    PathName = PathNames{i};
    
    %% add labels to electrodes
    text(trace(1,1), trace(1,2), trace(1,3) + 10,PathName, 'Color','red','FontSize',12)
    allinterconnects{i,1} = [trace(:,1:3) tracenormalup_smooth];
    
    clear tracenormalup_smooth
    figure(10)
    hold on
    %% plots all the interconnects and electrodes in Figure 3
    if i==1
    trimesh(dataref,'FaceColor','none','EdgeColor',[0,0,0]+0.8) %plot mesh
    axis equal; hold on; xlabel('x(mm)'); ylabel('y(mm)'); zlabel('z(mm)'); grid off  
    end
    plot3(allinterconnects{i}(:,1),allinterconnects{i}(:,2),allinterconnects{i}(:,3),'-m','LineWidth',5) %plot interconnects
    quiver3(allinterconnects{i}(:,1),allinterconnects{i}(:,2),allinterconnects{i}(:,3),allinterconnects{i}(:,4),allinterconnects{i}(:,5),allinterconnects{i}(:,6),'AutoScaleFactor',0.5) %plot normals

end

%% Plot and  export electrodes
electrodesexport = cell(nelectdisplay,1);
for i = 1:nelectdisplay % prepare date for exporting
    alltraces{i}(1,4:6) = alltraces{i}(1,4:6) / norm(alltraces{i}(1,4:6)); %creates normal vectors for electrode path based on the first row of the interconnect's normal
    electnormal = repmat(alltraces{i}(1,4:6), length(allelectrodes{i}), 1);  
    electrodesexport{i} = [allelectrodes{i} electnormal];
end
title('Interconnect traces and normal vectors')
%pose matching landmarks
% plot3(pt1(1),pt1(2),pt1(3),'or','MarkerSize',20,'LineWidth',2)  
% plot3(pt2(1),pt2(2),pt2(3),'*r','MarkerSize',20,'LineWidth',2)
% plot3(pt3(1),pt3(2),pt3(3),'*r','MarkerSize',20,'LineWidth',2)
%text(pt1(1),pt1(2),pt1(3)+10,'pt1','Color','red','FontSize',16)
%text(pt2(1),pt2(2),pt2(3)+10,'pt2','Color','red','FontSize',16)
%text(pt3(1),pt3(2),pt3(3)+10,'pt3','Color','red','FontSize',16)

Landmarks=[pt1';pt2';pt3';];
InterconnectElectrodePaths= [allinterconnects electrodesexport PathNames];

if strcmp(ExportFolder,'no')
else
    save(['Postprocessor/',ExportFolder,'/InterconnectElectrodePaths.mat'],'InterconnectElectrodePaths')
    save(['Postprocessor/',ExportFolder,'/HeadMesh.mat'],'dataref')
    save(['Postprocessor/',ExportFolder,'/Landmarks.mat'],'Landmarks')
    save(['Postprocessor/',ExportFolder,'/LandmarkNames.mat'],'LandmarkNames')
end

function plot_mesh(dataref, Fiducials)
    trimesh(dataref,'FaceColor','none','EdgeColor',[0,0,0]+0.8) %plot mesh
    hold on
    axis equal
    hold on
    plot3(Fiducials(1,1), Fiducials(1,2), Fiducials(1,3),'oy','MarkerSize',20,'LineWidth',2)
    plot3(Fiducials(2,1), Fiducials(2,2), Fiducials(2,3),'oy','MarkerSize',20,'LineWidth',2)
    plot3(Fiducials(3,1), Fiducials(3,2), Fiducials(3,3),'oy','MarkerSize',20,'LineWidth',2)
    plot3(Fiducials(4,1), Fiducials(4,2), Fiducials(4,3),'oy','MarkerSize',20,'LineWidth',2)
    title('Face fiducials')

end

function [EEGPts, EEGLab] = choose_329traces(EEGPts, EEGLab, EEGLab329_chosen, EEGLab_replace, EEGPts329, EEGLab329)

    EEGLab329_chosen_pts = [1:3; 1:3; 1:3; 1:3];
    for i = 1:length(EEGLab329)
        if (strcmp(EEGLab329(i), EEGLab329_chosen(1)))
            EEGLab329_chosen_pts(1,:) = EEGPts329(i,:);
        elseif (strcmp(EEGLab329(i), EEGLab329_chosen(2)))
            EEGLab329_chosen_pts(2,:) = EEGPts329(i,:);
        elseif (strcmp(EEGLab329(i), EEGLab329_chosen(3)))
            EEGLab329_chosen_pts(3,:) = EEGPts329(i,:);
        elseif (strcmp(EEGLab329(i), EEGLab329_chosen(4)))
            EEGLab329_chosen_pts(4,:) = EEGPts329(i,:);
        end
    end
    
    for i = 1:length(EEGLab)
        if (strcmp(EEGLab(i), EEGLab_replace(1)))
            EEGPts(i,:) = EEGLab329_chosen_pts(1,:);
            EEGLab(i,1) = cellstr(EEGLab329_chosen(1,1));
        elseif (strcmp(EEGLab(i), EEGLab_replace(2)))
            EEGPts(i,:) = EEGLab329_chosen_pts(2,:);
            EEGLab(i,1) = cellstr(EEGLab329_chosen(1,2));
        elseif (strcmp(EEGLab(i), EEGLab_replace(3)))
            EEGPts(i,:) = EEGLab329_chosen_pts(3,:);
            EEGLab(i,1) = cellstr(EEGLab329_chosen(1,3));
        elseif (strcmp(EEGLab(i), EEGLab_replace(4)))
            EEGPts(i,:) = EEGLab329_chosen_pts(4,:);
            EEGLab(i,1) = cellstr(EEGLab329_chosen(1,4));
        end
    end
end

function Fiducials = def_fiducials(filename)
    %% Define markers
    if strcmp(filename,'hair_mannequin.stl')
        Nasion = [6.79 -57.26 0.0]; 
        Inion = [6.79 101.63 -28.22];
        par = [-67.08 14.63 0.0];
        pal =[78.5 13.24 -0.1];
    end

    % Defines an array for all the fiducials in the specific order of Nasion,
    % Inion, left and right pre-auricular points
    Fiducials = [Nasion; Inion; pal; par];
end

function [Rx,Ry,Rz] = rotate_head(filename)
    %% ROTATE HEAD
    if strcmp(filename,'hair_mannequin.stl') % 2022/08/25
        % If STL file is not oriented correctly, define angles to rotate in order to get the head to be face forward
        Rz = rotz(0);
        Ry = roty(0);
        Rx = rotx(0);
    end 
end

function [pt1, pt2, pt3, LandmarkNames] = calc_landmarks(filename)
    %% Input Landmarks
    if strcmp(filename,'hair_mannequin.stl') %2022/08/30
        pt1 = [10.05 ; 61.49 ; 110.83];
        pt2 = [-35.57 ; 105.14 ; 45.64];
        pt3 = [-18.62 ; 109.23 ; -1.078];   
        LandmarkNames=["top"; "right ear"; "bottom"];
    end
end

