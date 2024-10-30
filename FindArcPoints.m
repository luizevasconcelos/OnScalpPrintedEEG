%% Paolo Giacometti
% Function calculates the coordinate position of points that lie in the arc
% formed by the intersection of a plane and a surface mesh. The plane
% intersecting the mesh is delimited by three fiducial positions. The
% function first calculates the plane intersecting the mesh and all the
% points that delimit the boundary of the intersection. Then the total 
% arc-length is calculated from those points. The points calculated are a
% percentage of the total arc-length distance apart from each other 
% determined by the input Percent. The points are calculated by selecting
% all the points in the arc, computing the segment lengths between arc
% points, and then identifying the arc points that lie in the arc before
% and after the point that would make the arc-length percentage distance
% equal to the point-to-point distance. A unit vector is calculated between
% these two points and also the length remaining to make up the percentage
% length from the first point. Using the unit vector for direction and the
% remaining length for magnitude, each percentage point is found.
% Inputs: Fiducials - [x,y,z] coordinates of 3 fiducials positions that
%         describe the plane.
%         Percent - Percent subdivision desired. 
%         MeshFaces - [:,:,:] indices of all the vertices in each
%         surface mesh face.
%         MeshVertices - [x,y,z] coordinates of all the vertices in a
%         surface mesh.
% Output: ArcPts - [x,y,z] coordinates of points in the arc displaced a
%         percentage distance.

function [ArcPts] = FindArcPoints(Fiducials, Percent, MeshFaces, MeshVertices)

    % Find points in plane intersecting mesh
    EdgePts = MeshPlaneIntersectPoints(Fiducials,MeshFaces,MeshVertices);
    
    % Find sorted points in arc
    if Percent<5 % For total arc 
        [~, ~, arcLength, PtSort] = ComputeArcLength(Fiducials,EdgePts);
    else % For half the arc
        [arcLength, PtSort] = ComputeArcLength(Fiducials,EdgePts);
    end
    % Compute arc segments
    PtSeg = sqrt(sum(diff(PtSort).^2,2)); % Starting from fiducial p1

    % Percentage arc-length
    pal = arcLength*Percent/100;

    % Percentages to compute points
    Perc = pal*(1:100/Percent)';

    % Cumulative sum of segments from points
    CSum = cumsum(PtSeg);

    % Find points that match each percentage subdivision
    PtPerc = zeros(size(Perc,1),3);
    for i=1:size(Perc,1)
        
        % Index of point just prior to percentage point
        Range = length(CSum(CSum<Perc(i)))+1;
        
        % Arc-length from start to range
        SegLength = sum(PtSeg(1:Range-1));
        
        if abs(SegLength-Perc(i))<0.0001
            % If percentage point lies on a mesh node
            PtPerc(i,:) = PtSort(Range,:);
        else
            
            % Arc points before and after percentage point
            EndPoints = PtSort(Range:Range+1,:); 

            % Unit vector from endpoints
            UVec = diff(EndPoints)/sqrt(sum(diff(EndPoints).^2));

            % Length from point before percentage point to percentage point
            RemLength = Perc(i)-SegLength;

            % Point position that matches percentage distance desired
            PtPerc(i,:) = EndPoints(1,:) + UVec*RemLength;
        end

    end

    % Add initial point to list
    ArcPts = [Fiducials(1,:); PtPerc];
end


%% Paolo Giacometti
% Function calculates the arc-length from a set of points along a surface.
% First, it sorts all the points so that they are ordered starting from the
% first fiducial position to the final fiducial position. It also adds the
% fiducial points to the list if they are not already there. Then, it 
% calculates the distance between all points in the arc and measures them 
% to get the total arc-length. Then, it cuts the data set to only include 
% points in plane that are within the arc. Finally, the function divides 
% the arc-length in the section before and after the middle fiducial
% position. 
% Inputs: Fiducials - [x,y,z] coordinates of 3 fiducials positions that
%         describe the plane.
%         EdgePts - [x,y,z] coordinates of points in contact to a plane 
%         intersecting a mesh.
% Outputs: arcLengthH - [d] half arc-length measured.
%          PtSortH - [x,y,z] coordinates of points in half the arc sorted 
%          from the first fiducial position to the second.
%          arcLength - [d] total arc-length measured.
%          PtSort - [x,y,z] coordinates of points in the arc sorted from
%          the first fiducial position to the last.

function [arcLengthH, PtSortH, arcLength, PtSort] = ComputeArcLength(Fiducials,EdgePts)


    % Fiducials - p1 - p2 - p3
    p1 = Fiducials(1,:);
    p2 = Fiducials(2,:);
    p3 = Fiducials(3,:);
          
    % Add fiducials to points
    [PtSort p2idx] = AddFiducials(Fiducials, EdgePts);
    
    % Compute arc segments
    PtSeg1 = sqrt(sum(diff(PtSort(1:p2idx,:)).^2,2));
    PtSeg2 = sqrt(sum(diff(PtSort(p2idx:end,:)).^2,2));
    PtSeg = [PtSeg1;PtSeg2];
     
    % Total arc-length
    arcLength = sum(PtSeg);
    % Half arc-length (from fiducial 1 to fiducial 2)
    arcLengthH = sum(PtSeg1);

    % Select points for half the arc (from fiducial 1 to fiducial 2)
    PtSortH = PtSort(1:p2idx,:);
    
end

%% Paolo Giacometti
% Function to add fiducial points in order of position given a sorted set
% of points. It also selects only the points present in the arc and removes
% all the extra points in the mesh edge.
% Inputs: Fiducials [x,y,z] coordinates of the fiducial positions
%         Points [x,y,z] coordinates of all other points to be sorted
% Output: PSort [x,y,z] coordinates of all points including the reference
%         point sorted in positional order starting with fiducial 1 and
%         ending in fiducial 3.
%         fid2 - index of middle point fiducial in Points
function [PSort,fid2] = AddFiducials(Fiducials, Points)

    % Check if fiducials are on the list of points   
    fid1 = find(sum(Points==repmat(Fiducials(1,:),size(Points,1),1),2)==3);
    fid2 = find(sum(Points==repmat(Fiducials(2,:),size(Points,1),1),2)==3);
    fid3 = find(sum(Points==repmat(Fiducials(3,:),size(Points,1),1),2)==3);
    
    % Add fiducials to list of points
    if isempty(fid1)
        Points = [Fiducials(1,:); Points];
        fid1 = 1;
    end
    
    if isempty(fid2)
        % Find points in Points closest to second fiducial
        DistP2 = sqrt(sum((Points-repmat(Fiducials(2,:),size(Points,1),1)).^2,2));
        mpt = min(DistP2);
        
        % Find indices of that point
        mpidx = find(DistP2==mpt);
        
        % Find point closest to that one
        mp2idx = find(DistP2==min([DistP2(mpidx-1) DistP2(mpidx+1)]));
        mp2idx = mp2idx(1);
        
        if mp2idx>mpidx
            Points = [Points(1:mpidx,:);Fiducials(2,:);Points(mpidx+1:end,:)];
        else
            Points = [Points(1:mpidx-1,:);Fiducials(2,:);Points(mpidx:end,:)];
        end
        
        % Get index of middle fiducial point in Points
        fid2 = find(sum(Points==repmat(Fiducials(2,:),size(Points,1),1),2)==3);
    end   
    
    if isempty(fid3)
        % Find points in Points closest to third fiducial
        DistP3 = sqrt(sum((Points-repmat(Fiducials(3,:),size(Points,1),1)).^2,2));
        mpt = min(DistP3);
        
        % Find indices of that point
        mpidx = find(DistP3==mpt);
        
        % Find point closest to that one
        mp3idx = find(DistP3==min([DistP3(mpidx-1) DistP3(mpidx+1)]));
        mp3idx = mp3idx(1);
        
        if mp3idx>mpidx
            Points = [Points(1:mpidx,:);Fiducials(3,:);Points(mpidx+1:end,:)];
        else
            Points = [Points(1:mpidx-1,:);Fiducials(3,:);Points(mpidx:end,:)];
        end
        fid3 = find(sum(Points==repmat(Fiducials(3,:),size(Points,1),1),2)==3);
    end
    
    % Remove points not present in the fiducial-to-fiducial arc
    if fid3>fid2
        PSort = Points(fid1:fid3,:);
    else
        resort = [1 size(Points,1):-1:2];
        PointsRS = Points(resort,:);
        fid3 = find(sum(PointsRS==repmat(Fiducials(3,:),size(PointsRS,1),1),2)==3);
        fid2 = find(sum(PointsRS==repmat(Fiducials(2,:),size(PointsRS,1),1),2)==3);
        PSort = PointsRS(fid1:fid3,:);
    end
        
end









