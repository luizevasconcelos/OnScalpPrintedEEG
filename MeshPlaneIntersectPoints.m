%% Paolo Giacometti
% % Function to take 3 points and compute a plane that passes through a mesh.
% % Then it computes all the points in the mesh intersecting the plane. First
% % it computes the plane equation. Then it find all the mesh faces that
% % intersect the plane. From this, there are 3 options of intersecion: 1) A
% % triangular face with a singleton point to one side of the plane and two
% % other points to the other side of the plane. 2) A triangular face with
% % one point in the plane, one point to one side of the plane and the third
% % point to the other side of the plane. 3) A triangular face with two point
% % in the plane and one point to either side of the plane. First the
% % function deals with option (1). It takes the points from those faces and
% % deduces the singleton point that lies to one side of the plane and the
% % pair of points that lies to the other side. From this, the function
% % calculates the intersection point between each pair point and the
% % singleton point. Similarly, for option (2), the intersection point 
% % between the two points at either side is computed, and then, along the 
% % ones that lie on the plane are added to the set. Finally, for option (3)
% % the points that lie exactly on the plane are added to the set. [There is
% % an option (4) where there is a single point that lies on the plane and
% % two points that are both at the same side of the plane. These are ignored
% % because they share the point that lies in the plane with faces in one of
% % the three other cases so they are already accounted for.]
% % Inputs: Fiducial - [x,y,z] coordinates of 3 fiducials positions that
% %         describe the plane.
% %         MeshFaces - [:,:,:] indices of all the vertices in each
% %         surface mesh face.
% %         MeshVertices - [x,y,z] coordinates of all the vertices in a
% %         surface mesh.
% % Output: EdgePts - [x,y,z] coordinates of points in the plane intersecting
% %         the mesh. These points are sorted positionally starting from
% %         fiducial point p1.
% % 
% % Example:
% Fiducial = [-10 -10 0;-10 10 0;10 10 0];
% MeshVertices = [-4 2 -3;-4 2 3;-2 2 -3;0 2 3;2 2 -3;4 2 3;4 2 -3;
%                 -4 -2 -3;-4 -2 3;-2 -2 -3;0 -2 3;2 -2 -3;4 -2 3;4 -2 -3];
% MeshFaces = [1 2 3; 2 3 4; 3 4 5; 4 5 6; 5 6 7; 8 9 10; 9 10 11; 10 11 12; 
%              11 12 13; 12 13 14; 1 8 9; 1 2 9; 7 13 14; 6 13 7];
% fv.Faces = MeshFaces;
% fv.Vertices = MeshVertices;
% figure
% plot3(MeshVertices(:,1),MeshVertices(:,2),MeshVertices(:,3),'.')
% hold on
% patch(fv,'FaceColor','w')
% [x,y]=meshgrid(-10:10);
% surf(x,y,zeros(21,21))
% EdgePts = MeshPlaneIntersectPoints(Fiducial, MeshFaces, MeshVertices);
% plot3(EdgePts(:,1),EdgePts(:,2),EdgePts(:,3),'r*')

function [EdgePts] = MeshPlaneIntersectPoints(Fiducial, MeshFaces, MeshVertices)

    % Compute plane that passes p1 - p2 - p3
    % Plane equation ax+by+cz+d=0
    % 3 point plane calculation 
    p1 = Fiducial(1,:);
    p2 = Fiducial(2,:);
    p3 = Fiducial(3,:);

    % Compute constants for plane equation
    a = (p2(1,2)-p1(1,2))*(p3(1,3)-p1(1,3))-(p3(1,2)-p1(1,2))*(p2(1,3)-p1(1,3));
    b = (p2(1,3)-p1(1,3))*(p3(1,1)-p1(1,1))-(p3(1,3)-p1(1,3))*(p2(1,1)-p1(1,1));
    c = (p2(1,1)-p1(1,1))*(p3(1,2)-p1(1,2))-(p3(1,1)-p1(1,1))*(p2(1,2)-p1(1,2));
    d = -(a*p1(1,1)+b*p1(1,2)+c*p1(1,3)); % This would work with p1, p2, or p3

    % Compute normalized normal to plane
    np = [a b c]/sqrt(a^2+b^2+c^2);
    
    % Compute distance to all mesh vertices
    Dist = (a*MeshVertices(:,1)+b*MeshVertices(:,2)+c*MeshVertices(:,3)+d)/sqrt(a^2+b^2+c^2);
    
    % Find sets of distances on mesh faces
    DistFace = Dist(MeshFaces);
    
    % Tabulate distances from plane (>0 in front, <0 behind, =0 on plane)
    LogicDist = (DistFace>0)-(DistFace<0);
    
    % Find which distances intersect the plane
    InterFaceDist = DistFace(abs(sum(LogicDist,2))==1&sum(abs(LogicDist),2)==3,:);
    
    % Find vertices of those faces
    InterFaceVert = MeshFaces(abs(sum(LogicDist,2))==1&sum(abs(LogicDist),2)==3,:);
    
    % Tabulate distances to plane only of chosen faces
    InterFaceLogic = LogicDist(abs(sum(LogicDist,2))==1&sum(abs(LogicDist),2)==3,:);
    
    % Find vertex indices of the single points at one side of the plane
    SingleLogic = InterFaceLogic;
    SingleLogic(sum(SingleLogic,2)<0,:)=-SingleLogic(sum(SingleLogic,2)<0,:);
    SingleLogic=SingleLogic<0;
    
    % Find vertices indices of the single points at one side of the plane
    IFV = InterFaceVert'; IFV = IFV(:);
    SL = SingleLogic'; SL = SL(:);
    SingleVertsidx = IFV(SL);
    
    % Find vertices of the single points at one side of the plane
    SingleVerts = MeshVertices(SingleVertsidx,:);
    
    % Find vertex indices of the pair points at other side of the plane
    PairLogic = ~SingleLogic;
    
    % Find vertices indices of the pair points at other side of the plane
    PL = PairLogic'; PL = PL(:);
    PairVertsidx = IFV(PL);
    PairVertsidx = [PairVertsidx(1:2:end) PairVertsidx(2:2:end)];

    % Find vertices of the pair points at other side of the plane
    PairVerts = MeshVertices(PairVertsidx(:,1),:);
    PairVerts(:,:,2) = MeshVertices(PairVertsidx(:,2),:);

    % Compute intersection point constant for P = P1+u(P2-P1) where
    % u=(ax1+by1+cz1+d)/(a(x1-x2)+b(y1-y2)+c(z1-z2)) from point pairs 
    % behind and in front of plane
    U = (a*SingleVerts(:,1)+b*SingleVerts(:,2)+c*SingleVerts(:,3)...
        +d)./(a*(SingleVerts(:,1)-PairVerts(:,1,1))+b*(SingleVerts(:,2)...
        -PairVerts(:,2,1))+c*(SingleVerts(:,3)-PairVerts(:,3,1)));

    U(:,2) = (a*SingleVerts(:,1)+b*SingleVerts(:,2)+c*SingleVerts(:,3)...
        +d)./(a*(SingleVerts(:,1)-PairVerts(:,1,2))+b*(SingleVerts(:,2)...
        -PairVerts(:,2,2))+c*(SingleVerts(:,3)-PairVerts(:,3,2))); 
    
    % Compute edge points from point pairs behind and in front of plane for
    % P = P1+u(P2-P1)
    InterPts = SingleVerts+repmat(U(:,1),1,3).*(PairVerts(:,:,1)-SingleVerts);
    InterPts(:,:,2) = SingleVerts+repmat(U(:,2),1,3).*(PairVerts(:,:,2)-SingleVerts);

    % Check that meshing errors havent created a pair of equal points    
    InterPts = round(InterPts*10000)/10000; % Remove computation error
    InterPts = InterPts(sum((InterPts(:,:,1)-InterPts(:,:,2))==0,2)~=3,:,:);    
    
    % Compute list of point indices that form the edge
    Edges = [(1:size(InterPts,1))' (size(InterPts,1)+1:2*size(InterPts,1))'];

    % Concatenate all interface points     
    InterPts = [InterPts(:,:,1);InterPts(:,:,2)];
    
    % Find faces with nodes that lie directly on plane (pick only those to
    % one side of the plane to avoid repeating edges)
    NodeFaceLogic = LogicDist((sum(LogicDist,2))==1&sum(abs(LogicDist),2)==1,:);
    if isempty(NodeFaceLogic)
    else
        % Find logical for points in the face lying on plane
        NodeFaceLogic = NodeFaceLogic<1;
        
        % Find indices of vertices of those faces (pick only those to one
        % side of the plane to avoid repeating edges)
        NodeFaceVert = MeshFaces((sum(LogicDist,2))==1&sum(abs(LogicDist),2)==1,:);
        
        % Find indices of vertices of the points in the plane
        NFV = NodeFaceVert'; NFV = NFV(:);
        NFL = NodeFaceLogic'; NFL = NFL(:);
        NodeVertsidx = NFV(NFL);
    
        % Find node coordinates of these points
        NodeVerts = MeshVertices(NodeVertsidx,:);

        % Compute list of point indices that form the node edges
        NodeEdge = [(1:2:size(NodeVertsidx,1))' (2:2:size(NodeVertsidx,1))'];
        
        % Concatenate edge points
        InterPts = [InterPts;NodeVerts];
        
        % Concatenate Edges to include these edges
        Edges = [Edges; NodeEdge+max(Edges(:))];
    end
    
    % Find faces with nodes that intersect the plane and have a zero (have
    % all 3 possibilities in one face (>0,<0,=0)
    % In this case the pair of points are at either side of the plane.
    % Need to identify these points (>0,<0), compute the intersection of
    % the segment between these points and the plane, and use the
    % intersection point as an edge point to the single point selected
    % that is already on plane (=0).      
    InterSingleLogic = LogicDist(abs(sum(LogicDist,2))==0&sum(abs(LogicDist),2)==2,:);
    if isempty(InterSingleLogic)
    else
        % Find logical for points in the face lying on plane
        InterSingleLogic = abs(InterSingleLogic);
        InterSingleLogic = InterSingleLogic<1;
        
        % Find indices of vertices of those faces
        InterNodeFaceVert = MeshFaces(abs(sum(LogicDist,2))==0&sum(abs(LogicDist),2)==2,:);
        
        % Find indices of vertices of the point in the plane
        INFV = InterNodeFaceVert'; INFV = INFV(:);
        ISL = InterSingleLogic'; ISL = ISL(:);
        InterNodeFaceVertsidx = INFV(ISL);
        
        % Find node coordinates of these points on plane
        InterSingleVerts = MeshVertices(InterNodeFaceVertsidx,:);
        
        % Find vertex indices of the pair points at either side of the plane
        InterPairLogic = ~InterSingleLogic;

        % Find vertices indices of the pair points at either side of the plane
        IPL = InterPairLogic'; IPL = IPL(:);
        InterPairVertsidx = INFV(IPL);
        InterPairVertsidx = [InterPairVertsidx(1:2:end) InterPairVertsidx(2:2:end)];   

        % Find vertices of the pair points at either side of the plane
        InterPairVerts = MeshVertices(InterPairVertsidx(:,1),:);
        InterPairVerts(:,:,2) = MeshVertices(InterPairVertsidx(:,2),:);

        % Compute intersection point constant for P = P1+u(P2-P1) where
        % u=(ax1+by1+cz1+d)/(a(x1-x2)+b(y1-y2)+c(z1-z2)) from point pairs 
        % behind and in front of plane
        Uo = (a*InterPairVerts(:,1,1)+b*InterPairVerts(:,2,1)+c*InterPairVerts(:,3,1)...
            +d)./(a*(InterPairVerts(:,1,1)-InterPairVerts(:,1,2))+b*(InterPairVerts(:,2,1)...
            -InterPairVerts(:,2,2))+c*(InterPairVerts(:,3,1)-InterPairVerts(:,3,2)));

        % Compute edge points from point pair behind and in front of plane for
        % P = P1+u(P2-P1)
        InterNodePts = InterPairVerts(:,:,1)+repmat(Uo,1,3).*(InterPairVerts(:,:,2)-InterPairVerts(:,:,1));        
        
        % Concatenate intersecting points with points in plane
        InterFacePts = [InterSingleVerts;InterNodePts];
        
        % Compute list of point indices that form the node edges
        InterNodeEdge = [(1:size(InterFacePts,1)/2)' (size(InterFacePts,1)/2+1:size(InterFacePts,1))'];
        
        % Concatenate edge points
        InterPts = [InterPts;InterFacePts];
        
        % Concatenate Edges to include these edges
        Edges = [Edges; InterNodeEdge+max(Edges(:))];
    end    

    % Round off points to 4th decimal point to ignore computational errors
    InterPts = round(InterPts*10000)/10000;

    % Compute list of points in plane positionally sorted
    EdgePts = SortEdgeNodes(Edges,InterPts,Fiducial);
    

end
 