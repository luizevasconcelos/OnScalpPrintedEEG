%% Paolo Giacometti
% Compute EEG positions given the faces and vertices of a mesh and the
% fiducial positions to base the positions off of. First, it makes a guess
% of where Cz position is based on the fiducials, and runs a minimization
% algorithm to optimize the position of Cz. Once it is found, it makes
% planes based on the fiducials and Cz in the same way as it is done in
% practice. One the main arcs are computed, arcs based on those computed
% points are also calculated to make up the entire set of EEG positions.
% Inputs: Fiducials [x,y,z] [Nz;Iz;M1;M2] coordinates of the 4 fiducials of
%         the head: the nasion, the inion, and the left and right 
%         preauricular points, in that order.
%         Faces [p1,p2,p3] indices of the points that make up each mesh
%         face.
%         Vertices [x,y,z] coordinates of all vertex points of the mesh.
%         System [1,2,3] the desired EEG system:
%                1) 10-5
%                2) 10-10
%                3) 10-20
%                4) 10-10 EEG and modified 10-10 NIRS system
%         CoordFlag [0,1] flag to choose to transform all coordinates to
%         align coordinate system with fiducials. 
% Outputs: EEGPos [x,y,z] coordinates of all points matching EEGLab.
%          EEGLab {lab} Cell array of names of EEG positions.
%          VerticesT [x,y,z] coordinates of all vertices of the mesh
%          transformed such that the origin lies at the midpoint of M1 and
%          M2 and Nz is oriented in the negative x-direction.

function [EEGPts,EEGLab,VerticesT] = ComputeEEGPos(Fiducials,Faces,Vertices,System,CoordFlag)
    
    
    % Transform all coordinates to match fiducials with axes where Nz
    % points in the negative x-direction and the origin lies in the
    % midpoint between M1 and M2.
    if CoordFlag
        % Transform vertices
        VerticesT = TransformCoords(Vertices,Fiducials(1,:),Fiducials(3,:),Fiducials(4,:));

        % Transform fiducials
        FiducialsT = TransformCoords(Fiducials,Fiducials(1,:),Fiducials(3,:),Fiducials(4,:));
    else
        VerticesT = Vertices;
        FiducialsT = Fiducials;
    end
    
    % Normalize Vertices to range from [0-100] to facilitate computation
    % and minimize error.
    Vmin = abs(min(VerticesT(:)));
    Vmax = max(Vertices(:));
    VerticesN = 100*(VerticesT + Vmin)/(Vmax+Vmin);
    FiducialsN = 100*(FiducialsT + Vmin)/(Vmax+Vmin);

    
    % Convert tetrahedral mesh into triangular mesh
    if size(Faces,2)==4        
        TRI = TriRep(Faces,VerticesN);
        FacesT = freeBoundary(TRI);
    else
        FacesT = Faces;       
    end
    
    % Fiducial positions
    Nz = FiducialsN(1,:); % Nasion
    Iz = FiducialsN(2,:); % Inion
    M1 = FiducialsN(3,:); % Left-Preauricular
    M2 = FiducialsN(4,:); % Right-Preauricular

    % Guess Cz (largest point in z in the mesh positioned in the x-y
    % average of M1 and M2)
    [~,CzIdx]=max(VerticesN(:,3));
    X = [mean([M1(1:2);M2(1:2)],1) VerticesN(CzIdx,3)];

    % Optimize Cz
    [Cz] = OptimizeCz(X, FiducialsN, FacesT, VerticesN);

    %% Compute Arc Points

    % Compute points from plane that passes Nz - Cz - Iz
    INPts = FindArcPoints([Nz;Cz;Iz],10,FacesT,VerticesN);
    INPts2 = FindArcPoints([Iz;Cz;Nz],10,FacesT,VerticesN);

    NFpz = INPts(2,:);
    Fpz = INPts(3,:);
    AFpz = INPts(4,:);
    AFz = INPts(5,:);
    AFFz = INPts(6,:);
    Fz = INPts(7,:);
    FFCz = INPts(8,:);
    FCz = INPts(9,:);
    FCCz = INPts(10,:);
    CCPz = INPts2(10,:);
    CPz = INPts2(9,:);
    CPPz = INPts2(8,:);
    Pz = INPts2(7,:);
    PPOz = INPts2(6,:);
    POz = INPts2(5,:);
    POOz = INPts2(4,:);
    Oz = INPts2(3,:);
    OIz = INPts2(2,:);
    
    % Compute points from plane that passes M1(T9) - Cz - M2(T10)
    PPts = FindArcPoints([M1;Cz;M2],10,FacesT,VerticesN);
    PPts2 = FindArcPoints([M2;Cz;M1],10,FacesT,VerticesN);
    
    T9h = PPts(2,:);
    T7 = PPts(3,:);
    T7h = PPts(4,:);
    C5 = PPts(5,:);
    C5h = PPts(6,:);
    C3 = PPts(7,:);
    C3h = PPts(8,:);
    C1 = PPts(9,:);
    C1h = PPts(10,:);
%    Cz = PPts(11,:);
    C2h = PPts2(10,:);
    C2 = PPts2(9,:);
    C4h = PPts2(8,:);
    C4 = PPts2(7,:);
    C6h = PPts2(6,:);
    C6 = PPts2(5,:);
    T8h = PPts2(4,:);
    T8 = PPts2(3,:);
    T10h = PPts2(2,:);
    
    % Compute points from plane that passes M1(T9) - Nz - M2(T10)
    FCircBot = FindArcPoints([M1;Nz;M2],10,FacesT,VerticesN);
    FCircBot2 = FindArcPoints([M2;Nz;M1],10,FacesT,VerticesN);
    
    FTT9 = FCircBot(2,:);
    FT9 = FCircBot(3,:);
    FFT9 = FCircBot(4,:);
    F9 = FCircBot(5,:);
    AFF9 = FCircBot(6,:);
    AF9 = FCircBot(7,:);
    AFp9 = FCircBot(8,:);
    N1 = FCircBot(9,:);
    N1h = FCircBot(10,:);
%    Nz = FCircBot(11,:);
    N2h = FCircBot2(10,:);
    N2 = FCircBot2(9,:);
    AFp10 = FCircBot2(8,:);
    AF10 = FCircBot2(7,:);
    AFF10 = FCircBot2(6,:);
    F10 = FCircBot2(5,:);
    FFT10 = FCircBot2(4,:);
    FT10 = FCircBot2(3,:);
    FTT10 = FCircBot2(2,:);
    
    % Compute points from plane that passes M1(T9) - Iz - M2(T10)
    BCircBot = FindArcPoints([M1;Iz;M2],10,FacesT,VerticesN);
    BCircBot2 = FindArcPoints([M2;Iz;M1],10,FacesT,VerticesN);
    
    TTP9 = BCircBot(2,:);
    TP9 = BCircBot(3,:);
    TPP9 = BCircBot(4,:);
    P9 = BCircBot(5,:);
    PPO9 = BCircBot(6,:);
    PO9 = BCircBot(7,:);
    POO9 = BCircBot(8,:);
    I1 = BCircBot(9,:);
    I1h = BCircBot(10,:);
%    Iz = BCircBot(11,:);
    I2h = BCircBot2(10,:);
    I2 = BCircBot2(9,:);
    POO10 = BCircBot2(8,:);
    PO10 = BCircBot2(7,:);
    PPO10 = BCircBot2(6,:);
    P10 = BCircBot2(5,:);
    TPP10 = BCircBot2(4,:);
    TP10 = BCircBot2(3,:);
    TTP10 = BCircBot2(2,:);
    
    % Compute points from plane that passes T9h - NFpz - T10h
    FCircMid = FindArcPoints([T9h;NFpz;T10h],10,FacesT,VerticesN);
    FCircMid2 = FindArcPoints([T10h;NFpz;T9h],10,FacesT,VerticesN);
    
    FTT9h = FCircMid(2,:);
    FT9h = FCircMid(3,:);
    FFT9h = FCircMid(4,:);
    F9h = FCircMid(5,:);
    AFF9h = FCircMid(6,:);
    AF9h = FCircMid(7,:);
    AFp9h = FCircMid(8,:);
    NFp1 = FCircMid(9,:);
    NFp1h = FCircMid(10,:);
%    NFpz = FCircMid(11,:);
    NFp2h = FCircMid2(10,:);
    NFp2 = FCircMid2(9,:);
    AFp10h = FCircMid2(8,:);
    AF10h = FCircMid2(7,:);
    AFF10h = FCircMid2(6,:);
    F10h = FCircMid2(5,:);
    FFT10h = FCircMid2(4,:);
    FT10h = FCircMid2(3,:);
    FTT10h = FCircMid2(2,:);
    
    % Compute points from plane that passes T9h - OIz - T10h
    BCircMid = FindArcPoints([T9h;OIz;T10h],10,FacesT,VerticesN);
    BCircMid2 = FindArcPoints([T10h;OIz;T9h],10,FacesT,VerticesN);
    
    TTP9h = BCircMid(2,:);
    TP9h = BCircMid(3,:);
    TPP9h = BCircMid(4,:);
    P9h = BCircMid(5,:);
    PPO9h = BCircMid(6,:);
    PO9h = BCircMid(7,:);
    POO9h = BCircMid(8,:);
    OI1 = BCircMid(9,:);
    OI1h = BCircMid(10,:);
%    OIz = BCircMid(11,:);
    OI2h = BCircMid2(10,:);
    OI2 = BCircMid2(9,:);
    POO10h = BCircMid2(8,:);
    PO10h = BCircMid2(7,:);
    PPO10h = BCircMid2(6,:);
    P10h = BCircMid2(5,:);
    TPP10h = BCircMid2(4,:);
    TP10h = BCircMid2(3,:);
    TTP10h = BCircMid2(2,:);

    % Compute points from plane that passes T7 - Fpz - T8
    FCircTop = FindArcPoints([T7;Fpz;T8],10,FacesT,VerticesN);
    FCircTop2 = FindArcPoints([T8;Fpz;T7],10,FacesT,VerticesN);
    
    FTT7 = FCircTop(2,:);
    FT7 = FCircTop(3,:);
    FFT7 = FCircTop(4,:);
    F7 = FCircTop(5,:);
    AFF7 = FCircTop(6,:);
    AF7 = FCircTop(7,:);
    AFp7 = FCircTop(8,:);
    Fp1 = FCircTop(9,:);
    Fp1h = FCircTop(10,:);
%    Fpz = FCircTop(11,:);
    Fp2h = FCircTop2(10,:);
    Fp2 = FCircTop2(9,:);
    AFp8 = FCircTop2(8,:);
    AF8 = FCircTop2(7,:);
    AFF8 = FCircTop2(6,:);
    F8 = FCircTop2(5,:);
    FFT8 = FCircTop2(4,:);
    FT8 = FCircTop2(3,:);
    FTT8 = FCircTop2(2,:);
    
    % Compute points from plane that passes T7 - Oz - T8
    BCircTop = FindArcPoints([T7;Oz;T8],10,FacesT,VerticesN);
    BCircTop2 = FindArcPoints([T8;Oz;T7],10,FacesT,VerticesN);
    
    TTP7 = BCircTop(2,:);
    TP7 = BCircTop(3,:);
    TPP7 = BCircTop(4,:);
    P7 = BCircTop(5,:);
    PPO7 = BCircTop(6,:);
    PO7 = BCircTop(7,:);
    POO7 = BCircTop(8,:);
    O1 = BCircTop(9,:);
    O1h = BCircTop(10,:);
%    Oz = BCircTop(11,:);
    O2h = BCircTop2(10,:);
    O2 = BCircTop2(9,:);
    POO8 = BCircTop2(8,:);
    PO8 = BCircTop2(7,:);
    PPO8 = BCircTop2(6,:);
    P8 = BCircTop2(5,:);
    TPP8 = BCircTop2(4,:);
    TP8 = BCircTop2(3,:);
    TTP8 = BCircTop2(2,:);

    % Compute points from plane that passes FTT7 - FCCz - FTT8
    FCCArc = FindArcPoints([FTT7;FCCz;FTT8],12.5,FacesT,VerticesN);
    FCCArc2 = FindArcPoints([FTT8;FCCz;FTT7],12.5,FacesT,VerticesN);
    
    FTT7h = FCCArc(2,:);
    FCC5 = FCCArc(3,:);
    FCC5h = FCCArc(4,:);
    FCC3 = FCCArc(5,:);
    FCC3h = FCCArc(6,:);
    FCC1 = FCCArc(7,:);
    FCC1h = FCCArc(8,:);
%    FCCz = FCCArc(9,:);
    FCC2h = FCCArc2(8,:);
    FCC2 = FCCArc2(7,:);
    FCC4h = FCCArc2(6,:);
    FCC4 = FCCArc2(5,:);
    FCC6h = FCCArc2(4,:);
    FCC6 = FCCArc2(3,:);
    FTT8h = FCCArc2(2,:);  
    
    % Compute points from plane that passes FT7 - FCz - FT8
    FCArc = FindArcPoints([FT7;FCz;FT8],12.5,FacesT,VerticesN);
    FCArc2 = FindArcPoints([FT8;FCz;FT7],12.5,FacesT,VerticesN);
    
    FT7h = FCArc(2,:);
    FC5 = FCArc(3,:);
    FC5h = FCArc(4,:);
    FC3 = FCArc(5,:);
    FC3h = FCArc(6,:);
    FC1 = FCArc(7,:);
    FC1h = FCArc(8,:);
%    FCz = FCArc(9,:);
    FC2h = FCArc2(8,:);
    FC2 = FCArc2(7,:);
    FC4h = FCArc2(6,:);
    FC4 = FCArc2(5,:);
    FC6h = FCArc2(4,:);
    FC6 = FCArc2(3,:);
    FT8h = FCArc2(2,:);
    
    % Compute points from plane that passes FFT7 - FFCz - FFT8
    FFCArc = FindArcPoints([FFT7;FFCz;FFT8],12.5,FacesT,VerticesN);
    FFCArc2 = FindArcPoints([FFT8;FFCz;FFT7],12.5,FacesT,VerticesN);
    
    FFT7h = FFCArc(2,:);
    FFC5 = FFCArc(3,:);
    FFC5h = FFCArc(4,:);
    FFC3 = FFCArc(5,:);
    FFC3h = FFCArc(6,:);
    FFC1 = FFCArc(7,:);
    FFC1h = FFCArc(8,:);
%    FFCz = FFCArc(9,:);
    FFC2h = FFCArc2(8,:);
    FFC2 = FFCArc2(7,:);
    FFC4h = FFCArc2(6,:);
    FFC4 = FFCArc2(5,:);
    FFC6h = FFCArc2(4,:);
    FFC6 = FFCArc2(3,:);
    FFT8h = FFCArc2(2,:);    
    
    % Compute points from plane that passes F7 - Fz - F8
    FArc = FindArcPoints([F7;Fz;F8],12.5,FacesT,VerticesN);
    FArc2 = FindArcPoints([F8;Fz;F7],12.5,FacesT,VerticesN);
    
    F7h = FArc(2,:);
    F5 = FArc(3,:);
    F5h = FArc(4,:);
    F3 = FArc(5,:);
    F3h = FArc(6,:);
    F1 = FArc(7,:);
    F1h = FArc(8,:);
%    Fz = FArc(9,:);
    F2h = FArc2(8,:);
    F2 = FArc2(7,:);
    F4h = FArc2(6,:);
    F4 = FArc2(5,:);
    F6h = FArc2(4,:);
    F6 = FArc2(3,:);
    F8h = FArc2(2,:);    
    
    % Compute points from plane that passes AFF7 - AFFz - AFF8
    AFFArc = FindArcPoints([AFF7;AFFz;AFF8],12.5,FacesT,VerticesN);
    AFFArc2 = FindArcPoints([AFF8;AFFz;AFF7],12.5,FacesT,VerticesN);
    
    AFF7h = AFFArc(2,:);
    AFF5 = AFFArc(3,:);
    AFF5h = AFFArc(4,:);
    AFF3 = AFFArc(5,:);
    AFF3h = AFFArc(6,:);
    AFF1 = AFFArc(7,:);
    AFF1h = AFFArc(8,:);
%    AFFz = AFFArc(9,:);
    AFF2h = AFFArc2(8,:);
    AFF2 = AFFArc2(7,:);
    AFF4h = AFFArc2(6,:);
    AFF4 = AFFArc2(5,:);
    AFF6h = AFFArc2(4,:);
    AFF6 = AFFArc2(3,:);
    AFF8h = AFFArc2(2,:);      
    
    % Compute points from plane that passes AF7 - AFz - AF8
    AFArc = FindArcPoints([AF7;AFz;AF8],12.5,FacesT,VerticesN);
    AFArc2 = FindArcPoints([AF8;AFz;AF7],12.5,FacesT,VerticesN);
    
    AF7h = AFArc(2,:);
    AF5 = AFArc(3,:);
    AF5h = AFArc(4,:);
    AF3 = AFArc(5,:);
    AF3h = AFArc(6,:);
    AF1 = AFArc(7,:);
    AF1h = AFArc(8,:);
%    AFz = AFArc(9,:);
    AF2h = AFArc2(8,:);
    AF2 = AFArc2(7,:);
    AF4h = AFArc2(6,:);
    AF4 = AFArc2(5,:);
    AF6h = AFArc2(4,:);
    AF6 = AFArc2(3,:);
    AF8h = AFArc2(2,:);     
    
    % Compute points from plane that passes AFp7 - AFpz - AFp8
    AFpArc = FindArcPoints([AFp7;AFpz;AFp8],25,FacesT,VerticesN);
    AFpArc2 = FindArcPoints([AFp8;AFpz;AFp7],25,FacesT,VerticesN);
    
    AFp5 = AFpArc(2,:);
    AFp3 = AFpArc(3,:);
    AFp1 = AFpArc(4,:);
%    AFpz = AFpArc(5,:);
    AFp2 = AFpArc2(4,:);
    AFp4 = AFpArc2(3,:);
    AFp6 = AFpArc2(2,:);
    
    % Compute points from plane that passes TTP7 - CCPz - TTP8
    CCPArc = FindArcPoints([TTP7;CCPz;TTP8],12.5,FacesT,VerticesN);
    CCPArc2 = FindArcPoints([TTP8;CCPz;TTP7],12.5,FacesT,VerticesN);
    
    TTP7h = CCPArc(2,:);
    CCP5 = CCPArc(3,:);
    CCP5h = CCPArc(4,:);
    CCP3 = CCPArc(5,:);
    CCP3h = CCPArc(6,:);
    CCP1 = CCPArc(7,:);
    CCP1h = CCPArc(8,:);
%    CCPz = CCPArc(9,:);
    CCP2h = CCPArc2(8,:);
    CCP2 = CCPArc2(7,:);
    CCP4h = CCPArc2(6,:);
    CCP4 = CCPArc2(5,:);
    CCP6h = CCPArc2(4,:);
    CCP6 = CCPArc2(3,:);
    TTP8h = CCPArc2(2,:);  

    % Compute points from plane that passes TP7 - CPz - TP8
    CPArc = FindArcPoints([TP7;CPz;TP8],12.5,FacesT,VerticesN);
    CPArc2 = FindArcPoints([TP8;CPz;TP7],12.5,FacesT,VerticesN);
    
    TP7h = CPArc(2,:);
    CP5 = CPArc(3,:);
    CP5h = CPArc(4,:);
    CP3 = CPArc(5,:);
    CP3h = CPArc(6,:);
    CP1 = CPArc(7,:);
    CP1h = CPArc(8,:);
%    CPz = CPArc(9,:);
    CP2h = CPArc2(8,:);
    CP2 = CPArc2(7,:);
    CP4h = CPArc2(6,:);
    CP4 = CPArc2(5,:);
    CP6h = CPArc2(4,:);
    CP6 = CPArc2(3,:);
    TP8h = CPArc2(2,:);  
    
    % Compute points from plane that passes TPP7 - CPPz - TPP8
    CPPArc = FindArcPoints([TPP7;CPPz;TPP8],12.5,FacesT,VerticesN);
    CPPArc2 = FindArcPoints([TPP8;CPPz;TPP7],12.5,FacesT,VerticesN);
    
    TPP7h = CPPArc(2,:);
    CPP5 = CPPArc(3,:);
    CPP5h = CPPArc(4,:);
    CPP3 = CPPArc(5,:);
    CPP3h = CPPArc(6,:);
    CPP1 = CPPArc(7,:);
    CPP1h = CPPArc(8,:);
%    CPPz = CPPArc(9,:);
    CPP2h = CPPArc2(8,:);
    CPP2 = CPPArc2(7,:);
    CPP4h = CPPArc2(6,:);
    CPP4 = CPPArc2(5,:);
    CPP6h = CPPArc2(4,:);
    CPP6 = CPPArc2(3,:);
    TPP8h = CPPArc2(2,:);     
    
    % Compute points from plane that passes P7 - Pz - P8
    PArc = FindArcPoints([P7;Pz;P8],12.5,FacesT,VerticesN);
    PArc2 = FindArcPoints([P8;Pz;P7],12.5,FacesT,VerticesN);
    
    P7h = PArc(2,:);
    P5 = PArc(3,:);
    P5h = PArc(4,:);
    P3 = PArc(5,:);
    P3h = PArc(6,:);
    P1 = PArc(7,:);
    P1h = PArc(8,:);
%    Pz = PArc(9,:);
    P2h = PArc2(8,:);
    P2 = PArc2(7,:);
    P4h = PArc2(6,:);
    P4 = PArc2(5,:);
    P6h = PArc2(4,:);
    P6 = PArc2(3,:);
    P8h = PArc2(2,:);     
    
    % Compute points from plane that passes PPO7 - PPOz - PPO8
    PPOArc = FindArcPoints([PPO7;PPOz;PPO8],12.5,FacesT,VerticesN);
    PPOArc2 = FindArcPoints([PPO8;PPOz;PPO7],12.5,FacesT,VerticesN);
    
    PPO7h = PPOArc(2,:);
    PPO5 = PPOArc(3,:);
    PPO5h = PPOArc(4,:);
    PPO3 = PPOArc(5,:);
    PPO3h = PPOArc(6,:);
    PPO1 = PPOArc(7,:);
    PPO1h = PPOArc(8,:);
%    PPOz = PPOArc(9,:);
    PPO2h = PPOArc2(8,:);
    PPO2 = PPOArc2(7,:);
    PPO4h = PPOArc2(6,:);
    PPO4 = PPOArc2(5,:);
    PPO6h = PPOArc2(4,:);
    PPO6 = PPOArc2(3,:);
    PPO8h = PPOArc2(2,:);
        
    % Compute points from plane that passes PO7 - POz - PO8
    POArc = FindArcPoints([PO7;POz;PO8],12.5,FacesT,VerticesN);
    POArc2 = FindArcPoints([PO8;POz;PO7],12.5,FacesT,VerticesN);
    
    PO7h = POArc(2,:);
    PO5 = POArc(3,:);
    PO5h = POArc(4,:);
    PO3 = POArc(5,:);
    PO3h = POArc(6,:);
    PO1 = POArc(7,:);
    PO1h = POArc(8,:);
%    POz = POArc(9,:);
    PO2h = POArc2(8,:);
    PO2 = POArc2(7,:);
    PO4h = POArc2(6,:);
    PO4 = POArc2(5,:);
    PO6h = POArc2(4,:);
    PO6 = POArc2(3,:);
    PO8h = POArc2(2,:);    
    
    % Compute points from plane that passes POO7 - POOz - POO8
    POOArc = FindArcPoints([POO7;POOz;POO8],25,FacesT,VerticesN);
    POOArc2 = FindArcPoints([POO8;POOz;POO7],25,FacesT,VerticesN);
    
    POO5 = POOArc(2,:);
    POO3 = POOArc(3,:);
    POO1 = POOArc(4,:);
%    POOz = POOArc(5,:);
    POO2 = POOArc2(4,:);
    POO4 = POOArc2(3,:);
    POO6 = POOArc2(2,:);    
    
    
    
    if System ==1
        
        % EEG 10-5 position labels
        EEGLab = {'Nz';'N1h';'NFp1h';'NFpz';'NFp2h';'N2h';'N1';'NFp1';'Fp1';
                  'Fp1h';'Fpz';'Fp2h';'Fp2';'NFp2';'N2';'AFp9';'AFp9h';'AFp7';
                  'AFp5';'AFp3';'AFp1';'AFpz';'AFp2';'AFp4';'AFp6';'AFp8';
                  'AFp10h';'AFp10';'AF9';'AF9h';'AF7';'AF7h';'AF5';'AF5h';
                  'AF3';'AF3h';'AF1';'AF1h';'AFz';'AF2h';'AF2';'AF4h';'AF4';
                  'AF6h';'AF6';'AF8h';'AF8';'AF10h';'AF10';'AFF9';'AFF9h';
                  'AFF7';'AFF7h';'AFF5';'AFF5h';'AFF3';'AFF3h';'AFF1';'AFF1h';
                  'AFFz';'AFF2h';'AFF2';'AFF4h';'AFF4';'AFF6h';'AFF6';'AFF8h';
                  'AFF8';'AFF10h';'AFF10';'F9';'F9h';'F7';'F7h';'F5';'F5h';
                  'F3';'F3h';'F1';'F1h';'Fz';'F2h';'F2';'F4h';'F4';'F6h';
                  'F6';'F8h';'F8';'F10h';'F10';'FFT9';'FFT9h';'FFT7';'FFT7h';
                  'FFC5';'FFC5h';'FFC3';'FFC3h';'FFC1';'FFC1h';'FFCz';'FFC2h';
                  'FFC2';'FFC4h';'FFC4';'FFC6h';'FFC6';'FFT8h';'FFT8';'FFT10h';
                  'FFT10';'FT9';'FT9h';'FT7';'FT7h';'FC5';'FC5h';'FC3';'FC3h';
                  'FC1';'FC1h';'FCz';'FC2h';'FC2';'FC4h';'FC4';'FC6h';'FC6';
                  'FT8h';'FT8';'FT10h';'FT10';'FTT9';'FTT9h';'FTT7';'FTT7h';
                  'FCC5';'FCC5h';'FCC3';'FCC3h';'FCC1';'FCC1h';'FCCz';'FCC2h';
                  'FCC2';'FCC4h';'FCC4';'FCC6h';'FCC6';'FTT8h';'FTT8';'FTT10h';
                  'FTT10';'T9(M1)';'T9h';'T7';'T7h';'C5';'C5h';'C3';'C3h';
                  'C1';'C1h';'Cz';'C2h';'C2';'C4h';'C4';'C6h';'C6';'T8h';
                  'T8';'T10h';'T10(M2)';'TTP9';'TTP9h';'TTP7';'TTP7h';'CCP5';
                  'CCP5h';'CCP3';'CCP3h';'CCP1';'CCP1h';'CCPz';'CCP2h';'CCP2';
                  'CCP4h';'CCP4';'CCP6h';'CCP6';'TTP8h';'TTP8';'TTP10h';
                  'TTP10';'TP9';'TP9h';'TP7';'TP7h';'CP5';'CP5h';'CP3';'CP3h';
                  'CP1';'CP1h';'CPz';'CP2h';'CP2';'CP4h';'CP4';'CP6h';'CP6';
                  'TP8h';'TP8';'TP10h';'TP10';'TPP9';'TPP9h';'TPP7';'TPP7h';
                  'CPP5';'CPP5h';'CPP3';'CPP3h';'CPP1';'CPP1h';'CPPz';'CPP2h';
                  'CPP2';'CPP4h';'CPP4';'CPP6h';'CPP6';'TPP8h';'TPP8';'TPP10h';
                  'TPP10';'P9';'P9h';'P7';'P7h';'P5';'P5h';'P3';'P3h';'P1';
                  'P1h';'Pz';'P2h';'P2';'P4h';'P4';'P6h';'P6';'P8h';'P8';
                  'P10h';'P10';'PPO9';'PPO9h';'PPO7';'PPO7h';'PPO5';'PPO5h';
                  'PPO3';'PPO3h';'PPO1';'PPO1h';'PPOz';'PPO2h';'PPO2';'PPO4h';
                  'PPO4';'PPO6h';'PPO6';'PPO8h';'PPO8';'PPO10h';'PPO10';'PO9';
                  'PO9h';'PO7';'PO7h';'PO5';'PO5h';'PO3';'PO3h';'PO1';'PO1h';
                  'POz';'PO2h';'PO2';'PO4h';'PO4';'PO6h';'PO6';'PO8h';'PO8';
                  'PO10h';'PO10';'POO9';'POO9h';'POO7';'POO5';'POO3';'POO1';
                  'POOz';'POO2';'POO4';'POO6';'POO8';'POO10h';'POO10';'I1';
                  'OI1';'O1';'O1h';'Oz';'O2h';'O2';'OI2';'I2';'I1h';'OI1h';
                  'OIz';'OI2h';'I2h';'Iz'};

        EEGPts = [Nz;N1h;NFp1h;NFpz;NFp2h;N2h;N1;NFp1;Fp1;Fp1h;Fpz;Fp2h;Fp2;
                  NFp2;N2;AFp9;AFp9h;AFp7;AFp5;AFp3;AFp1;AFpz;AFp2;AFp4;AFp6;
                  AFp8;AFp10h;AFp10;AF9;AF9h;AF7;AF7h;AF5;AF5h;AF3;AF3h;AF1;
                  AF1h;AFz;AF2h;AF2;AF4h;AF4;AF6h;AF6;AF8h;AF8;AF10h;AF10;
                  AFF9;AFF9h;AFF7;AFF7h;AFF5;AFF5h;AFF3;AFF3h;AFF1;AFF1h;
                  AFFz;AFF2h;AFF2;AFF4h;AFF4;AFF6h;AFF6;AFF8h;AFF8;AFF10h;
                  AFF10;F9;F9h;F7;F7h;F5;F5h;F3;F3h;F1;F1h;Fz;F2h;F2;F4h;F4;
                  F6h;F6;F8h;F8;F10h;F10;FFT9;FFT9h;FFT7;FFT7h;FFC5;FFC5h;
                  FFC3;FFC3h;FFC1;FFC1h;FFCz;FFC2h;FFC2;FFC4h;FFC4;FFC6h;
                  FFC6;FFT8h;FFT8;FFT10h;FFT10;FT9;FT9h;FT7;FT7h;FC5;FC5h;
                  FC3;FC3h;FC1;FC1h;FCz;FC2h;FC2;FC4h;FC4;FC6h;FC6;FT8h;FT8;
                  FT10h;FT10;FTT9;FTT9h;FTT7;FTT7h;FCC5;FCC5h;FCC3;FCC3h;
                  FCC1;FCC1h;FCCz;FCC2h;FCC2;FCC4h;FCC4;FCC6h;FCC6;FTT8h;
                  FTT8;FTT10h;FTT10;M1;T9h;T7;T7h;C5;C5h;C3;C3h;C1;C1h;Cz;
                  C2h;C2;C4h;C4;C6h;C6;T8h;T8;T10h;M2;TTP9;TTP9h;TTP7;TTP7h;
                  CCP5;CCP5h;CCP3;CCP3h;CCP1;CCP1h;CCPz;CCP2h;CCP2;CCP4h;
                  CCP4;CCP6h;CCP6;TTP8h;TTP8;TTP10h;TTP10;TP9;TP9h;TP7;TP7h;
                  CP5;CP5h;CP3;CP3h;CP1;CP1h;CPz;CP2h;CP2;CP4h;CP4;CP6h;CP6;
                  TP8h;TP8;TP10h;TP10;TPP9;TPP9h;TPP7;TPP7h;CPP5;CPP5h;CPP3;
                  CPP3h;CPP1;CPP1h;CPPz;CPP2h;CPP2;CPP4h;CPP4;CPP6h;CPP6;
                  TPP8h;TPP8;TPP10h;TPP10;P9;P9h;P7;P7h;P5;P5h;P3;P3h;P1;P1h;
                  Pz;P2h;P2;P4h;P4;P6h;P6;P8h;P8;P10h;P10;PPO9;PPO9h;PPO7;
                  PPO7h;PPO5;PPO5h;PPO3;PPO3h;PPO1;PPO1h;PPOz;PPO2h;PPO2;
                  PPO4h;PPO4;PPO6h;PPO6;PPO8h;PPO8;PPO10h;PPO10;PO9;PO9h;PO7;
                  PO7h;PO5;PO5h;PO3;PO3h;PO1;PO1h;POz;PO2h;PO2;PO4h;PO4;PO6h;
                  PO6;PO8h;PO8;PO10h;PO10;POO9;POO9h;POO7;POO5;POO3;POO1;
                  POOz;POO2;POO4;POO6;POO8;POO10h;POO10;I1;OI1;O1;O1h;Oz;
                  O2h;O2;OI2;I2;I1h;OI1h;OIz;OI2h;I2h;Iz];                    
        
    elseif System == 2
        
        % EEG 10-10 position labels
        EEGLab = {'Nz';'N1';'Fp1';'Fpz';'Fp2';'N2';'AF9';'AF7';'AF3';'AFz';
            	  'AF4';'AF8';'AF10';'F9';'F7';'F5';'F3';'F1';'Fz';'F2';'F4';
                  'F6';'F8';'F10';'FT9';'FT7';'FC5';'FC3';'FC1';'FCz';'FC2';
                  'FC4';'FC6';'FT8';'FT10';'T9(M1)';'T7';'C5';'C3';'C1';'Cz';
                  'C2';'C4';'C6';'T8';'T10(M2)';'TP9';'TP7';'CP5';'CP3';'CP1';
                  'CPz';'CP2';'CP4';'CP6';'TP8';'TP10';'P9';'P7';'P5';'P3';
                  'P1';'Pz';'P2';'P4';'P6';'P8';'P10';'PO9';'PO7';'PO3';'POz';
                  'PO4';'PO8';'PO10';'I1';'O1';'Oz';'O2';'I2';'Iz'};

        EEGPts = [Nz;N1;Fp1;Fpz;Fp2;N2;AF9;AF7;AF3;AFz;AF4;AF8;AF10;F9;F7;
                  F5;F3;F1;Fz;F2;F4;F6;F8;F10;FT9;FT7;FC5;FC3;FC1;FCz;FC2;
                  FC4;FC6;FT8;FT10;M1;T7;C5;C3;C1;Cz;C2;C4;C6;T8;M2;TP9;TP7;
                  CP5;CP3;CP1;CPz;CP2;CP4;CP6;TP8;TP10;P9;P7;P5;P3;P1;Pz;P2;
                  P4;P6;P8;P10;PO9;PO7;PO3;POz;PO4;PO8;PO10;I1;O1;Oz;O2;I2;
                  Iz];
              
    elseif System == 3
        % EEG 10-20 position labels
        EEGLab = {'Fp1'; 'Fp2'; 'F7'; 'F3'; 'Fz'; 'F4'; 'F8'; 'T7'; 'C3';
               'Cz'; 'C4'; 'T8'; 'P7'; 'P3'; 'Pz'; 'P4'; 'P8'; 'O1'; 'O2'};

        EEGPts = [Fp1;Fp2;F7;F3;Fz;F4;F8;T7;C3;Cz;C4;T8;P7;P3;Pz;P4;P8;O1;O2];
        
    elseif System == 4
        % EEG 10-10 position labels and NIRS modified 10-10 position labels
        
        % EEG 10-10 position labels
        EEGLab = {'Nz';'Fp1';'Fpz';'Fp2';'AF7';'AF3';'AFz';'AF4';'AF8';'F7';
                  'F5';'F3';'F1';'Fz';'F2';'F4';'F6';'F8';'FT7';'FC5';'FC3';
                  'FC1';'FCz';'FC2';'FC4';'FC6';'FT8';'T9(M1)';'T7';'C5';
                  'C3';'C1';'Cz';'C2';'C4';'C6';'T8';'T10(M2)';'TP7';'CP5';
                  'CP3';'CP1';'CPz';'CP2';'CP4';'CP6';'TP8';'P7';'P5';'P3';
                  'P1';'Pz';'P2';'P4';'P6';'P8';'PO7';'PO3';'POz';'PO4';
                  'PO8';'O1';'Oz';'O2';'Iz'};

        EEGPts = [Nz;Fp1;Fpz;Fp2;AF7;AF3;AFz;AF4;AF8;F7;
                  F5;F3;F1;Fz;F2;F4;F6;F8;FT7;FC5;FC3;FC1;FCz;FC2;
                  FC4;FC6;FT8;M1;T7;C5;C3;C1;Cz;C2;C4;C6;T8;M2;TP7;
                  CP5;CP3;CP1;CPz;CP2;CP4;CP6;TP8;P7;P5;P3;P1;Pz;P2;
                  P4;P6;P8;PO7;PO3;POz;PO4;PO8;O1;Oz;O2;
                  Iz];
              
        NIRSPts = [NFp1h;NFp2h;AFp9h;AFp1;AFp2;AFp10h;AFF9h;AFF5;AFF1;AFF2;
                  AFF6;AFF10h;FFT9h;FFT7h;FFC5h;FFC3h;FFC1h;FFC2h;FFC4h;
                  FFC6h;FFT8h;FFT10h;FTT9h;FTT7h;FCC5h;FCC3h;FCC1h;FCC2h;
                  FCC4h;FCC6h;FTT8h;FTT10h;TTP9h;TTP7h;CCP5h;CCP3h;CCP1h;
                  CCP2h;CCP4h;CCP6h;TTP8h;TTP10h;TPP9h;TPP7h;CPP5h;CPP3h;
                  CPP1h;CPP2h;CPP4h;CPP6h;TPP8h;TPP10h;PPO9h;PPO5;PPO1;PPO2;
                  PPO6;PPO10h;POO9h;POO1;POO2;POO10h;OI1h;OI2h];
              
        NIRSLab = {'Fpd1';'Fpd2';'Fps3';'Fps1';'Fps2';'Fps4';'Fd9';'Fd7';
                   'Fd3';'Fd4';'Fd8';'Fd10';'Fs9';'Fs7';'Fs5';'Fs3';'Fs1';
                   'Fs2';'Fs4';'Fs6';'Fs8';'Fs10';'Td9';'Td7';'Cd5';'Cd3';
                   'Cd1';'Cd2';'Cd4';'Cd6';'Td8';'Td10';'Ts9';'Ts7';'Cs5';
                   'Cs3';'Cs1';'Cs2';'Cs4';'Cs6';'Ts8';'Ts10';'Pd9';'Pd7';
                   'Pd5';'Pd3';'Pd1';'Pd2';'Pd4';'Pd6';'Pd8';'Pd10';'Ps9';
                   'Ps7';'Ps3';'Ps4';'Ps8';'Ps10';'Od3';'Od1';'Od2';'Od4';
                   'Os1';'Os2'};
               
        EEGPts = [EEGPts;NIRSPts];
        EEGLab = [EEGLab;NIRSLab];
        
    else
        error('Matlab error :: System input case unavailable, try (1),(2),(3), or (4)')
    end
    
    
    % Un-Normalize Vertices from range of [0-100] to return to original
    % units.
    VerticesT = (VerticesN*(Vmax+Vmin)/(100))-Vmin;
    EEGPts = (EEGPts*(Vmax+Vmin)/(100))-Vmin;
    
end


%% Paolo Giacometti
% Function to optimize the Cz position in a surface mesh. First, the
% positions are calculated for the inion-nasion and the
% preauricular arcs from the desired Cz. Then, the error in the positioning
% of Cz as the absolute difference between the desired position and the 
% calculated position. It optimizes the position of Cz so that it is at the
% center of both arcs.
% Inputs: Cz - [x,y,z] coordinates of initial guess of center point.
%         Fiducials - [x,y,z] coordinates of fiducial points for both arcs.
%         MeshFaces - [:,:,:] indices of all the vertices in each
%         surface mesh face.
%         MeshVertices - [x,y,z] coordinates of all the vertices in a
%         surface mesh.
% Output: OptCz - [x,y,z] coordinates of optimized center point.
%         error - [20x1] optimization error for each iteration alternating
%         from inion nasion to preauricular arc error. 

function [OptCz,error] = OptimizeCz(Cz, Fiducials, MeshFaces, MeshVertices)
    
    I = Fiducials(1,:);
    N = Fiducials(2,:);
    P1 = Fiducials(3,:);
    P2 = Fiducials(4,:);

    Percent = 2;
    
    maxiter = 10;
    error = zeros(2*maxiter,1);
    CzVec = zeros(2*maxiter,3);
    for i=1:maxiter

        % Set fiducials for arc
        PointsIN = [I; Cz; N]; % Inion-Nasion arc

        % Calculate points in inion-nasion arc at percentages
        PtPercIN = FindArcPoints(PointsIN, Percent, MeshFaces, MeshVertices);

        % Get Cz from calculated points (will change index (26) if
        % percentage is changed)
        CzIN = PtPercIN(26,:);
        
        % Compute Error
        error(2*i-1,1) = sqrt(sum((CzIN-Cz).^2));
        
        % Generate new Cz as the mean of previous two
        newCz = mean([CzIN;Cz]);
        
        % Add new Cz to list
        CzVec(2*i-1,:) = newCz;

        % Set fiducials for arc
        PointsP = [P1; newCz; P2]; % Preauricular arc.    

        % Calculate points in preauricular arc at percentages
        PtPercP = FindArcPoints(PointsP, Percent, MeshFaces, MeshVertices);
        
        % Get Cz from calculated points (will change index (26) if
        % percentage is changed)
        CzP = PtPercP(26,:);
        
        % Calculate error from distance of new Cz to previous Cz
        error(2*i,1) = sqrt(sum((CzP-newCz).^2));
        
        % Generate new Cz as the mean of previous two
        newCz2 = mean([CzP;newCz]);
        
        % Add new Cz to list
        CzVec(2*i,:) = newCz2;
                
        % Make new Cz current Cz
        Cz = newCz2;
        
    end
    
    % Find minimum error
    [minerr minerridx] = min(error);
    
    % Return Cz with minimum error
    OptCz = CzVec(minerridx,:);

end


%% Paolo Giacometti
% Function inputs a set of 3 points that define a plane. It uses these
% points to transform a whole set of coordinate points such that the origin
% lies in the midpoint between p2 and p3 and p1 lies normal to the segment
% p2-p3 in the negative x-direction.
% Inputs: xyz [x,y,z] coordinates of points to be transformed
%         p1 - [x,y,z] coordinates of point that defines a plane
%         p2 - [x,y,z] coordinates of point that defines a plane
%         p3 - [x,y,z] coordinates of point that defines a plane
% Output: XYZT [x,y,z] transformed coordinates

function [XYZT]=TransformCoords(xyz,p1,p2,p3)

    % Compute plane defined by p1,p2, and p3
    % Compute constants for plane equation [ax+by+cz+d=0]
    a = (p2(1,2)-p1(1,2))*(p3(1,3)-p1(1,3))-(p3(1,2)-p1(1,2))*(p2(1,3)-p1(1,3));
    b = (p2(1,3)-p1(1,3))*(p3(1,1)-p1(1,1))-(p3(1,3)-p1(1,3))*(p2(1,1)-p1(1,1));
    c = (p2(1,1)-p1(1,1))*(p3(1,2)-p1(1,2))-(p3(1,1)-p1(1,1))*(p2(1,2)-p1(1,2));
    d = -(a*p1(1,1)+b*p1(1,2)+c*p1(1,3)); % This would work with p1, p2, or p3

    % Compute normalized normal to plane
    np = [a b c]/sqrt(a^2+b^2+c^2);    
    
    % Compute fourth point to define transformation
    p4 = mean([p2;p3])+np;
    
    % Compute distance between p2 and p3
    d1 = abs(sqrt(sum(diff([p2;p3]).^2)));
    
    % Compute distance between p1 and new origin
    d2 = abs(sqrt(sum(diff([p1;mean([p2;p3])]).^2)));
    
    % Compute matrix Pb for transformation Pa = T*Pb
    Pb = [[p1' p2' p3' p4']; [1 1 1 1]];
        
    % Compute matrix Pa for transformation Pa = T*Pb
    Pa = [[0 -d1/2 d1/2 0];[d2 0 0 0];[0 0 0 1];[1 1 1 1]];

    % Compute transformation matrix for Pa = T*Pb 
    T=Pa*Pb'*(Pb*Pb')^(-1);
    
    % Compute matrix for xyz points for global transformation
    XYZ = [xyz'; ones(1,size(xyz,1))];
    
    % Compute transformation of points
    XYZT = T*XYZ;
    
    % Transpose XYZT and ignore 4th row to obtain transformed coordinates
    XYZT = XYZT(1:3,:)';
    
end



