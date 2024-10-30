# OnScalpPrintedEEG
Authors: Eric Li, Sandhya Tiku, Luize Scalco de Vasconcelos - Nanshu Lu Research Group - University of Texas at Austin

Code to generate a custom EEG electrode and interconnect layout on a 3D STL head model and corresponding G-code for digital printing using 5-axis robot (3 translation+2 rotation axes) and ink dispensing (ON/OFF)

To generate the custom electrode and interconnect layout:

1. Start Matlab
2. Open the file "GenerateCustomLayout.m"
3. Run

To generate the printer control instructions:

1. Start Matlab
2. Open the file "PostProcessor/PrinterControl_GcodeGenerator.m"
3. Run


The Mesh2EEG functions ("ComputeEEGPos.m", "SortEdgeNodes.m", "FindArcPoints.m", "MeshPlaneIntersectPoints.m") for automatic calculation of the International 10-20, 10-10, and 10-5 scalp coordinates of EEG electrodes on a boundary element mesh of a human head are by Giacometti, P., Perdue, K.L., Diamond, S.G. from "Algorithm to find high density EEG scalp coordinates and analysis of their correspondence to structural and functional regions of the brain. J Neurosci Methods. 229:84-96. (2014)"
