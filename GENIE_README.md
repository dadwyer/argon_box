Slightly modified version of Argon Box script. Takes in gRooTracker files, places final state particles in a box of liquid Argon


Usage:
$python argon_box_genie.py --nevents=1 --source=gRooTrackerInput.root --output=out.root --enable_edepsim --detX=25 --detY=5 --detZ=5 --shift=1800

Extra/Changed Options:
nevents -> set to 0 for all events, else n events

det[X,Y,Z] -> half of block dimension (m)

shift -> G4 requires the positions to be relative to the center of the detector, so shift takes the 'real world' X position, 
         translates to detector coordinates, and then back into original coordinates. Could work also put in Y,Z shifts..
         Units are cm -- sorry for inconsistency!!



Added an output of the gRooTracker event code. To grab from the output file:

$ root out.root
root: char code[51] 
root: argon->SetBranchAddress("code",&code)
root: argon->GetEntry(0)
root: print code


