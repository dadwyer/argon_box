#!/usr/bin/env python
#
# Geant4-based simulation of particle interactions in a box of liquid argon
#
#  For use in quick studies of particle interactions and energy
#  deposition in liquid argon.  Uses standard G4 physics lists, and
#  includes a primitive option to generate a 3D map of the ionization
#  energy deposition from the particles.
#
#  Includes two particle generator options:
#   - Monoenergetic particles of a specific type (e.g. e-, n, p, mu+)
#   - HepEVT-formatted input data
#
#  Output is a plain ROOT tree containing information on each Geant4
#  track step.  If the primitive energy deposition simulation is
#  enabled, then a list of energy depositions and their coordinates
#  are also included.
# 
# dadwyer@lbl.gov, Aug. 18, 2017

# Import necessary Geant4 classes
# Generators
from Geant4 import G4VUserPrimaryGeneratorAction, G4ParticleGun 
from Geant4 import G4ParticleTable, G4PrimaryParticle, G4PrimaryVertex
# User-defined actions
from Geant4 import G4UserRunAction, G4UserEventAction, G4UserSteppingAction
# Physics processes
from Geant4 import FTFP_BERT, QGSP_BERT, QGSP_BERT_HP
# Run Manager and UI commands
from Geant4 import gRunManager, gApplyUICommand, G4UImanager
# Physical Units
from Geant4 import MeV, cm, mm, g, mole, cm3, deg, GeV, m, c_light, ns
# Materials
from Geant4 import G4Material, G4NistManager, G4Element, G4Material
# Geometry
from Geant4 import G4VUserDetectorConstruction, G4LogicalVolume, G4PVPlacement
from Geant4 import G4Box, G4Tubs, G4SubtractionSolid, G4UnionSolid
from Geant4 import G4ThreeVector, G4RotationMatrix
# Visualization
from Geant4 import gVisManager
# Random numbers
from Geant4 import G4RandomDirection, HepRandom
# G4 Configuration
from Geant4 import G4UserLimits
# Pre-defined NIST materials
import g4py.NISTmaterials

# General tools
from copy import copy
from time import sleep

# General math
from math import sqrt, cos, sin, acos

# ROOT
from array import array
from ROOT import TTree, TFile

# Command-line arguments
from argparse import ArgumentParser

import numpy as np

class MySimulation:
    "My Simple Simulation"

    def __init__(self):
        '''Constructor'''
        self._ofilename = 'argon_sim.root'
        self._random_seed = 23
        self._ofile = None
        self._otree = None
#        self._gtree = None
        self._treebuffer = None
        self._materials = {}
        self._geom = None
        self._physlist_name = 'QGSP_BERT'
        self._physics_list = None
        self._source = None
        self._shift = None
        self._detX = None
        self._detY = None
        self._detZ = None
        self._energies = None
        self._source_position = None
        self._generator = None
        self._include_edepsim = False
        self._edep_step = 2.0 * mm
        return

    def initialize(self):
        '''Initialize the simulation'''
        # Parse command-line arguments
        self._parse_args()
        # Prepare output file
        self._prepare_output()
        # Initialize random number generator
        self._init_random()
        # Prepare materials
        self._init_materials()
        # Prepare geometry
        self._init_geometry()
        # Prepare physics list
        self._init_physics_list()
        # Prepare generator
        self._init_generator()
        # Prepare user actions
        self._init_user_actions()
        # Prepare run manager
        self._init_run_manager()
        return

    def run(self):
        '''Run the simulation'''
        if self._nevents > 0:
          print "Provided events ", self._nevents
          sleep(2)
          gRunManager.BeamOn(self._nevents)
        else:
          nFile = TFile(self._source,'READ')
          t = nFile.Get("gRooTracker")
          print "Got entries from tree ", t.GetEntries()
          sleep(2)
          gRunManager.BeamOn(t.GetEntries())
        return

    def finalize(self):
        '''Finalize the simulation'''
        self._close_output()
        return

    ######################

    def _parse_args(self):
        '''Parse command line arguments'''
        parser = ArgumentParser()
        parser.add_argument('--nevents', type=int, help='Number of events',
                            default=1)
        parser.add_argument('--source', help='Particle name or generator',
                            default='e-')
        parser.add_argument('--energy', type=float, help='Particle energy',
                            default=2.0 * GeV)
        parser.add_argument('--output', help='Output filename')
        parser.add_argument('--seed', type=int, help='Random seed')
        parser.add_argument('--enable_edepsim',
                            help='Enable energy deposition in simulation',
                            action='store_true')
        parser.add_argument('--edep_step', type=float,
                            help='Maximum step size in energy dep simulation',
                            default=2.0 * mm)
        parser.add_argument('--physlist',
                            help='G4 Physics List: FTFP_BERT, QGSP_BERT, QGSP_BERT_HP',
                            default='QGSP_BERT')
        parser.add_argument('--shift', type=float,
                            help='Shift from GENIE geom -> centered argon box in cm',
                            default='0.0')
        parser.add_argument('--detX', type=float,
                            help='1/2 full Detector X dimension (m)',
                            default='0.0')
        parser.add_argument('--detY', type=float,
                            help='1/2 full Detector Y dimension (m)',
                            default='0.0')
        parser.add_argument('--detZ', type=float,
                            help='1/2 full Detector Z dimension (m)',
                             default='0.0')

        self._args = parser.parse_args()
        if self._args.nevents is not None:
            self._nevents = self._args.nevents
        if self._args.source is not None:
            self._source = self._args.source
        if self._args.energy is not None:
            self._energies = [self._args.energy * GeV,]
        if self._args.output is not None:
            self._ofilename = self._args.output
        if self._args.seed is not None:
            self._random_seed = self._args.seed
        if self._args.shift is not None:
            self._shift = self._args.shift
        if self._args.detX is not None:
            self._detX = self._args.detX
        if self._args.detY is not None:
            self._detY = self._args.detY
        if self._args.detZ is not None:
            self._detZ = self._args.detZ
        if self._args.enable_edepsim:
            self._include_edepsim = True
        if self._args.edep_step:
            self._edep_step = self._args.edep_step
        print "Configuration:"
        print "  nevents:",self._nevents
        print "   source:",self._source
        print "   shift:",self._shift
        print "   detX:",self._detX
        print "   detY:",self._detY
        print "   detZ:",self._detZ
        print "   energies:",self._energies
        print "   output:",self._ofilename
        print " physlist:",self._physlist_name
        print "     seed:",self._random_seed
        print "  edepsim:",self._include_edepsim
        print " edepstep:",self._edep_step
        return

    def _prepare_output(self):
        '''Prepare ROOT tree for output data'''
        ofile = TFile(self._ofilename,'RECREATE')

        otree = TTree('argon','Argon Simulation')
        tb = TreeBuffer()
        tb.nevents = self._nevents
        tb.maxInit = 100
        tb.maxTracks = 1000000
        tb.maxNQ = 10000000
        tb.shift = self._shift
        tb.ev = array('i',[0])
        # Ancestor particles (e.g. neutrino)
        #   Note: For simplicity, assume only one potential ancestor
        tb.pida = array('i',[0])
        tb.xa = array('d',[0])
        tb.ya = array('d',[0])
        tb.za = array('d',[0])
        tb.ta = array('d',[0])
        tb.pxa = array('d',[0])
        tb.pya = array('d',[0])
        tb.pza = array('d',[0])
        tb.ekina = array('d',[0])
        tb.ma = array('d',[0])        
        tb.code = bytearray(51) 
        tb.codemap = {}
        # Geant4 initial state particles
        tb.ni = array('i',[0])
        tb.pidi = array('i',[0]*tb.maxInit)
        tb.xi = array('d',[0]*tb.maxInit)
        tb.yi = array('d',[0]*tb.maxInit)
        tb.zi = array('d',[0]*tb.maxInit)
        tb.ti = array('d',[0]*tb.maxInit)
        tb.pxi = array('d',[0]*tb.maxInit)
        tb.pyi = array('d',[0]*tb.maxInit)
        tb.pzi = array('d',[0]*tb.maxInit)
        tb.ekini = array('d',[0]*tb.maxInit)
        tb.mi = array('d',[0]*tb.maxInit)
        # Geant4 track step data
        tb.nstep = array('i',[0])
        tb.tid = array('i',[0]*tb.maxTracks)
        tb.pid = array('i',[0]*tb.maxTracks)
        tb.parid = array('i',[0]*tb.maxTracks)
        tb.xs = array('d',[0]*tb.maxTracks)
        tb.ys = array('d',[0]*tb.maxTracks)
        tb.zs = array('d',[0]*tb.maxTracks)
        tb.xe = array('d',[0]*tb.maxTracks)
        tb.ye = array('d',[0]*tb.maxTracks)
        tb.ze = array('d',[0]*tb.maxTracks)
        tb.te = array('d',[0]*tb.maxTracks)
        
        tb.pxs = array('d',[0]*tb.maxTracks)
        tb.pys = array('d',[0]*tb.maxTracks)
        tb.pzs = array('d',[0]*tb.maxTracks)
        tb.pxe = array('d',[0]*tb.maxTracks)
        tb.pye = array('d',[0]*tb.maxTracks)
        tb.pze = array('d',[0]*tb.maxTracks)
        
        tb.ekin = array('d',[0]*tb.maxTracks)
        tb.edep = array('d',[0]*tb.maxTracks)
        # Geant4 energy deposition data
        tb.nq = array('i',[0])
        tb.tidq = array('i',[0]*tb.maxNQ)
        tb.pidq = array('i',[0]*tb.maxNQ)
        tb.sidq = array('i',[0]*tb.maxNQ)
        tb.dq = array('d',[0]*tb.maxNQ)
        tb.xq = array('d',[0]*tb.maxNQ)
        tb.yq = array('d',[0]*tb.maxNQ)
        tb.zq = array('d',[0]*tb.maxNQ)
        # Hook for ancestor particle
        tb.ancestor = None
       
        otree.Branch('ev',tb.ev,'ev/I')
        otree.Branch('pida',tb.pida,'pida/I')
        otree.Branch('xa',tb.xa,'xa/D')
        otree.Branch('ya',tb.ya,'ya/D')
        otree.Branch('za',tb.za,'za/D')
        otree.Branch('ta',tb.ta,'ta/D')
        otree.Branch('pxa',tb.pxa,'pxa/D')
        otree.Branch('pya',tb.pya,'pya/D')
        otree.Branch('pza',tb.pza,'pza/D')
        otree.Branch('ekina',tb.ekina,'ekina/D')
        otree.Branch('ma',tb.ma,'ma/D')
        otree.Branch('code',tb.code,'code[51]/C')
        otree.Branch('ni',tb.ni,'ni/I')
        otree.Branch('pidi',tb.pidi,'pidi[ni]/I')
        otree.Branch('xi',tb.xi,'xi[ni]/D')
        otree.Branch('yi',tb.yi,'yi[ni]/D')
        otree.Branch('zi',tb.zi,'zi[ni]/D')
        otree.Branch('ti',tb.ti,'ti[ni]/D')
        otree.Branch('pxi',tb.pxi,'pxi[ni]/D')
        otree.Branch('pyi',tb.pyi,'pyi[ni]/D')
        otree.Branch('pzi',tb.pzi,'pzi[ni]/D')
        otree.Branch('ekini',tb.ekini,'ekini[ni]/D')
        otree.Branch('mi',tb.mi,'mi[ni]/D')
        otree.Branch('nstep',tb.nstep,'nstep/I')
        otree.Branch('tid',tb.tid,'tid[nstep]/I')
        otree.Branch('pid',tb.pid,'pid[nstep]/I')
        otree.Branch('parid',tb.parid,'parid[nstep]/I')
        otree.Branch('ekin',tb.ekin,'ekin[nstep]/D')
        otree.Branch('edep',tb.edep,'edep[nstep]/D')
        otree.Branch('xs',tb.xs,'xs[nstep]/D')
        otree.Branch('ys',tb.ys,'ys[nstep]/D')
        otree.Branch('zs',tb.zs,'zs[nstep]/D')
        otree.Branch('xe',tb.xe,'xe[nstep]/D')
        otree.Branch('ye',tb.ye,'ye[nstep]/D')
        otree.Branch('ze',tb.ze,'ze[nstep]/D')
        otree.Branch('te',tb.te,'te[nstep]/D')
        otree.Branch('pxs',tb.pxs,'pxs[nstep]/D')
        otree.Branch('pys',tb.pys,'pys[nstep]/D')
        otree.Branch('pzs',tb.pzs,'pzs[nstep]/D')
        otree.Branch('pxe',tb.pxe,'pxe[nstep]/D')
        otree.Branch('pye',tb.pye,'pye[nstep]/D')
        otree.Branch('pze',tb.pze,'pze[nstep]/D')
        otree.Branch('nq',tb.nq,'nq/I')
        otree.Branch('tidq',tb.tidq,'tidq[nq]/I')
        otree.Branch('pidq',tb.pidq,'pidq[nq]/I')
        otree.Branch('sidq',tb.sidq,'sidq[nq]/I')
        otree.Branch('dq',tb.dq,'dq[nq]/D')
        otree.Branch('xq',tb.xq,'xq[nq]/D')
        otree.Branch('yq',tb.yq,'yq[nq]/D')
        otree.Branch('zq',tb.zq,'zq[nq]/D')

        self._ofile = ofile
        self._otree = otree
        self._treebuffer = tb

        return
    
    def _close_output(self):
        '''Write ROOT data, and close file'''
        self._ofile.cd()
        self._otree.Write()
        self._ofile.Close()
        return

    def _init_random(self):
        '''Initialize random number generator'''
        HepRandom.setTheSeed(self._random_seed)
        return

    def _init_materials(self):
        '''Prepare list of materials'''
        g4py.NISTmaterials.Construct()
        # Define liquid argon (name,z,a,density)
        LAr_mat = G4Material("liquidArgon", 18., 39.95*g/mole, 1.390*g/cm3)
        self._materials['liquidArgon'] = LAr_mat
        return

    def _init_geometry(self):
        '''Initialize geant geometry'''
        self._geom = MyDetectorConstruction(materials=self._materials,detX=self._detX,detY=self._detY,detZ=self._detZ)
        gRunManager.SetUserInitialization(self._geom)        
        return

    def _init_physics_list(self):
        '''Initialize the physics list'''
        # Use standard physics list
        if self._physlist_name == 'FTFP_BERT':
            self._physics_list = FTFP_BERT()
        elif self._physlist_name == 'QGSP_BERT':
            self._physics_list = QGSP_BERT()
        elif self._physlist_name == 'QGSP_BERT_HP':
            self._physics_list = QGSP_BERT_HP()
        else:
            raise ValueError("Invalid physics list: \'%r\'.  Please choose from:FTFP_BERT, QGSP_BERT, QGSP_BERT_HP" % self._physlist_name)
        #self._physics_list.RegisterPhysics(G4MyStepLimiterPhysics())
        gRunManager.SetUserInitialization(self._physics_list)
        return

    def _init_generator(self):
        '''Initialize particle generator'''
        self._generator = MyGenieEvtGeneratorAction(
            genieEvtFilename=self._source,
            treebuffer=self._treebuffer,
            ofile=self._ofile)
        gRunManager.SetUserAction(self._generator)
        return

    def _init_user_actions(self):
        '''Initialize user actions'''
        self._run_action = MyRunAction()
        self._event_action = MyEventAction(treebuffer=self._treebuffer,
                                           outputtree=self._otree)
        self._step_action = MySteppingAction(treebuffer=self._treebuffer,
                                             geometry=self._geom,
                                             materials=self._materials,
                                             edepsim=self._include_edepsim,
                                             edepstep=self._edep_step)
        gRunManager.SetUserAction(self._run_action)
        gRunManager.SetUserAction(self._event_action)
        gRunManager.SetUserAction(self._step_action)
        return
            
    def _init_run_manager(self):
        '''Initialize the Geant4 run manager'''
        #gApplyUICommand('/run/setCut 0.1 mm')
        gApplyUICommand('/run/initialize')
        return
        
###################################################################

# Tree Buffer
class TreeBuffer:
    '''Dummy class for collecting tree data buffers'''
    pass

# Tree Buffer
class GTreeBuffer:
    '''Dummy class for collecting tree data buffers'''
    pass

# Detector Construction
class MyDetectorConstruction(G4VUserDetectorConstruction):
    "My Detector Construction"

    def __init__(self, materials={},detX=100.,detY=100.,detZ=100.):
        G4VUserDetectorConstruction.__init__(self)
        self.vols = []
        self.rots = []
        self.materials = materials
        # Make table of indices for each volume region 
        self.volumeidx = {
            'World':0,
        }
        self.detX = detX
        self.detY = detY
        self.detZ = detZ
        print self.detX, self.detY, self.detZ
        
    def Construct(self):
        '''Construct geometry'''
        # World (box of liquid argon)
        #world_s = G4Box("World", 19.5*m, 1.5*m, 2.5*m)
        world_s = G4Box("World", self.detX*m, self.detY*m, self.detZ*m)
        world_l = G4LogicalVolume(world_s, self.materials['liquidArgon'],
                                  "World")
        world_p = G4PVPlacement(None,                  #no rotation
                                G4ThreeVector(),       #at (0,0,0)
                                world_l,               #its logical volume
                                "World",               #its name
                                None,                  #its mother  volume
                                False,                 #no boolean operation
                                0)                     #copy number
        self.vols += [world_s, world_l, world_p]

        # Restrict step size to ensure proper charge voxelization?
        maxStep = 3 * mm
        self._userLimits = G4UserLimits(maxStep)
        world_l.SetUserLimits(self._userLimits)
        
        # Return world volume
        return world_p


###################################################################
# GENIE port 
class MyGenieEvtGeneratorAction(G4VUserPrimaryGeneratorAction):
    "Generator based on RooTracker input data"
    
    def __init__(self, genieEvtFilename, treebuffer, ofile):
        G4VUserPrimaryGeneratorAction.__init__(self)
        self.isInitialized = False
        self.inputFile = genieEvtFilename
        self.hepEvts = None
        self.currentEventIdx = 0
        self._tb = treebuffer
        self._ofile = ofile
        return

    def initialize(self):
        # Prepare generator
        ifile = TFile(self.inputFile, 'READ')
        datatree = ifile.Get('gRooTracker')
        self._ofile.cd()
        print 'GENIE Generator: Loaded Tree with %d events from %s' % (
            datatree.GetEntries(), self.inputFile)
        self.hepEvts = self.parseGenieEvts(datatree) 
        self.isInitialized = True
        return

    def parseGenieEvts(self,datatree):
        '''Parse GENIE input data'''
        hepEvts = []
        entry_count = 0
        while entry_count < self._tb.nevents:
            datatree.GetEntry(entry_count)
            currentEvent = {}
            currentEvent['event_num'] = datatree.EvtNum
            n_particles = datatree.StdHepN
            currentEvent['particles'] = []

            #event position in meters
            evt_pos = [datatree.EvtVtx[0],
                       datatree.EvtVtx[1],
                       datatree.EvtVtx[2]]

            #retuns lists: info[particle]
            status = datatree.StdHepStatus
            pdg = datatree.StdHepPdg
            first_parent = datatree.StdHepFm
            first_child = datatree.StdHepFd
            last_parent = datatree.StdHepLm
            last_child = datatree.StdHepLd

            #funky nested list stuff
            p4 = datatree.StdHepP4
            p4_shaped = np.reshape(p4, (len(p4)/4,4) )

            q4 = datatree.StdHepX4
            q4_shaped = np.reshape(q4, (len(p4)/4,4) )
            
            tempcode = bytearray(51)
            tempcode[:51] = ""
#            self._tb.code[:51] = "" 
            event_code = str(datatree.EvtCode).split(';')
            counter = 0
            for part in event_code:
              if 'proc' in part:
                this_code = part + ';' + event_code[counter+1]
#                self._tb.code[:len(this_code)+1] = part + ';' + event_code[counter + 1]
                tempcode[:len(this_code)+1] = part + ';' + event_code[counter + 1]
                self._tb.codemap[datatree.EvtNum] = tempcode
#                print self._tb.code 
              counter = counter + 1

            for npart in range(datatree.StdHepN):
                if ( status[npart] > 1 ):
                    n_particles = n_particles - 1 
                    continue
                if ( status[npart] == 0 ):
                    if ( abs(pdg[npart]) == 12 or abs(pdg[npart]) == 14 ):
                        print 'Found intial state neutrino'
                    else:      
                        n_particles = n_particles - 1
                        print 'initial state', abs(pdg[npart])
                        continue
                px = p4_shaped[npart][0]
                py = p4_shaped[npart][1]
                pz = p4_shaped[npart][2]
                E = p4_shaped[npart][3]
                p2 = px**2 + py**2 + pz**2              
                     
                x = (evt_pos[0]*10**3 + q4_shaped[npart][0]*(10**-9) + self._tb.shift*10.)*mm #shifting over for events
                y = (evt_pos[1]*10**3 + q4_shaped[npart][1]*(10**-9))*mm
                z = (evt_pos[2]*10**3 + q4_shaped[npart][2]*(10**-9))*mm
                t = q4_shaped[npart][3]*(10**-9)*mm/c_light 


                if (E*E - p2) < 0:
                    m = 0
                else:
                    m = sqrt(E*E - p2)

                particle = {
                    'index':npart, #is this right?
                    'status':status[npart],
                    'pdgid':pdg[npart],
                    'first_parent':first_parent[npart],
                    'first_child':first_child[npart],
                    'last_parent':last_parent[npart],
                    'last_child':last_child[npart],
                    'momentum':[px*GeV, py*GeV, pz*GeV],
                    'energy':E*GeV,
                    'mass':m*GeV,
                    'position':[x, y, z],
                    'time':t,
                }
                currentEvent['n_particles'] = n_particles
                currentEvent['particles'].append(particle)
            hepEvts.append( copy(currentEvent) )
            entry_count = entry_count + 1 

#        for entry in datatree:
#            currentEvent = {}
#            self._tb.code = ""
#            currentEvent['event_num'] = entry.EvtNum
#            n_particles = entry.StdHepN
#            currentEvent['particles'] = []
#
#            #event position in meters
#            evt_pos = [entry.EvtVtx[0],
#                       entry.EvtVtx[1],
#                       entry.EvtVtx[2]]
#
#            #retuns lists: info[particle]
#            status = entry.StdHepStatus
#            pdg = entry.StdHepPdg
#            first_parent = entry.StdHepFm
#            first_child = entry.StdHepFd
#            last_parent = entry.StdHepLm
#            last_child = entry.StdHepLd
#
#            #funky nested list stuff
#            p4 = entry.StdHepP4
#            p4_shaped = np.reshape(p4, (len(p4)/4,4) )
#
#            q4 = entry.StdHepX4
#            q4_shaped = np.reshape(q4, (len(p4)/4,4) )
#
#            event_code = str(entry.EvtCode).split(';')
#            counter = 0
#            for part in event_code:
#              if 'proc' in part:
#                self._tb.code = part + ';' + event_code[counter + 1]
#                print part + ';' + event_code[counter + 1]
#                print self._tb.code
#              counter = counter + 1
#
#            for npart in range(entry.StdHepN):
##                print npart, status[npart], pdg[npart]
#                if ( status[npart] > 1 ):
##                    print 'Other state', pdg[npart]
##                    if ( (abs(pdg[npart]) != 14)  or (abs(pdg[npart]) != 12) ):                    
##                       print 'not inital neutrino or final state. skipping'
#                    n_particles = n_particles - 1 
#                    continue
#                if ( status[npart] == 0 ):
#                    if ( abs(pdg[npart]) == 12 or abs(pdg[npart]) == 14 ):
#                        print 'Found intial state neutrino'
#                    else:      
#                        n_particles = n_particles - 1
#                        print 'initial state', abs(pdg[npart])
#                        continue
#                px = p4_shaped[npart][0]
#                py = p4_shaped[npart][1]
#                pz = p4_shaped[npart][2]
#                E = p4_shaped[npart][3]
#                p2 = px**2 + py**2 + pz**2              
#                     
#                x = (evt_pos[0]*10**3 + q4_shaped[npart][0]*(10**-9) + self._tb.shift*10.)*mm #shifting over for events
#                y = (evt_pos[1]*10**3 + q4_shaped[npart][1]*(10**-9))*mm
#                z = (evt_pos[2]*10**3 + q4_shaped[npart][2]*(10**-9))*mm
#                t = q4_shaped[npart][3]*(10**-9)*mm/c_light 
#
#
#                if (E*E - p2) < 0:
#                    m = 0
#                else:
#                    m = sqrt(E*E - p2)
#
#                particle = {
#                    'index':npart, #is this right?
#                    'status':status[npart],
#                    'pdgid':pdg[npart],
#                    'first_parent':first_parent[npart],
#                    'first_child':first_child[npart],
#                    'last_parent':last_parent[npart],
#                    'last_child':last_child[npart],
#                    'momentum':[px*GeV, py*GeV, pz*GeV],
#                    'energy':E*GeV,
#                    'mass':m*GeV,
#                    'position':[x, y, z],
#                    'time':t,
#                }
#                currentEvent['n_particles'] = n_particles
#                currentEvent['particles'].append(particle)
#            hepEvts.append( copy(currentEvent) )
        return hepEvts


    def GeneratePrimaries(self, event):
        if not self.isInitialized:
            self.initialize()
        if self.currentEventIdx > len(self.hepEvts):
            print ('Error: No more HepEVT data! (max events=%d)'
                   % self.currentEventIdx)
            return event
        currentEvt = self.hepEvts[self.currentEventIdx]
        self._tb.ancestor = None
        # Set default position and time to zero 
        time = 0
        finalStateParticles = []
        for heppart in currentEvt['particles']:
            status = heppart['status']
            if status not in [0,1]:
                # ignore all but initial, final state particles
                #  FIXME: should probably throw a warning
                continue
            heppos = heppart['position']
            hepmom = heppart['momentum']
            if status == 0:
                # Initial state particle.  Add info to event tree
                self._tb.ancestor = heppart
            elif status == 1:
                # Final state particle.  Add to event
                position = G4ThreeVector(heppos[0],heppos[1],heppos[2])
                particle = G4PrimaryParticle(heppart['pdgid'])
                particle.Set4Momentum(hepmom[0],hepmom[1],hepmom[2],
                                      heppart['energy'])
                # Ensure mass is exact
                particle.SetMass(heppart['mass'])
                finalStateParticles.append([particle, position, time])
        for [particle, position, time] in finalStateParticles:
            vertex = G4PrimaryVertex(position, time)
            vertex.SetPrimary(particle)
            event.AddPrimaryVertex(vertex)
        self.currentEventIdx += 1

# ------------------------------------------------------------------
class MyRunAction(G4UserRunAction):
    "My Run Action"

    def EndOfRunAction(self, run):
        print "*** End of Run"
        print "- Run sammary : (id= %d, #events= %d)" \
            % (run.GetRunID(), run.GetNumberOfEventToBeProcessed())

# ------------------------------------------------------------------
class MyEventAction(G4UserEventAction):
    "My Event Action"
    def __init__(self, treebuffer, outputtree):
        '''Constructor'''
        G4UserEventAction.__init__(self)
        self._tb = treebuffer
        self._otree = outputtree
    
    def BeginOfEventAction(self, event):
        self._tb.ev[0] = -1
        self._tb.pida[0] = 0
        self._tb.xa[0] = 0
        self._tb.ya[0] = 0
        self._tb.za[0] = 0
        self._tb.ta[0] = 0
        self._tb.pxa[0] = 0
        self._tb.pya[0] = 0
        self._tb.pza[0] = 0
        self._tb.ekina[0] = 0
        self._tb.ma[0] = 0
        self._tb.edep[0] = 0
        self._tb.ni[0] = 0
        self._tb.nstep[0] = 0
        self._tb.nq[0] = 0
        return
    
    def EndOfEventAction(self, event):
        '''Record event'''
        tb = self._tb
        # Event ID
        tb.ev[0] = event.GetEventID()
        tb.code[:51] = ""
        tb.code[:len(tb.codemap[event.GetEventID()])+1] = tb.codemap[event.GetEventID()]        
        print tb.code
        print len(tb.code)
        # Ancestor particle (e.g. neutrino)
        if tb.ancestor != None:
            # Log info
            heppart = tb.ancestor
            tb.pida[0] = heppart['pdgid']
            heppos = heppart['position']
            hepmom = heppart['momentum']
            tb.xa[0] = (heppos[0] - tb.shift*10.)/ cm
            tb.ya[0] = heppos[1] / cm
            tb.za[0] = heppos[2] / cm
            tb.ta[0] = heppart['time'] / ns
            tb.pxa[0] = hepmom[0] / GeV
            tb.pya[0] = hepmom[1] / GeV
            tb.pza[0] = hepmom[2] / GeV
            tb.ekina[0] = (heppart['energy'] - heppart['mass']) / GeV
            tb.ma[0] = heppart['mass'] / GeV
        # Primary particles
        ni = 0
        for pv_idx in range(event.GetNumberOfPrimaryVertex()):
            pvtx = event.GetPrimaryVertex(pv_idx)
            for pp_idx in range(pvtx.GetNumberOfParticle()):
                ppart = pvtx.GetPrimary(pp_idx)
                tb.pidi[ni] = ppart.GetPDGcode()
                tb.xi[ni] = (pvtx.GetX0() - tb.shift*10.) / cm
                tb.yi[ni] = pvtx.GetY0() / cm
                tb.zi[ni] = pvtx.GetZ0() / cm
                tb.ti[ni] = pvtx.GetT0() / ns
                tb.pxi[ni] = ppart.GetPx() / GeV
                tb.pyi[ni] = ppart.GetPy() / GeV
                tb.pzi[ni] = ppart.GetPz() / GeV
                # Ick, G4py interface doesn't provide particle energy func :/
                pmomSq = ppart.GetPx()**2 + ppart.GetPy()**2 + ppart.GetPz()**2
                pmass = ppart.GetMass()
                penergy = sqrt(pmomSq + pmass**2)
                tb.ekini[ni] = (penergy - pmass) / GeV
                tb.mi[ni] = pmass / GeV
                ni += 1
                if ni == tb.maxInit:
                    print "Warning: Reached max number of initial particles."
                    break
        tb.ni[0] = ni
        self._otree.Fill()
        return

# ------------------------------------------------------------------
class MySteppingAction(G4UserSteppingAction):
    "My Stepping Action"
    def __init__(self, treebuffer, geometry, materials, edepsim, edepstep):
        '''Constructor'''
        G4UserSteppingAction.__init__(self)
        self._tb = treebuffer
        self._geom = geometry
        self._materials = materials
        self._include_edepsim = edepsim
        self._edep_step = edepstep
    
    def UserSteppingAction(self, step):
        '''Collect data for current simulation step'''
        tb = self._tb
        istp = tb.nstep[0]
        if istp>=tb.maxTracks:
            print 'Reached maximum tracks:',istp
            return
        track = step.GetTrack()
        prestep = step.GetPreStepPoint()
        poststep = step.GetPostStepPoint()
        tb.tid[istp] = track.GetTrackID()
        tb.pid[istp] = track.GetDefinition().GetPDGEncoding()
        tb.parid[istp] = track.GetParentID()
        tb.ekin[istp] = prestep.GetKineticEnergy() / GeV
        # Sum simulated energy deposition by volume
        tb.edep[istp] = step.GetTotalEnergyDeposit() / GeV
        # Capture step position
        prepos = prestep.GetPosition()
        postpos = poststep.GetPosition()
        tb.xs[istp] = (prepos.x - tb.shift*10.) / cm
        tb.ys[istp] = prepos.y / cm
        tb.zs[istp] = prepos.z / cm 
        tb.xe[istp] = (postpos.x - tb.shift*10.) / cm
        tb.ye[istp] = postpos.y / cm
        tb.ze[istp] = postpos.z / cm
        tb.te[istp] = poststep.GetGlobalTime() / ns
        # Add in momentum 
        premom = prestep.GetMomentum()
        postmom = poststep.GetMomentum()
        tb.pxs[istp] = premom.x / GeV
        tb.pys[istp] = premom.y / GeV
        tb.pzs[istp] = premom.z / GeV 
        tb.pxe[istp] = postmom.x / GeV
        tb.pye[istp] = postmom.y / GeV
        tb.pze[istp] = postmom.z / GeV

        tb.nstep[0] += 1
        if not self._include_edepsim:
            # Stop here, unless simple energy deposition simulation is enabled
            return
        # Arrrgghhh! no python hook to non-ionizing energy deposit...
        #eIon = step.GetTotalEnergyDeposit()-step.GetNonIonizingEnergyDeposit()
        eIon = step.GetTotalEnergyDeposit()
        if track.GetDefinition().GetPDGCharge() != 0.0 and eIon > 0:
            deltaPos = postpos - prepos
            iq = tb.nq[0]
            dnq = int(deltaPos.mag() / self._edep_step)           
            drq = G4ThreeVector(deltaPos)
            drq.setMag(self._edep_step)
            for qidx in range(dnq):
                curpos = prepos + qidx * drq
                tb.tidq[iq+qidx] = track.GetTrackID()
                tb.pidq[iq+qidx] = track.GetDefinition().GetPDGEncoding()
                tb.sidq[iq+qidx] = istp
                tb.dq[iq+qidx] = (eIon / dnq) / MeV
                tb.xq[iq+qidx] = (curpos.x - tb.shift*10.) / cm
                tb.yq[iq+qidx] = curpos.y / cm
                tb.zq[iq+qidx] = curpos.z / cm
            tb.nq[0] += dnq
        return

    
###################################################################

if '__main__' == __name__:
    # Run the simulation
    print "Running simulation"
    mysim = MySimulation()
    mysim.initialize()
    mysim.run()
    mysim.finalize()
    
