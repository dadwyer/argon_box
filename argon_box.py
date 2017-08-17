#!/usr/bin/env python
#
# Geant4-based simulation of interactions in box of argon
#


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
from Geant4 import MeV, cm, mm, g, mole, cm3, deg, GeV, m, c_light
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

# General math
from math import sqrt, cos, sin, acos

# ROOT
from array import array
from ROOT import TTree, TFile

# Command-line arguments
from argparse import ArgumentParser


class MySimulation:
    "My Simple Simulation"

    def __init__(self):
        '''Constructor'''
        self._ofilename = 'argon_sim.root'
        self._random_seed = 23
        self._ofile = None
        self._otree = None
        self._treebuffer = None
        self._materials = {}
        self._geom = None
        self._physics_list = None
        self._source = None
        self._energies = None
        self._source_position = None
        self._generator = None
        self._include_edepsim = True
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
        gRunManager.BeamOn(self._nevents)
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
                            default=2.0)
        parser.add_argument('--output', help='Output filename')
        parser.add_argument('--seed', type=int, help='Random seed')
        parser.add_argument('--disable_edepsim',
                            help='Disable energy deposition in simulation',
                            action='store_true')
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
        if self._args.disable_edepsim:
            self._include_edepsim = False
        print "Configuration:"
        print "  nevents:",self._nevents
        print "   source:",self._source
        print "   energies:",self._energies
        print "   output:",self._ofilename
        print "     seed:",self._random_seed
        print "  edepsim:",self._include_edepsim
        return

    def _prepare_output(self):
        '''Prepare ROOT tree for output data'''
        ofile = TFile(self._ofilename,'RECREATE')
        otree = TTree('argon','Argon Simulation')
        tb = TreeBuffer()

        tb.maxTracks = 100000
        tb.maxNQ = 10000000
        tb.ev = array('i',[0])
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
        tb.ekin = array('d',[0]*tb.maxTracks)
        tb.edep = array('d',[0]*tb.maxTracks)
        tb.nq = array('i',[0])
        tb.tidq = array('i',[0]*tb.maxNQ)
        tb.pidq = array('i',[0]*tb.maxNQ)
        tb.sidq = array('i',[0]*tb.maxNQ)
        tb.dq = array('d',[0]*tb.maxNQ)
        tb.xq = array('d',[0]*tb.maxNQ)
        tb.yq = array('d',[0]*tb.maxNQ)
        tb.zq = array('d',[0]*tb.maxNQ)
        
        otree.Branch('ev',tb.ev,'ev/I')
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
        self._geom = MyDetectorConstruction(materials=self._materials)
        gRunManager.SetUserInitialization(self._geom)        
        return

    def _init_physics_list(self):
        '''Initialize the physics list'''
        # Use standard physics list
        self._physics_list = QGSP_BERT()
        #self._physics_list = QGSP_BERT_HP()
        #self._physics_list = FTFP_BERT()
        #self._physics_list.RegisterPhysics(G4MyStepLimiterPhysics())
        gRunManager.SetUserInitialization(self._physics_list)
        return

    def _init_generator(self):
        '''Initialize particle generator'''
        if self._source.endswith('hepevt'):
            self._generator = MyHepEvtGeneratorAction(
                hepEvtFilename=self._source)
        else:
            self._source_pos = G4ThreeVector(0,0,0)
            self._generator = MyParticleGeneratorAction(
                particleName=self._source,
                energies=self._energies,
                position=self._source_pos)
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
                                             edepsim=self._include_edepsim)
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

# Detector Construction
class MyDetectorConstruction(G4VUserDetectorConstruction):
    "My Detector Construction"

    def __init__(self, materials={}):
        G4VUserDetectorConstruction.__init__(self)
        self.vols = []
        self.rots = []
        self.materials = materials
        # Make table of indices for each volume region 
        self.volumeidx = {
            'World':0,
        }
        
    def Construct(self):
        '''Construct geometry'''
        # World (box of liquid argon)
        world_s = G4Box("World", 100*m, 100*m, 100*m)
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

# Particle interaction generator
class MyParticleGeneratorAction(G4VUserPrimaryGeneratorAction):
    "Generator for single type of particles (e.g. e-, gammas, etc)"

    def __init__(self, particleName='e-', energies=[1.0*GeV,],
                 position=G4ThreeVector(0,0,0)):
        G4VUserPrimaryGeneratorAction.__init__(self)
        #self.particleGun= G4ParticleGun(1)
        self.isInitialized = False
        self.particleName = particleName
        self.energies = energies
        self.position = position
        self.particleDef = None
        pass

    def initialize(self):
        # Prepare generator
        particleTable = G4ParticleTable.GetParticleTable()
        self.particleDef = particleTable.FindParticle(self.particleName)
        self.isInitialized = True
        return
    
    def GeneratePrimaries(self, event):
        #self.particleGun.GeneratePrimaryVertex(event)
        if not self.isInitialized:
            self.initialize()
        position = self.GenerateVertexPosition()
        time = 0.
        vertex = G4PrimaryVertex(position, time)
        mass = self.particleDef.GetPDGMass()
        for enr in self.energies:
            particle = G4PrimaryParticle(self.particleDef.GetPDGEncoding())
            particle.Set4Momentum(sqrt(enr**2-mass**2),0,0,enr)
            # Ensure mass is exact
            particle.SetMass(mass)
            # Set direction
            particleDirection = self.GenerateParticleDirection()
            particle.SetMomentumDirection(particleDirection)
            vertex.SetPrimary(particle)
        event.AddPrimaryVertex(vertex)

    def GenerateVertexPosition(self):
        '''Generate a vertex position'''
        return G4ThreeVector(self.position)

    def GenerateParticleDirection(self):
        '''Generate a particle direction'''
        return G4ThreeVector(0,0,1)

# Particle interaction generator
class MyHepEvtGeneratorAction(G4VUserPrimaryGeneratorAction):
    "Generator based on HepEVT input data"

    def __init__(self, hepEvtFilename):
        G4VUserPrimaryGeneratorAction.__init__(self)
        #self.particleGun= G4ParticleGun(1)
        self.isInitialized = False
        self.inputFile = hepEvtFilename
        self.hepEvts = None
        self.currentEventIdx = 0
        pass

    def initialize(self):
        # Prepare generator
        datalines = open(self.inputFile).readlines()
        self.hepEvts = self.parseHepEvts(datalines) 
        print 'HepEVT Generator: Loaded %d events from %s' % (
            len(self.hepEvts), self.inputFile)
        self.isInitialized = True
        return

    def parseHepEvts(self,data):
        '''Parse HepMC HEPEVT input data'''
        hepEvts = []
        in_block = False
        idx = 0
        while idx < len(data):
            line = data[idx]
            line = line.strip()
            print 'line:',line
            if line.startswith('___'):
                # Ignore extra formatting lines
                idx += 1
                continue
            if line.startswith('***** HEPEVT'):
                # Found a new HEPEVT block.  Process.
                currentEvent = {}
                tokens = line.split()
                event_num = int((tokens[4].strip())[:-1])
                n_particles = int(tokens[5])
                currentEvent['event_num'] = event_num
                currentEvent['n_particles'] = n_particles
                currentEvent['particles'] = []
                for npart in range(n_particles + 1):
                    # Parse info from each particle
                    first_line = data[idx + 5 + npart*3]
                    second_line = data[idx + 5 + npart*3 + 1]
                    first_line = first_line.replace('(','').replace(')','').replace(',','')
                    second_line = second_line.replace('(','').replace(')','').replace(',','')
                    #print '  first:',first_line
                    #print '  second:',second_line
                    first_toks = first_line.strip().split()
                    second_toks = second_line.strip().split()
                    particle = {
                        'index':int(first_toks[0]),
                        'status':int(first_toks[1]),
                        'pdgid':int(second_toks[0]),
                        'first_parent':int(first_toks[2]),
                        'first_child':int(first_toks[3]),
                        'last_parent':int(second_toks[1]),
                        'last_child':int(second_toks[2]),
                        'momentum':[float(first_toks[4])*GeV,
                                    float(first_toks[5])*GeV,
                                    float(first_toks[6])*GeV],
                        'energy':float(first_toks[7])*GeV,
                        'mass':float(first_toks[8])*GeV,
                        'position':[float(second_toks[3])*mm,
                                    float(second_toks[4])*mm,
                                    float(second_toks[5])*mm],
                        'time':float(second_toks[6])*mm/c_light,
                    }
                    currentEvent['particles'].append(particle)
                hepEvts.append( copy(currentEvent) )
                idx += (5 + (n_particles+1)*3 + 1)
                print ' HepEVT Loaded:', len(hepEvts)
                #if len(hepEvts) > 20: break
        return hepEvts
        
    def GeneratePrimaries(self, event):
        #self.particleGun.GeneratePrimaryVertex(event)
        if not self.isInitialized:
            self.initialize()
        if self.currentEventIdx > len(self.hepEvts):
            print ('Error: No more HepEVT data! (max events=%d)'
                   % self.currentEventIdx)
            return event
        currentEvt = self.hepEvts[self.currentEventIdx]
        # Set default position and time to zero 
        position = G4ThreeVector(0,0,0)
        time = 0
        finalStateParticles = []
        print " Curevent: ",self.currentEventIdx
        for heppart in currentEvt['particles']:
            status = heppart['status']
            if status == 0:
                # initial state particle.  get position and time
                heppos = heppart['position']
                position = G4ThreeVector(heppos[0],heppos[1],heppos[2])
                time = heppart['time']
            elif status == 1:
                # final state particle.  add to event
                hepmom = heppart['momentum']
                particle = G4PrimaryParticle(heppart['pdgid'])
                #particle.SetTrackID(heppart['index'])
                particle.Set4Momentum(hepmom[0],hepmom[1],hepmom[2],
                                      heppart['energy'])
                # Ensure mass is exact
                particle.SetMass(heppart['mass'])
                finalStateParticles.append(particle)
                print "      adding: ",heppart['pdgid']
            else:
                # ignore others
                continue
        vertex = G4PrimaryVertex(position, time)
        for particle in finalStateParticles:
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
        self._tb.edep[0] = 0
        self._tb.nstep[0] = 0
        self._tb.nq[0] = 0
        return
    
    def EndOfEventAction(self, event):
        '''Record event'''
        self._tb.ev[0] = event.GetEventID()
        self._otree.Fill()
        return

# ------------------------------------------------------------------
class MySteppingAction(G4UserSteppingAction):
    "My Stepping Action"
    def __init__(self, treebuffer, geometry, materials, edepsim):
        '''Constructor'''
        G4UserSteppingAction.__init__(self)
        self._tb = treebuffer
        self._geom = geometry
        self._materials = materials
        self._include_edepsim = edepsim
    
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
        tb.xs[istp] = prepos.x / cm
        tb.ys[istp] = prepos.y / cm
        tb.zs[istp] = prepos.z / cm 
        tb.xe[istp] = postpos.x / cm
        tb.ye[istp] = postpos.y / cm
        tb.ze[istp] = postpos.z / cm
        tb.nstep[0] += 1
        if not self._include_edepsim:
            # Stop here, unless simple energy deposition simulation is enabled
            return
        # Distribute ionization charge along track
        deltaHit = 2 * mm
        # Arrrgghhh! no python hook to non-ionizing energy deposit...
        #eIon = step.GetTotalEnergyDeposit()-step.GetNonIonizingEnergyDeposit()
        eIon = step.GetTotalEnergyDeposit()
        if track.GetDefinition().GetPDGCharge() != 0.0 and eIon > 0:
            deltaPos = postpos - prepos
            iq = tb.nq[0]
            dnq = int(deltaPos.mag() / deltaHit)           
            drq = G4ThreeVector(deltaPos)
            drq.setMag(deltaHit)
            for qidx in range(dnq):
                curpos = prepos + qidx * drq
                tb.tidq[iq+qidx] = track.GetTrackID()
                tb.pidq[iq+qidx] = track.GetDefinition().GetPDGEncoding()
                tb.sidq[iq+qidx] = istp
                tb.dq[iq+qidx] = (eIon / dnq) / MeV
                tb.xq[iq+qidx] = curpos.x / cm
                tb.yq[iq+qidx] = curpos.y / cm
                tb.zq[iq+qidx] = curpos.z / cm
            tb.nq[0] += dnq
        return

    
###################################################################

if '__main__' == __name__:
    # Run the simulation

    mysim = MySimulation()
    mysim.initialize()
    mysim.run()
    mysim.finalize()
    
