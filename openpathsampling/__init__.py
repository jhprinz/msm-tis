from pathsimulator import (
    PathSimulator, FullBootstrapping, Bootstrapping, PathSampling, MCStep
)

from ensemble import (
    Ensemble, EnsembleCombination, EnsembleFactory, EntersXEnsemble,
    EmptyEnsemble, ExitsXEnsemble, FullEnsemble, PartInXEnsemble,
    AllInXEnsemble, AllOutXEnsemble, WrappedEnsemble,
    SuffixTrajectoryEnsemble, PrefixTrajectoryEnsemble,
    PartOutXEnsemble, LengthEnsemble, NegatedEnsemble,
    ReversedTrajectoryEnsemble, SequentialEnsemble, VolumeEnsemble,
    SequentialEnsemble, IntersectionEnsemble, UnionEnsemble,
    SymmetricDifferenceEnsemble, RelativeComplementEnsemble,
    SingleFrameEnsemble, MinusInterfaceEnsemble, TISEnsemble,
    OptionalEnsemble, join_ensembles
)

from snapshot import Snapshot, Momentum, ToySnapshot, AbstractSnapshot
from snapshot_content import Configuration, Momentum
from trajectory import Trajectory
from sample import Sample, SampleSet

from collectivevariable import CV_Function, CV_MDTraj_Function, CV_MSMB_Featurizer, \
    CV_Volume, CollectiveVariable

from pathmover import (
    RandomChoiceMover, PathMover, ConditionalSequentialMover,
    PartialAcceptanceSequentialMover, BackwardShootMover, ForwardShootMover,
    BackwardExtendMover, ForwardExtendMover, MinusMover,
    SingleReplicaMinusMover, PathMoverFactory, PathReversalMover,
    ReplicaExchangeMover, EnsembleHopMover, ReplicaIDChangeMover,
    SequentialMover, ConditionalMover,
    PathSimulatorMover, PathReversalSet, NeighborEnsembleReplicaExchange,
    SampleMover, StateSwapMover, FinalSubtrajectorySelectMover, EngineMover,
    FirstSubtrajectorySelectMover, MultipleSetMinusMover,
    OneWayShootingMover, RandomSubtrajectorySelectMover, SubPathMover,
    EnsembleFilterMover, SelectionMover, FirstAllowedMover,
    LastAllowedMover, OneWayExtendMover, SubtrajectorySelectMover
)

from shooting import ShootingPointSelector, UniformSelector, \
    GaussianBiasSelector, FirstFrameSelector, FinalFrameSelector

from dynamics_engine import DynamicsEngine

from openmm_engine import OpenMMEngine

# from lammps_engine import LammpsEngine

from volume import (Volume, VolumeCombination, VolumeFactory, VoronoiVolume, 
    EmptyVolume, FullVolume, CVRangeVolume, CVRangeVolumePeriodic,
    IntersectionVolume, UnionVolume, SymmetricDifferenceVolume,
    RelativeComplementVolume, join_volumes
)

from tools import empty_snapshot_from_openmm_topology, snapshot_from_pdb, \
    to_openmm_topology, trajectory_from_mdtraj

from tools import simtk_units_from_md_snapshot

from topology import ToyTopology, Topology, MDTrajTopology

from toy_dynamics.toy_pes import Gaussian, HarmonicOscillator, LinearSlope, \
    OuterWalls, Toy_PES, Toy_PES_Add, Toy_PES_Sub

from toy_dynamics.toy_engine import ToyEngine

from toy_dynamics.toy_integrators import LangevinBAOABIntegrator, \
    LeapfrogVerletIntegrator

from analysis.tis_analysis import (
    TISTransition, Transition, TPSTransition
)

from analysis.move_scheme import MoveScheme, DefaultScheme, LockedMoveScheme

from analysis.network import (
    MSTISNetwork, TransitionNetwork, MISTISNetwork, TPSNetwork
)

from analysis.replica_network import (
    ReplicaNetwork, trace_ensembles_for_replica,
    trace_replicas_for_ensemble, condense_repeats,
    ReplicaNetworkGraph
)

from live_visualization import LiveVisualization

from pathmover import Details, MoveDetails, SampleDetails

from pathmovechange import (
    EmptyPathMoveChange, ConditionalSequentialPathMoveChange,
    PathMoveChange, PartialAcceptanceSequentialPathMoveChange,
    RandomChoicePathMoveChange, SamplePathMoveChange,
    SequentialPathMoveChange,  KeepLastSamplePathMoveChange,
    FilterSamplesPathMoveChange,
    PathSimulatorPathMoveChange, AcceptedSamplePathMoveChange,
    RejectedSamplePathMoveChange, SubPathMoveChange,
    FilterByEnsemblePathMoveChange
)

from storage.storage import Storage, AnalysisStorage

def git_HEAD(): # pragma: no cover
    from subprocess import check_output
    import os.path
    git_dir = os.path.dirname(os.path.realpath(__file__))
    return check_output(["git", "-C", git_dir, "rev-parse", "HEAD"])[:-1]
    # chops the newline at the end

try:
    import version
except ImportError: # pragma: no cover
    pass
