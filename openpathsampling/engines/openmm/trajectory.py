import mdtraj as md
import numpy as np
import simtk.unit as u

from tools import FileEngine
from openpathsampling.netcdfplus import StorableObject, LoaderProxy
import openpathsampling as paths
from openpathsampling.engines import Trajectory
from snapshot import Snapshot
from topology import Topology, MDTrajTopology

from uuid import UUID

class ExternalMDTrajectory(Trajectory):
    """
    A trajectory stored in an external file

    The main difference is that these snapshots are only available in a
    block and snapshots will be generated if needed. Also snapshots
    are not saved, but the UUID counter is advanced by as many counts
    as necessary. Therefore you cannot use these trajectories without
    UUID support
    """

    def __init__(self, filename):
        super(ExternalMDTrajectory, self).__init__()
        self.filename = filename

        f = md.load(filename)

        # advance the UUID for each virtual snapshot
        StorableObject.CREATION_COUNT += 2 * len(f)

        self.extend([
                        LoaderProxy(self, UUID(int=int(self.__uuid__) + 2 + idx * 2))
                        for idx in range(len(f))
                        ])

        self._f = f
        self._topology = MDTrajTopology(f.topology)


    def __getitem__(self, index):
        # Allow for numpy style selection using lists
        if hasattr(index, '__iter__'):
            ret = [list.__getitem__(self, i) for i in index]
        elif isinstance(index, UUID):
            ret = list.__getitem__(self, (int(index) - int(self.__uuid__) - 2) / 2)
        else:
            ret = list.__getitem__(self, index)

        if type(ret) is list:
            ret = Trajectory(ret)
        elif hasattr(ret, '_idx'):
            f = self._f
            velocities = np.zeros(f.xyz[0].shape)

            snapshot = Snapshot.construct(
                coordinates=u.Quantity(f.xyz[0], u.nanometers),
                box_vectors=u.Quantity(f.unitcell_vectors[0], u.nanometers),
                velocities=u.Quantity(velocities, u.nanometers / u.picoseconds),
                engine=FileEngine(self._topology, self.filename)
            )

            snapshot.__uuid__ = ret._idx

            return snapshot

        return ret

