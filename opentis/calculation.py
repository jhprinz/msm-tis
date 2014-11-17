from __future__ import print_function
from __future__ import absolute_import

import os
import sys
from .globalstate import GlobalState
from .pathmover import PathMover, MoveDetails, ReplicaExchange
from .trajectory import Sample

class Calculation(object):

    calc_name = "Calculation"

    def __init__(self, storage, simulator=None, ensembles=None, 
                 movers=None):
        self.storage = storage
        self.simulator = simulator
        self.movers = movers
        self.globalstate = GlobalState(ensembles)

    # TODO: this should be a property
    def set_replicas(self, replicas):
        if type(replicas) is dict:
            for ensemble, trajectory in dict.iteritems():
                self.globalstate[ensemble] = trajectory
        elif type(replicas) is list:
            ensembles = self.globalstate.ensembles
            for idx, trajectory in enumerate(replicas):
                self.globalstate[ensembles[idx]] = trajectory

    def run(self, nsteps):
        print("Running an empty calculation? Try a subclass, maybe!")

class BootstrapEnsembleChangeMove(PathMover):
    def move(self, trajectory, ensemble):
        details = MoveDetails()
        details.inputs = [trajectory]
        details.mover = self
        details.final = trajectory
        details.success = True
        details.acceptance = 1.0
        details.result = trajectory

        sample = Sample(
            trajectory=details.result,
            mover=self,
            ensemble=ensemble,
            details=details
        )

        return sample

class Bootstrapping(Calculation):
    """The ensembles for the Bootstrapping calculation must be one ensemble
    set, in increasing order."""

    calc_name = "Bootstrapping"


    def run(self, nsteps):

        swapmove = ReplicaExchange()
        bootstrapmove = BootstrapEnsembleChangeMove()

        ens_num = 0
        failsteps = 0
        # if we fail nsteps times in a row, kill the job
        while ens_num < self.globalstate.size - 1 and failsteps < nsteps:
            print("Trying move in ensemble", ens_num)
            # Generate Samples

            samples = [ self.movers[ens_idx].move(self.globalstate[ens_idx]) for ens_idx in range(ens_num, ens_num + 1) ]

#            if ens_num > 1:
#                ex_samples = swapmove.move(self.globalstate[ens_num - 1], self.globalstate[ens_num - 2], self.globalstate.ensembles[ens_num - 1], self.globalstate.ensembles[ens_num - 2])
#                samples = samples + ex_samples

            # Generate new globalstate using only the one sample
            globalstate = self.globalstate.move(samples)

            # Now save all samples
            self.globalstate.save_samples(self.storage)

            # update to new globalstate
            self.globalstate = globalstate

            if ens_num < self.globalstate.size:
                # We can try to switch to the next ensemble

                # Check if the new trajectory (still in the old ensemble) would
                # fir in the next ensemble
                if self.globalstate.ensembles[ens_num + 1](self.globalstate[ens_num]):
                    # Yes, so apply the BootStrapMove and generate a new sample in the next ensemble
                    sample = bootstrapmove.move(self.globalstate[ens_num], self.globalstate.ensembles[ens_num + 1])
                    globalstate = self.globalstate.move([sample])

                    # Now save all samples
                    self.globalstate.save_samples(self.storage)

                    # update to new globalstate
                    self.globalstate = globalstate

                    failsteps = 0
                    ens_num += 1

            failsteps += 1