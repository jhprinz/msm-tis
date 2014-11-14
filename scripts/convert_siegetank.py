'''
Created on 03.09.2014
@author: jan-hendrikprinz
'''


import base64 # This may require Python 3
import os
import gzip
import siegetank
import datetime as dt

import mdtraj as md
import os

import simtk.unit as u
import numpy as np
from opentis.snapshot import Snapshot, Configuration, Momentum
from opentis.trajectory import Trajectory

from opentis.storage import Storage

def data_from_target(target):
    data = {
        'steps_per_frame': target.options['steps_per_frame'],
        'description' : target.options['description'],
        'date' : dt.datetime.fromtimestamp(target.creation_date),
        'date_string' : dt.datetime.fromtimestamp(target.creation_date).strftime("%B %d, %Y"),
        'engines' : target.engines,
        'stage' : target.stage,
        'weight' : target.weight,
        'id' : target.id,
        'owner' : target.owner,
        'shards' : target.shards,
        'uri' : target.uri
        }
    return data

def data_from_stream(stream):
    data = {
        'active' : stream.active,
        'error_count' : stream.error_count,
        'frames' : stream.frames,
        'status' : stream.status,
        'id' : stream.id,
        'uri' : stream.uri
        }
    return data

def mdframe_to_snapshot(frame):
    c = Configuration(coordinates=u.Quantity(frame._xyz[0], u.nanometers), topology=frame.topology)
    s = Snapshot(configuration=c)
    return s

def md_to_trajectory(mdtrajectory):
    mdtrajectory = md.Trajectory()
    t = Trajectory()
    for frame in mdtrajectory:
        c = Configuration(
            coordinates=u.Quantity(frame._xyz[0], u.nanometers),
            box_vectors=u.Quantity(frame.unitcell_vectors, u.nanometers),
            topology=frame.topology
        )
        s = Snapshot(configuration=c)
        t.append(s)

    return t

if __name__ == '__main__':
    # Need a more secure way to store and load this.
    my_token = os.environ["SIEGETANK_TOKEN"]
    target_token = '309216fc-08d1-4b46-822a-27d9b92fadab'

    TEMP_FOLDER = 'temp'

    siegetank.login(my_token)

    TOP_FILE = os.path.join('/Users','jan-hendrikprinz','Studium','git','cragganmore', 'siegetank', '2D6', 'RUNS', 'RUN1', 'minimized.pdb')
    top_file = md.load(TOP_FILE)

    heme_indices = [ a.index for a in top_file.topology.atoms if a.residue.name == 'HEM']
    counterion_indices = [ a.index for a in top_file.topology.atoms if a.residue.name == 'Na+']
    solvent_indices = [ a.index for a in top_file.topology.atoms if a.residue.name == 'HOH']
    protein_indices = [ a.index for a in top_file.topology.atoms if a.residue.name != 'Na+' and a.residue.name != 'HEM' and a.residue.name != 'HOH']

    relevant_indices = protein_indices + heme_indices

    top_relevant = md.load(TOP_FILE)
    top_relevant.restrict_atoms(relevant_indices)

    solute_topology = top_relevant.restrict_atoms(relevant_indices)

    target = siegetank.load_target(target_token)

    storage = Storage('siegetank.nc', mode='w', topology_file=TOP_FILE)

    for stream in target.streams[0]:
        data_folder = os.path.join(TEMP_FOLDER, stream.id+'_data')
        stream.sync(data_folder)

        parts_folders = os.walk(data_folder).next()[1]
        sorted_folders = sorted(parts_folders, key=lambda value: int(value))

        if len(sorted_folders) > 0:
            stream_files = [os.path.join(data_folder, folder, 'frames.xtc') for folder in sorted_folders]

            traj = md.load(stream_files, top = top_file)
            traj.restrict_atoms(relevant_indices)

            storage.trajectory.save(md_to_trajectory(traj))

    pass