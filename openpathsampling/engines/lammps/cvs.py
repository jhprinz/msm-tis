import openpathsampling as paths

class LAMMPSComputeCV(paths.FunctionCV):
    """CV based on a LAMMPS ``compute`` statement"""
    # this is ridiculously ugly, but couldn't get it to work with normal
    # _eval and inherit from CollectiveVariable
    def __init__(self, name, groupid_style_args, extract_style,
                 extract_type, engine, cv_time_reversible=False):
        super(LAMMPSComputeCV, self).__init__(
            name=name,
            f=self._lammps_eval,
            cv_time_reversible=cv_time_reversible
        )
        self.compute_command = "compute {name} {groupid_style_args}".format(
            name=name,
            groupid_style_args=groupid_style_args
        )
        self.groupid_style_args = groupid_style_args
        self.engine = engine
        self.extract_style = extract_style
        self.extract_type = extract_type
        self.engine.command(self.compute_command)
        self._first_run = False

    def to_dict(self):
        return {
            'name': self.name,
            'groupid_style_args': self.groupid_style_args,
            'extract_style': self.extract_style,
            'extract_type': self.extract_type,
            'engine': self.engine,
            'cv_time_reversible': self.cv_time_reversible
        }

    @classmethod
    def from_dict(cls, dct):
        args = ["name", "groupid_style_args", "extract_style",
                "extract_type", "engine", "cv_time_reversible"]
        return cls(**{key: dct[key] for key in args})

    def _lammps_eval(self, snapshot):
        cur_snap = self.engine.current_snapshot
        if snapshot != cur_snap:
            self.engine.current_snapshot = snapshot
            self.engine.command("run 0")
        elif not self._first_run:
            self.engine.command("run 0")
            self._first_run = True

        result = self.engine.lammps.extract_compute(self.name,
                                                    self.extract_style,
                                                    self.extract_type)

        if snapshot != cur_snap:
            self.engine.current_snapshot = cur_snap

        return result
