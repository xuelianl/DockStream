# from dockstream.utils.enums.AutodockVina_enums import AutodockVinaExecutablesEnum
from dockstream.utils.execute_external.execute import ExecutorBase
import os


# EE = AutodockVinaExecutablesEnum()


class Dock37Executor(ExecutorBase):

    def __init__(self, prefix_execution=None, binary_location=None):
        super().__init__(prefix_execution=prefix_execution, binary_location=binary_location)

    def execute(self, command: str, arguments: list, check=True, location=None):
        # check, whether a proper executable is provided
        # if command not in [EE.VINA]:
        #     raise ValueError("Parameter command must be an dictionary of the internal AutoDock Vina executable list.")

        return super().execute(command=command,
                               arguments=arguments,
                               check=check,
                               location=location)

    def is_available(self):
        if os.path.exists('/pubhome/soft/ucsfdock/3.7/bela30c'):  #### /pubhome/xli02/scripts/ucsf_dock37/bela30c
            return True
        else:
            return False
