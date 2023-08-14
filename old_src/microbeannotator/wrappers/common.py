import sys
from abc import ABC, abstractmethod
from dataclasses import dataclass
from subprocess import Popen
from typing import Dict, List


class WrapperAlgorithmBase(ABC):
    """
    Abstract base class for all algorithm classes.

    Attributes:
        flags (Dict): A dictionary that maps flags to attributes.
    Methods:
        to_string: generate a shell command for this program.
    """

    @property
    @abstractmethod
    def flags(self) -> Dict:
        """
        Get dictionary with execution flags.

        Returns:
            Dict: Execution flags.
        """

    def to_string(self) -> str:
        """
        Convert command and flags to string.

        Returns:
            str: String with execution command and flags.
        """


@dataclass
class ExecutableWrapper:
    """
    Base class for all executable wrappers.

    Methods:
        _start_cmd: Executes a command using subprocess.Popen
    """

    def _start_cmd(self, cmds: str, detach: bool = False) -> Popen:
        """
        Executes a command using subprocess.Popen.

        Args:
            cmds (str): Command to be executed.
            detach (bool, optional): Detach from process. Defaults to False.

        Returns:
            Popen: Process executed.
        """
        print(cmds)
        process: Popen = Popen(cmds, stdout=sys.stdout, universal_newlines=True, shell=True)
        if not detach:
            process.wait()

        return process


def get_cmd_with_flags(flags: Dict, program_name: str = "") -> str:
    """
    Helper function to convert a dictionary into a command string.

    Args:
        flags (Dict): Dictionary with flags to add to the command.
        program_name (str, optional): Name of program to execute. Defaults to "".

    Returns:
        str: Command string to execute.
    """
    _cmds: List[str] = [program_name] if program_name != "" else []
    for flag, value in flags.items():
        if isinstance(value, bool):
            _cmds.append(flag)
        else:
            _cmds.append(flag)
            _cmds.append(str(value))

    return " ".join(_cmds)
