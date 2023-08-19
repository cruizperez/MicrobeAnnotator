import logging
import sys
from typing import Any, List, Optional

INFO_FORMAT = "%(asctime)s [%(levelname)s]: %(message)s"
WARNING_FORMAT = "%(asctime)s %(name)s: Function: %(funcName)s - Line: %(lineno)s [%(levelname)s]: %(message)s"


class LoggerFilter(logging.Filter):
    """Filter class that ensures handlers will log information only pertient to their level: INFO, WARNING, ERROR"""

    def __init__(self, level: int) -> None:
        self.__level = level

    def filter(self, logRecord: Any) -> bool:
        """
        Create a filter that will only log information that is at or below the level of the handler

        Args:
            logRecord (Any): Log record to be filtered

        Returns:
            bool: True if logRecord.levelno <= self.__level, False otherwise
        """
        return logRecord.levelno <= self.__level


class MicrobeAnnotatorLogger(logging.Logger):
    """Logger class for MicrobeAnnotator"""

    def __init__(self, name: str, logfile: Optional[str] = None) -> None:
        """
        Initialize a MicrobeAnnotatorLogger object

        Args:
            name (str): Name of logger
            logfile (Optional[str], optional): Path to logfile. Defaults to None.
        """
        # Validate logger name
        if not isinstance(name, str):
            raise TypeError(f"Logger name must be a string, not {type(name)}")

        super().__init__(name)

        handlers = self._create_stream_handler()
        if logfile:
            if not isinstance(logfile, str):
                raise TypeError(f"Logfile path must be a string, not {type(logfile)}")
            handlers.extend(self._create_file_handler(logfile))

        self.handlers: List[logging.Handler] = handlers
        self.level = logging.INFO

    def _create_stream_handler(self) -> List[logging.Handler]:
        """
        Create a list of stream handlers for the logger

        Returns:
            List[logging.Handler]: List of stream handlers
        """
        # Define stream handler for information
        info_handler = logging.StreamHandler(sys.stdout)
        info_handler.setLevel(logging.INFO)
        info_handler.setFormatter(logging.Formatter(INFO_FORMAT))
        info_handler.addFilter(LoggerFilter(logging.INFO))

        # Define stream handler for warnings
        warning_handler = logging.StreamHandler(sys.stdout)
        warning_handler.setLevel(logging.WARNING)
        warning_handler.setFormatter(logging.Formatter(WARNING_FORMAT))
        warning_handler.addFilter(LoggerFilter(logging.WARNING))

        # Define stream handler for errors and critical
        error_handler = logging.StreamHandler(sys.stderr)
        error_handler.setLevel(logging.ERROR)
        error_handler.setFormatter(logging.Formatter(WARNING_FORMAT))

        return [info_handler, warning_handler, error_handler]

    def _create_file_handler(self, logfile: str) -> List[logging.Handler]:
        """
        Create a list of file handlers for the logger

        Args:
            logfile (str): Path to logfile

        Returns:
            List[logging.Handler]: List of file handlers
        """
        # Define file handler for information
        info_handler = logging.FileHandler(logfile)
        info_handler.setLevel(logging.INFO)
        info_handler.setFormatter(logging.Formatter(INFO_FORMAT))
        info_handler.addFilter(LoggerFilter(logging.INFO))

        # Define file handler for warnings
        warning_handler = logging.FileHandler(logfile)
        warning_handler.setLevel(logging.WARNING)
        warning_handler.setFormatter(logging.Formatter(WARNING_FORMAT))
        warning_handler.addFilter(LoggerFilter(logging.WARNING))

        # Define file handler for errors and critical
        error_handler = logging.FileHandler(logfile)
        error_handler.setLevel(logging.ERROR)
        error_handler.setFormatter(logging.Formatter(WARNING_FORMAT))

        return [info_handler, warning_handler, error_handler]
