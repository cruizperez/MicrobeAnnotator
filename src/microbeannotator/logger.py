import logging
import sys
from typing import Any, List, Optional

from microbeannotator.constants import INFO_FORMAT, WARNING_FORMAT
from microbeannotator.errors import AttributeTypeError


class LoggingFilter(logging.Filter):
    """
    Filter class to ensure handlers will only log information relevant to their level:
    INFO, WARNING, ERROR.
    """

    def __init__(self, level: int) -> None:
        self.__level = level

    def filter(self, logRecord: Any) -> bool:
        """
        Creates a filter level so handlers will only log information pertinent to their own level.

        Args:
            logRecord (Any): Record to be processed.

        Returns:
            bool: Record should be processed (True if record level is below filter level).
        """
        return logRecord.levelno <= self.__level


class MicrobeAnnotatorLogger(logging.Logger):
    def __init__(self, logger_name: str, logfile: Optional[str] = None) -> None:
        """
        Calls the constructor of logging.Logger to create a new logger.
        Modifies the handlers and the level.

        Args:
            logger_name (str): Name of logger.
            logfile (Optional[str], optional): File to store logs. Defaults to None.

        Raises:
            AttributeTypeError: Exception raised when logger_name or logfile are not str.
        """
        if not isinstance(logger_name, str):
            raise AttributeTypeError(logger_name, str)

        super().__init__(logger_name)

        handlers = self._create_stream_handler()
        if logfile:
            if not isinstance(logfile, str):
                raise AttributeTypeError(logfile, str)
            handlers.extend(self._create_file_handler(logfile))

        self.handlers = handlers
        self.level = logging.INFO

    def _create_stream_handler(self) -> List[logging.Handler]:
        """
        Creates list of handlers that output to stdout/stderr.

        Returns:
            List[logging.Handler]: List of StreamHandler (stdout & stderr).
        """
        # Define stream handler for INFO
        info_handler = logging.StreamHandler(sys.stdout)
        info_handler.setLevel(logging.INFO)
        info_handler.setFormatter(logging.Formatter(INFO_FORMAT))
        info_handler.addFilter(LoggingFilter(logging.INFO))

        # Define stream handler for WARNING
        warning_handler = logging.StreamHandler(sys.stdout)
        warning_handler.setLevel(logging.WARNING)
        warning_handler.setFormatter(logging.Formatter(WARNING_FORMAT))
        warning_handler.addFilter(LoggingFilter(logging.WARNING))

        # Define stream handler for ERROR
        error_handler = logging.StreamHandler(sys.stdout)
        error_handler.setLevel(logging.ERROR)
        error_handler.setFormatter(logging.Formatter(WARNING_FORMAT))
        error_handler.addFilter(LoggingFilter(logging.ERROR))

        return [info_handler, warning_handler, error_handler]

    def _create_file_handler(self, logfile: str) -> List[logging.Handler]:
        """
        Creates list of handlers that output to logfile.

        Args:
            logfile (str): Path to file where logs are stored.

        Returns:
            List[logging.Handler]: List of FileHandler (stores in logfile).
        """
        # Define file handler for INFO
        info_handler = logging.FileHandler(logfile, mode="a")
        info_handler.setLevel(logging.INFO)
        info_handler.setFormatter(logging.Formatter(INFO_FORMAT))
        info_handler.addFilter(LoggingFilter(logging.INFO))

        # Define file handler for WARNING
        warning_handler = logging.FileHandler(logfile, mode="a")
        warning_handler.setLevel(logging.WARNING)
        warning_handler.setFormatter(logging.Formatter(WARNING_FORMAT))
        warning_handler.addFilter(LoggingFilter(logging.WARNING))

        # Define file handler for ERROR
        error_handler = logging.FileHandler(logfile, mode="a")
        error_handler.setLevel(logging.ERROR)
        error_handler.setFormatter(logging.Formatter(WARNING_FORMAT))
        error_handler.addFilter(LoggingFilter(logging.ERROR))

        return [info_handler, warning_handler, error_handler]
