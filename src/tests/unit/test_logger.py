import logging
from pathlib import Path
from typing import Union

import pytest

from microbeannotator.errors import AttributeTypeError
from microbeannotator.logger import LoggingFilter, MicrobeAnnotatorLogger


@pytest.fixture()
def logger() -> MicrobeAnnotatorLogger:
    """
    Setup a new mock logger
    """
    return MicrobeAnnotatorLogger("test_logger")


class TestLoggingFilter:
    @pytest.mark.parametrize(
        "logging_level, expected_level",
        [
            (logging.NOTSET, 0),
            (logging.DEBUG, 10),
            (logging.INFO, 20),
            (logging.WARNING, 30),
            (logging.ERROR, 40),
            (logging.CRITICAL, 50),
            (15, 15),
            (25, 25),
            (36, 36),
        ],
    )
    def test_logging_filter_level(self, logging_level: int, expected_level: int) -> None:
        """
        Test expected logging level

        Args:
            logging_level (int): Logging level to test.
            expected_level (int): Expected value for logging level.
        """
        filter = LoggingFilter(logging_level)
        assert filter._LoggingFilter__level == expected_level  # type: ignore

    @pytest.mark.parametrize(
        "logger_filter_level, logrecord_level, process_logrecord",
        [
            (logging.INFO, logging.INFO, True),
            (logging.INFO, logging.WARNING, False),
            (logging.INFO, logging.ERROR, False),
            (logging.WARNING, logging.INFO, True),
            (logging.WARNING, logging.WARNING, True),
            (logging.WARNING, logging.ERROR, False),
            (logging.ERROR, logging.INFO, True),
            (logging.ERROR, logging.WARNING, True),
            (logging.ERROR, logging.ERROR, True),
        ],
    )
    def test_logging_filter_process(
        self,
        logger_filter_level: int,
        logrecord_level: int,
        process_logrecord: bool,
    ) -> None:
        # Arrange
        mock_logrecord = logging.LogRecord(
            name="test",
            level=logrecord_level,
            pathname="",
            lineno=0,
            msg="test message",
            args=(),
            exc_info=None,
        )
        # Act
        logfilter = LoggingFilter(logger_filter_level)
        # Assert
        assert logfilter.filter(mock_logrecord) == process_logrecord


class TestMicrobeAnnotatorLogger:
    def test_create_stream_handler(self, logger: MicrobeAnnotatorLogger) -> None:
        """Test the number, level and type of the stream handlers created by the logger"""
        stream_handler_list = logger._create_stream_handler()
        assert len(stream_handler_list) == 3
        for handler, expected_level in zip(stream_handler_list, [20, 30, 40]):
            assert type(handler) == logging.StreamHandler
            assert handler.level == expected_level

    def test_create_file_handler(self, tmp_path: Path, logger: MicrobeAnnotatorLogger) -> None:
        """Test the number, level and type of the file handlers created by the logger"""
        expected_logfile = f"{tmp_path}/test_logfile.txt"
        stream_handler_list = logger._create_file_handler(expected_logfile)
        assert Path(expected_logfile).is_file()
        assert len(stream_handler_list) == 3
        for handler, expected_level in zip(stream_handler_list, [20, 30, 40]):
            assert type(handler) == logging.FileHandler
            assert handler.level == expected_level

    def test_logger_creation(self, tmp_path: Path) -> None:
        """Assert number of handlers increases when log_file is passed"""
        logger = MicrobeAnnotatorLogger("test")
        assert len(logger.handlers) == 3
        logger = MicrobeAnnotatorLogger("test", logfile=f"{tmp_path}/test_logfile.txt")
        assert len(logger.handlers) == 6

    @pytest.mark.parametrize(
        "logname, log_file",
        [
            (10000, "/tmp/test_file.txt"),
            ("logger_name", True),
            (10000, False),
        ],
    )
    def test_logger_creation_bad_name(
        self, logname: Union[str, int], log_file: Union[str, bool]
    ) -> None:
        """Test if error is raised when a non-string is passed as the logger name."""
        with pytest.raises(AttributeTypeError):
            MicrobeAnnotatorLogger(logname, logfile=log_file)  # type: ignore
