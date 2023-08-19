import logging
from pathlib import Path
from unittest.mock import Mock

import pytest

from microbeannotator.logger import LoggerFilter, MicrobeAnnotatorLogger


class TestLoggerFilter:
    """Class to test the LoggerFilter class"""

    def test_filter_below_level(self) -> None:
        """Test that a log record with a level below the filter's level passes through the filter"""
        # Create a mock log record with a level below the filter's level
        log_record = Mock()
        log_record.levelno = logging.INFO  # Log level INFO

        filter_instance = LoggerFilter(level=logging.WARNING)  # Filter level set to WARNING
        result = filter_instance.filter(log_record)

        assert result is True  # Log record should pass through the filter

    def test_filter_above_level(self) -> None:
        """Test that a log record with a level above the filter's level is filtered out"""
        # Create a mock log record with a level above the filter's level
        log_record = Mock()
        log_record.levelno = logging.ERROR  # Log level ERROR

        filter_instance = LoggerFilter(level=logging.WARNING)  # Filter level set to WARNING
        result = filter_instance.filter(log_record)

        assert result is False  # Log record should be filtered out

    def test_filter_equal_level(self) -> None:
        """Test that a log record with a level equal to the filter's level passes through the filter"""
        # Create a mock log record with a level equal to the filter's level
        log_record = Mock()
        log_record.levelno = logging.WARNING  # Log level WARNING

        filter_instance = LoggerFilter(level=logging.WARNING)  # Filter level set to WARNING
        result = filter_instance.filter(log_record)

        assert result is True  # Log record should pass through the filter

    def test_filter_higher_level(self) -> None:
        """Test that a log record with a higher log level than the filter's level is filtered out"""
        # Create a mock log record with a higher log level than the filter's level
        log_record = Mock()
        log_record.levelno = logging.ERROR  # Log level ERROR

        filter_instance = LoggerFilter(level=logging.WARNING)  # Filter level set to WARNING
        result = filter_instance.filter(log_record)

        assert result is False  # Log record should be filtered out


class TestMicrobeAnnotatorLogger:
    def test_init_with_valid_name(self) -> None:
        """Test that a MicrobeAnnotatorLogger instance is created with a valid name"""
        logger = MicrobeAnnotatorLogger("test_logger")
        assert isinstance(logger, MicrobeAnnotatorLogger)
        assert logger.name == "test_logger"

    def test_init_with_invalid_name_type(self) -> None:
        """Test that a TypeError is raised when a MicrobeAnnotatorLogger is created with an invalid name type"""
        with pytest.raises(TypeError):
            MicrobeAnnotatorLogger(123)  # type: ignore

    def test_init_with_valid_logfile(self, tmp_path: Path) -> None:
        """Test that a MicrobeAnnotatorLogger instance is created with a valid logfile"""
        log_file = tmp_path / "test.log"
        logger = MicrobeAnnotatorLogger("test_logger", logfile=str(log_file))
        assert len(logger.handlers) == 6
        assert isinstance(logger.handlers[0], logging.StreamHandler)
        assert isinstance(logger.handlers[-1], logging.FileHandler)
        assert logger.handlers[-1].baseFilename == str(log_file)

    def test_init_with_invalid_logfile_type(self) -> None:
        """Test that a TypeError is raised when a MicrobeAnnotatorLogger is created with an invalid logfile type"""
        with pytest.raises(TypeError):
            MicrobeAnnotatorLogger("test_logger", logfile=123)  # type: ignore

    def test_create_stream_handler(self) -> None:
        """Test that a MicrobeAnnotatorLogger instance is created with a valid logfile"""
        logger = MicrobeAnnotatorLogger("test_logger")
        handlers = logger._create_stream_handler()
        assert len(handlers) == 3
        assert all([isinstance(handler, logging.StreamHandler)] for handler in handlers)

    def test_create_file_handler(self, tmp_path: Path) -> None:
        """Test that a MicrobeAnnotatorLogger instance is created with a valid logfile"""
        log_file = tmp_path / "test.log"
        logger = MicrobeAnnotatorLogger("test_logger")
        handlers = logger._create_file_handler(str(log_file))
        assert len(handlers) == 3
        assert all([isinstance(handler, logging.FileHandler)] for handler in handlers)
