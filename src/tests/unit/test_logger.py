import logging
from unittest.mock import Mock

from microbeannotator.logger import LoggerFilter


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
