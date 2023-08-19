from unittest.mock import Mock, patch

import pytest
from bs4 import BeautifulSoup

from microbeannotator.utils.scraper import get_web_content


class TestUtils:
    @pytest.fixture
    def mock_requests_get(self) -> Mock:
        """Create a mock response for requests.get"""
        mock_response = Mock()
        mock_response.status_code = 200
        mock_response.content = b"<html><body>Mock HTML Content</body></html>"
        return mock_response

    def test_get_web_content_success(self, mock_requests_get: Mock) -> None:
        """Test that get_web_content returns a BeautifulSoup object and that the content is correct."""
        with patch("requests.get", return_value=mock_requests_get):
            url = "http://example.com"
            result = get_web_content(url)

            assert isinstance(result, BeautifulSoup)
            assert str(result) == "<html><body>Mock HTML Content</body></html>"

    def test_get_web_content_failure(self, mock_requests_get: Mock) -> None:
        """Test that get_web_content raises a ConnectionError when the response status code is not 200."""
        mock_requests_get.status_code = 404  # Simulate a failure response
        with patch("requests.get", return_value=mock_requests_get):
            url = "http://example.com"
            with pytest.raises(ConnectionError) as error:
                get_web_content(url)
                assert "Bad Request - The request was malformed or invalid." in str(error.value)
