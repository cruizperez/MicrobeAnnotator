from enum import Enum

import requests
from bs4 import BeautifulSoup

from microbeannotator import logger


class HttpStatusCode(Enum):
    OK = 200
    CREATED = 201
    NO_CONTENT = 204
    MOVED_PERMANENTLY = 301
    FOUND = 302
    SEE_OTHER = 303
    NOT_MODIFIED = 304
    BAD_REQUEST = 400
    UNAUTHORIZED = 401
    FORBIDDEN = 403
    NOT_FOUND = 404
    METHOD_NOT_ALLOWED = 405
    TOO_MANY_REQUESTS = 429
    INTERNAL_SERVER_ERROR = 500
    BAD_GATEWAY = 502
    SERVICE_UNAVAILABLE = 503
    GATEWAY_TIMEOUT = 504

    def meaning(self) -> str:
        meanings = {
            self.OK: "OK - The request was successful.",
            self.CREATED: "Created - The request was successful, and a new resource was created.",
            self.NO_CONTENT: "No Content - The request was successful, but there is no content to send.",
            self.MOVED_PERMANENTLY: "Moved Permanently - The requested resource has been permanently moved.",
            self.FOUND: "Found - The requested resource has been temporarily moved.",
            self.SEE_OTHER: "See Other - The requested resource has been temporarily moved.",
            self.NOT_MODIFIED: "Not Modified - The client's cached copy of the resource is still valid.",
            self.BAD_REQUEST: "Bad Request - The request was malformed or invalid.",
            self.UNAUTHORIZED: "Unauthorized - The request requires user authentication.",
            self.FORBIDDEN: "Forbidden - The client does not have permission to access the resource.",
            self.NOT_FOUND: "Not Found - The requested resource was not found.",
            self.METHOD_NOT_ALLOWED: "Method Not Allowed - The HTTP method is not allowed for the resource.",
            self.TOO_MANY_REQUESTS: "Too Many Requests - Client has sent too many requests in a short period.",
            self.INTERNAL_SERVER_ERROR: "Internal Server Error - An unexpected error occurred on the server.",
            self.BAD_GATEWAY: "Bad Gateway - The server received an invalid response from an upstream server.",
            self.SERVICE_UNAVAILABLE: "Service Unavailable - Server is temporarily unable to handle the request.",
            self.GATEWAY_TIMEOUT: "Gateway Timeout - Server didn't receive timely response.",
        }

        return meanings.get(self.value)


def get_web_content(url: str) -> BeautifulSoup:
    """Returns the content of a web page as a BeautifulSoup object."""
    response = requests.get(url)
    if response.status_code == 200:
        content = response.content
        return BeautifulSoup(content, "html.parser")
    else:
        try:
            status_code_info = HttpStatusCode(response.status_code)
            logger.error(f"Connection Error {response.status_code}: {status_code_info.meaning()}")
        except ValueError:
            logger.error(f"Connection Error {response.status_code}: Unknown status code.")
        raise ConnectionError
