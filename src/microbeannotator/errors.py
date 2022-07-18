from typing import Any, Optional


class AttributeTypeError(Exception):
    """
    Exception raised when an object is not of the expected type.
    """
    def __init__(self, attribute: Any, attr_type: type, message: Optional[str] = None) -> None:
        if not message:
            self.message = f"Attribute <{attribute}> is not of type <{attr_type}>"
        else:
            self.message = message
        super().__init__(self.message)