from dataclasses import dataclass

@dataclass
class ParserException(Exception):
    """
    Exception for parser errors.

    Args:
        rule: The rule that caused the error.
        parameter: The parameter that caused the error.
        message: The error message.
    """

    rule: str
    parameter: str
    message: str
    
    def __str__(self):
        """String representation of the exception."""
        return f"{self.rule}: {self.message} (parameter: {self.parameter})"