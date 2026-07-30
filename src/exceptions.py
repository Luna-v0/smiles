from dataclasses import dataclass, field


@dataclass
class ParserException(Exception):
    """
    Exception for parser errors.

    Args:
        rule: The rule that caused the error.
        parameter: The parameter that caused the error.
        message: The error message.
        position: Character index into the input at the failure point, when known.
        expected: Legal tokens at the failure point, when known.
    """

    rule: str
    parameter: str
    message: str
    position: int | None = None
    expected: set[str] | None = None

    def __str__(self):
        """String representation of the exception."""
        return f"{self.rule}: {self.message} (parameter: {self.parameter})"


@dataclass
class RingSemanticsException(ParserException):
    """
    Exception for ring-bond semantics violations (OpenSMILES §3.4/§3.6).

    Ring-number matching needs unbounded state, so these rules live outside
    the context-free grammar: unclosed rings, self ring-bonds, duplicate
    bonds via ring closure, and mismatched ring-closure bond orders.
    """
