from .structure import Chemistry
from typing import List
import re


def _parse_electron_configuration(electron_config: List[str]) -> dict[tuple[int, str], int]:
    """
    Parses a list of electron configuration strings (e.g., ['1s2', '2s2', '2p6'])
    into a dictionary mapping (principal_quantum_number, subshell_type) to electron count.
    """
    parsed_config = {}
    for orbital_str in electron_config:
        match = re.match(r'(\d+)([spdf])(\d+)', orbital_str)
        if match:
            n = int(match.group(1))
            subshell = match.group(2)
            electrons = int(match.group(3))
            parsed_config[(n, subshell)] = electrons
    return parsed_config


chemistry = Chemistry()