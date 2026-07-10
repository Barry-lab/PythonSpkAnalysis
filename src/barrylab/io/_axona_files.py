"""Low-level readers for the Axona / DACQ file family.

This module is the byte-level foundation of the DACQ front-end loader. Every
Axona file (``.set``, ``.pos``, tetrode ``.N``, ``.eeg``/``.egf``) is a text
header of ``key value`` lines terminated by the literal marker ``data_start``,
optionally followed by a big-endian binary payload ending before ``data_end``.

Ported from the MATLAB ``read_binary_data.m``, ``find_word.m`` and
``key_value.m``. Kept deliberately small and dependency-light (numpy only) so
it is easy to audit against the originals.
"""

from typing import List, Optional, Tuple

import numpy as np

_DATA_START = b"data_start"
_DATA_END = b"data_end"

# A parsed header: an ordered list of (key, value) string pairs, mirroring the
# 2-column cell array the MATLAB toolbox passes around.
Header = List[Tuple[str, str]]


def _parse_header_text(text: str) -> Header:
    """Split raw header text into ordered (key, value) pairs.

    Each line is ``key`` followed by whitespace and the remainder as the value
    (values may contain spaces, e.g. ``trial_date Tuesday, 14 Dec 2010``).
    """
    header = []  # type: Header
    for line in text.replace("\r\n", "\n").replace("\r", "\n").split("\n"):
        if not line.strip():
            continue
        parts = line.split(None, 1)
        key = parts[0]
        value = parts[1] if len(parts) > 1 else ""
        header.append((key, value.rstrip()))
    return header


def read_axona_file(path: str, has_binary: bool) -> Tuple[Header, Optional[bytes]]:
    """Read an Axona file, returning its header and (optional) binary payload.

    Parameters
    ----------
    path
        Full path including extension.
    has_binary
        ``True`` for ``.pos``/tetrode/``.eeg`` files (header + binary segment),
        ``False`` for ``.set`` (pure text, no binary).

    Returns
    -------
    (header, data)
        ``data`` is the raw bytes strictly between the ``data_start`` marker and
        the ``\\r\\n`` that precedes ``data_end`` (``None`` for text-only files).
    """
    with open(path, "rb") as fh:
        raw = fh.read()

    if not has_binary:
        return _parse_header_text(raw.decode("latin-1")), None

    i0 = raw.find(_DATA_START)
    if i0 < 0:
        raise ValueError("could not find data_start marker in %s" % path)
    header = _parse_header_text(raw[:i0].decode("latin-1"))

    # Binary begins immediately after the 10-character 'data_start' token and
    # ends 2 bytes before 'data_end' (dropping the preceding CR/LF), matching
    # read_binary_data.m's ``dsmarker+1 : demarker-3`` slice.
    data_begin = i0 + len(_DATA_START)
    i1 = raw.find(_DATA_END, data_begin)
    data_end = (i1 - 2) if i1 >= 0 else len(raw)
    return header, raw[data_begin:data_end]


def key_value(header, key, cast="str", exact=False):
    """Look up a header value, mirroring MATLAB ``key_value.m``.

    Parameters
    ----------
    header
        Parsed header (list of (key, value) pairs).
    key
        Key to look up.
    cast
        ``"str"`` returns the raw string; ``"num"`` parses the first
        whitespace-delimited token as a float (e.g. ``"50.0 hz" -> 50.0``,
        ``"365" -> 365.0``).
    exact
        If ``True`` require an exact key match; otherwise match by prefix
        (the MATLAB default), returning the first match.

    Returns
    -------
    The value, or ``None`` if the key is absent.
    """
    match = None
    for k, v in header:
        if (k == key) if exact else k.startswith(key):
            match = v
            break
    if match is None:
        return None
    if cast == "num":
        token = match.strip().split(" ")[0].strip() if match.strip() else ""
        return float(token) if token else None
    return match


def big_endian_view(u8, dtype):
    """Interpret a contiguous uint8 array as a big-endian ``dtype`` array.

    The final axis is collapsed: ``(..., n_bytes)`` becomes
    ``(..., n_bytes / itemsize)``.
    """
    contiguous = np.ascontiguousarray(u8, dtype=np.uint8)
    return contiguous.view(np.dtype(dtype).newbyteorder(">"))


def little_endian_view(u8, dtype):
    """Interpret a contiguous uint8 array as a little-endian ``dtype`` array.

    Unlike every other DACQ file (big-endian), ``.egf`` is written
    little-endian (``getegf.m`` opens it with MATLAB's ``'ieee-le'``).
    """
    contiguous = np.ascontiguousarray(u8, dtype=np.uint8)
    return contiguous.view(np.dtype(dtype).newbyteorder("<"))
