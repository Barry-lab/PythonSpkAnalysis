"""Front-end loaders: format-specific readers that emit a common Recording."""

from .axona import load_dacq

__all__ = ["load_dacq"]
