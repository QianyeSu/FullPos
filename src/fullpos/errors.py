class FullposError(Exception):
    """Base exception for fullpos."""


class FullposNotImplementedError(FullposError, NotImplementedError):
    """Raised for planned behavior that is not implemented yet."""


class FullposBackendError(FullposError, RuntimeError):
    """Raised when the required native FULLPOS/ECTRANS backend fails."""
