#!/usr/bin/env python3
"""
R Integration Compatibility Layer for MetExtract II

R/rpy2 has been removed. This module is kept for backward compatibility
but all functions now indicate that R is not available.
"""

import logging

# Global flag – R is no longer used
R_AVAILABLE = False
R_ERROR_MESSAGE = "R/rpy2 has been removed. Peak picking uses native Python."


def check_r_availability():
    """R is no longer required – always returns False."""
    return False


def get_r_interface():
    """R is no longer available – always returns None."""
    return None


def is_r_available():
    """R is no longer available – always returns False."""
    return False


def get_r_error_message():
    """Return message explaining R removal."""
    return R_ERROR_MESSAGE


def require_r(func):
    """Decorator – always raises since R was removed."""

    def wrapper(*args, **kwargs):
        raise RuntimeError(f"R has been removed from MetExtract II. Function {func.__name__} is no longer available.")

    return wrapper


def optional_r(fallback_func=None):
    """Decorator – always uses fallback since R was removed."""

    def decorator(func):
        def wrapper(*args, **kwargs):
            if fallback_func:
                logging.info(f"R removed, using fallback for {func.__name__}")
                return fallback_func(*args, **kwargs)
            else:
                logging.warning(f"R removed and no fallback for {func.__name__}")
                return None

        return wrapper

    return decorator
