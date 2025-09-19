#!/usr/bin/env python3
"""
R Integration Compatibility Layer for MetExtract II
This module provides graceful fallbacks when R/rpy2 is not available
"""

import os
import sys
import logging

# Global flag to track R availability
R_AVAILABLE = False
R_ERROR_MESSAGE = ""

def check_r_availability():
    """Check if R and rpy2 are available"""
    global R_AVAILABLE, R_ERROR_MESSAGE
    
    try:
        import rpy2.robjects as ro
        r = ro.r
        v = r("R.Version()$version.string")
        R_AVAILABLE = True
        logging.info(f"R is available: {v}")
        return True
    except ImportError as e:
        R_ERROR_MESSAGE = f"rpy2 module not installed: {e}"
        logging.warning(R_ERROR_MESSAGE)
    except Exception as e:
        R_ERROR_MESSAGE = f"R runtime error: {e}"
        logging.warning(R_ERROR_MESSAGE)
    
    R_AVAILABLE = False
    return False

def get_r_interface():
    """Get R interface if available, None otherwise"""
    if not R_AVAILABLE:
        return None
    
    try:
        import rpy2.robjects as ro
        return ro.r
    except:
        return None

def is_r_available():
    """Check if R is available"""
    return R_AVAILABLE

def get_r_error_message():
    """Get the error message if R is not available"""
    return R_ERROR_MESSAGE

def require_r(func):
    """Decorator to mark functions that require R"""
    def wrapper(*args, **kwargs):
        if not R_AVAILABLE:
            raise RuntimeError(f"R is required for this function but is not available: {R_ERROR_MESSAGE}")
        return func(*args, **kwargs)
    return wrapper

def optional_r(fallback_func=None):
    """Decorator for functions that can work without R"""
    def decorator(func):
        def wrapper(*args, **kwargs):
            if R_AVAILABLE:
                return func(*args, **kwargs)
            elif fallback_func:
                logging.info(f"R not available, using fallback for {func.__name__}")
                return fallback_func(*args, **kwargs)
            else:
                logging.warning(f"R not available and no fallback for {func.__name__}")
                return None
        return wrapper
    return decorator

# Check R availability on module import
check_r_availability()