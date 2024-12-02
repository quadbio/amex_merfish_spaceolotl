# From: https://github.com/scverse/spatialdata/blob/main/src/spatialdata/_logging.py

import logging
from rich.console import Console
from rich.logging import RichHandler
import os

def _setup_logger() -> logging.Logger:
    
    """Setup the logger with rich formatting."""
    
    logger = logging.getLogger(__name__)

    # This is necessary because the sub-scripts import _setup_logger and the
    # quant_unsupervised.py imports _setup_logger as well as the sub-scripts,
    # which adds multiple handlers
    if logger.hasHandlers():
        return logger
        
    level = os.environ.get("LOGLEVEL", logging.INFO)
    logger.setLevel(level=level)
    console = Console(force_terminal=True)
    
    if console.is_jupyter:
        console.is_jupyter = False
    
    ch = RichHandler(show_path=False, console=console, show_time=logger.level == logging.DEBUG)
    logger.addHandler(ch)

    # Prevents double outputs
    logger.propagate = False
    return logger