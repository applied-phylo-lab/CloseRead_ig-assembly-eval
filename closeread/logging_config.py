# logging_config.py
from __future__ import annotations
import logging
from pathlib import Path
import os
from datetime import datetime

FORMAT = "%(asctime)s | %(levelname)-8s | %(name)s | %(message)s"
DATEFMT = "%Y-%m-%d %H:%M:%S"

def setup_global_logging(log_dir: Path, level: str = "INFO", console: bool = True) -> None:
    log_dir = Path(log_dir)
    log_dir.mkdir(parents=True, exist_ok=True)
    suffix = datetime.now().strftime("%Y%m%d_%H%M%S")
    logfile = os.path.join(log_dir, f"pipeline_{suffix}.log")

    root = logging.getLogger()
    # Clear existing handlers (important for re-runs)
    for h in list(root.handlers):
        root.removeHandler(h)

    root.setLevel(getattr(logging, level.upper()))

    file_handler = logging.FileHandler(logfile, mode="w")  # truncate each run
    file_handler.setFormatter(logging.Formatter(FORMAT, DATEFMT))
    root.addHandler(file_handler)

    if console:
        sh = logging.StreamHandler()
        sh.setFormatter(logging.Formatter(FORMAT, DATEFMT))
        root.addHandler(sh)

    root.info("Logging initialized -> %s", logfile)
