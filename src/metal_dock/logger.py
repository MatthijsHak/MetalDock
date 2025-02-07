import sys
import logging
from pathlib import Path

class MetalDockLogger:
    _instance = None
    _initialized = False

    def __new__(cls):
        if cls._instance is None:
            cls._instance = super(MetalDockLogger, cls).__new__(cls)
        return cls._instance

    def __init__(self):
        if not MetalDockLogger._initialized:
            self.logger = logging.getLogger('MetalDock')
            self.logger.setLevel(logging.INFO)
            
            # Create formatters
            console_formatter = logging.Formatter('%(message)s')
            file_formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
            
            # Console handler
            console_handler = logging.StreamHandler(sys.stdout)
            console_handler.setFormatter(console_formatter)
            self.logger.addHandler(console_handler)
            
            # Store handlers for later use
            self.handlers = {'console': console_handler}
            
            MetalDockLogger._initialized = True

    def setup_file_logger(self, log_dir: Path):
        """
        Setup file logging based on the method.

        Args:
            log_dir (Path): The directory to store the log file.
        """
        # Define log file path
        log_file_path = log_dir / f'MetalDock.log'

        # Remove existing log file if it exists
        if log_file_path.exists():
            log_file_path.unlink()

        # Create file handler
        file_handler = logging.FileHandler(log_file_path)
        file_formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        file_handler.setFormatter(file_formatter)
        
        self.logger.addHandler(file_handler)
        self.handlers['file'] = file_handler

    def set_level(self, level: int):
        """
        Set logging level.

        Args:
            level (int): The logging level.
        """
        self.logger.setLevel(level)

    def info(self, msg: str):
        """
        Log info message.

        Args:
            msg (str): The message to log.
        """
        self.logger.info(msg)

    def debug(self, msg: str):
        """
        Log debug message.

        Args:
            msg (str): The message to log.
        """
        self.logger.debug(msg)

    def warning(self, msg: str):
        """
        Log warning message.

        Args:
            msg (str): The message to log.
        """
        self.logger.warning(msg)

    def error(self, msg: str):
        """
        Log error message.

        Args:
            msg (str): The message to log.
        """
        self.logger.error(msg)

    def critical(self, msg: str):
        """
        Log critical message.

        Args:
            msg (str): The message to log.
        """
        self.logger.critical(msg)
