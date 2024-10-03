# utils.py

import logging
import os

def setup_logging(step_name, out_dir):
    """
    Set up logging for a specific step and create log files inside the /log folder.
    
    :param step_name: Name of the step (used to name the log file).
    :param out_dir: Output directory where logs will be stored in a /log folder.
    :return: A configured logger instance.
    """
    log_dir = os.path.join(out_dir, "log")
    if not os.path.exists(log_dir):
        os.makedirs(log_dir)
    
    log_file = os.path.join(log_dir, f"{step_name}.log")
    
    logging.basicConfig(
        filename=log_file,
        filemode='a',  # Append mode
        format='%(asctime)s - %(levelname)s - %(message)s',
        level=logging.INFO
    )
    logger = logging.getLogger(step_name)
    return logger
