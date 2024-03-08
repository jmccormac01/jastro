"""
Custom logging for jastro
"""
import sys
import logging

def custom_logger(logger_name, logger_type, level=logging.DEBUG):
    """
    Method to return a custom logger with the given name and level

    Parameters
    ----------
    logger_name : string
        name of logfile for this object
    level : logging.level
        which serverity level to display beyond

    Returns
    -------
    logger : logging.logger
        used to log info to correct files during reduction

    Raises
    ------
    None
    """
    logger = logging.getLogger(logger_name)
    logger.setLevel(level)
    format_string = ("%(asctime)s — %(levelname)s — %(funcName)s:"
                     "%(lineno)d — %(message)s")
    log_format = logging.Formatter(format_string)
    # Creating and adding the console handler
    if logger_type in ("stdout", "both"):
        console_handler = logging.StreamHandler(sys.stdout)
        console_handler.setFormatter(log_format)
        logger.addHandler(console_handler)
    # Creating and adding the file handler
    if logger_type in ("file", "both"):
        file_handler = logging.FileHandler(logger_name, mode='a')
        file_handler.setFormatter(log_format)
        logger.addHandler(file_handler)
    return logger
