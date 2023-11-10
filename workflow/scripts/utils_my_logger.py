#! /usr/bin/env/python
# coding: utf-8

## -----
# In-house decorator functions
# Author: G. Molano, LA (gonmola@hotmail.es)
# Last modified:
## -----
import logging

class My_logger:
    def __init__(self, log_filename, logger_name=None, level=None):
        self.logger_name = logger_name

        # create logger
        if logger_name:
            logger = logging.getLogger(logger_name)
        else:
            logger = logging.getLogger(__name__)
        
        if level:
            if level == "INFO":     
                logger.setLevel(logging.INFO)
            elif level == "WARNING":
                logger.setLevel(logging.WARNING)
            else:
                logger.setLevel(logging.DEBUG)
        else:
            logger.setLevel(logging.DEBUG)
        # configure the handler and formatter as needed
        handler = logging.FileHandler(log_filename, mode='a')
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        # add formatter to handler
        handler.setFormatter(formatter)
        # add handler to logger
        logger.addHandler(handler)
    
    def get_logger(self):
        logger = logging.getLogger(self.logger_name)
        return logger
