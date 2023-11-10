#! /usr/bin/env/python
# coding: utf-8

## -----
# In-house decorator functions
# Author: G. Molano, LA (gonmola@hotmail.es)
# Last modified:
## -----

import functools
import time
from typing import Union
import logging

class MyLogger:
    def __init__(self):
        logging.basicConfig(level=logging.DEBUG, filename = "log.log", filemode = "a",
                            format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')

    def get_logger(self, name=None):
        return logging.getLogger(name)


def get_default_logger():
    return MyLogger().get_logger()

def get_logger(my_logger):
    """
    Handle three different types of logger:
    1. No logger is passed: decorators creates a logger.
    2. `My_logger` is passed: an instance of custome logger class is passed.
    3. `logging.logger` is passed: logger itsefs is passed.
    """
    # Logger options
    if my_logger is None:
        logger = get_default_logger()
    else:
        if isinstance(my_logger, MyLogger):
            logger = my_logger.get_logger(__name__)
        else:
            logger = my_logger
    return logger

def timer(_func=None, *, my_logger: Union[MyLogger, logging.Logger] = None):
    """Print the runtime of the decorated function

    """
    def decorator_log(func):
        @functools.wraps(func)
        def wrapper_timer(*args, **kwargs):
            # Logger options
            logger = get_logger(my_logger)

            start_time = time.perf_counter()    # 1
            value = func(*args, **kwargs)
            end_time = time.perf_counter()      # 2
            run_time = end_time - start_time    # 3
            logger.debug(f"Finished {func.__name__!r} in {run_time:.4f} secs")
            return value
        return wrapper_timer
    
    if _func is None:
        return decorator_log
    else:
        return decorator_log(_func)

def debuger(_func=None, *, my_logger: Union[MyLogger, logging.Logger] = None):
    """Print the function signature and return value
    """
    def decorator_log(func):
        @functools.wraps(func)
        def wrapper_debug(*args, **kwargs):
            # Logger options
            logger = get_logger(my_logger)

            args_repr = [repr(a) for a in args]                      # 1
            kwargs_repr = [f"{k}={v!r}" for k, v in kwargs.items()]  # 2
            signature = ", ".join(args_repr + kwargs_repr)           # 3

            logger.debug(f"function {func.__name__} called with args {signature}")
            try:
                value = func(*args, **kwargs)
                logger.debug(f"{func.__name__!r} returned {value!r}")           # 4
                return value
            except Exception as e:
                logger.exception(f"Exception raised in {func.__name__}. exception: {str(e)}")
                raise e
        return wrapper_debug
   
    if _func is None:
        return decorator_log
    else:
        return decorator_log(_func)
