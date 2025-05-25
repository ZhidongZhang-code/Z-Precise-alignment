#!/usr/bin/env python3

import logging

def setup_logging(log_file_path=None, log_level=logging.INFO):
    """
    配置日志系统。

    :param log_file_path: 日志文件的路径。如果未指定，则日志输出到控制台。
    :param log_level: 日志级别，默认为 logging.INFO。
    :return: 配置好的 logger 对象。
    """
    logger = logging.getLogger('MyScriptLogger')
    logger.setLevel(log_level)

    # 移除已存在的所有处理器，避免重复日志
    if logger.handlers:
        for handler in logger.handlers:
            logger.removeHandler(handler)

    # 根据是否指定了日志文件路径来决定日志的输出方式
    if log_file_path:
        file_handler = logging.FileHandler(log_file_path)
    else:
        file_handler = logging.StreamHandler()

    file_handler.setLevel(log_level)
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    return logger
