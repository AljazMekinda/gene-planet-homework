from datetime import datetime
import os
# from utils.ftplib_ie import FTPClient
import yaml
import json
import logging
import logging.config
import io

logger = logging.getLogger("utils")


def load_config(file_path: str, type="config", key=None):
    with open(file_path, "rb") as f:
        file_contents = f.read()
    file_type = file_path.split('.')[-1]
    if file_type == 'json':
        config = json.loads(file_contents)
    elif file_type == 'yaml':
        config = yaml.safe_load(file_contents)
    else:
        raise ValueError(f"Unrecognized file type: {file_type}")
    if type == "config":
        if key is None:
            return config
        else:
            return config[key]


def setup_logging(logging_config_path, filename_prefix=None, logger_name=None,
                  logs_timestamp=None):
    """
    Setup logging configuration. Creates a folder for logs if it does not exist.
    If filename_prefix is not None, it will create a folder with that name and
    save the logs there.
    :param logging_config_path:
    :param filename_prefix:
    :param logger_name:
    :return:
    """
    create_if_not_exists("logs")
    create_if_not_exists(os.path.join("logs", filename_prefix))
    if logs_timestamp is None:
        logs_timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    create_if_not_exists(os.path.join("logs", filename_prefix, logs_timestamp))

    logging_config = load_config(logging_config_path)

    if filename_prefix is not None:
        for handler, handler_config in logging_config["handlers"].items():
            if "filename" in handler_config.keys():
                handler_config[
                    "filename"] = os.path.join("logs", filename_prefix,
                                               logs_timestamp,
                                               f"{handler_config['filename']}")
                logging_config["handlers"].update({handler: handler_config})

    logging.config.dictConfig(logging_config)
    if logger_name is not None:
        return logging.getLogger(logger_name)


def get_saving_log(logger):
    try:
        with open(os.path.join("logs", "gui_saving.log"), 'r') as file:
            log_contents = file.read()
    except FileNotFoundError:
        logger.info("No saving log found")
        log_contents = "No saving log found"
    return log_contents


def create_if_not_exists(dir):
    """
    Creates the folder on the local disk if it does not exist.

    :param dir: Folder path
    :type dir: str
    """
    if not os.path.exists(os.path.normpath(dir)):
        os.makedirs(os.path.normpath(dir))


