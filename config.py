import logging
import time


class LoggingConfig:
    def initialize(self):
        # Force logging timestamps to UTC
        logging.Formatter.converter = time.gmtime

        LOG_FORMAT = (
            "%(asctime)s [%(filename)s:%(funcName)s] %(levelname)s - %(message)s"
        )

        logging.basicConfig(
            level=logging.INFO,
            format=LOG_FORMAT,
        )
