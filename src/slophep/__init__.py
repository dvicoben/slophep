
import logging

handler = logging.StreamHandler()
handler.setLevel(logging.DEBUG)
formatter = logging.Formatter(
    "%(levelname)s:%(name)s: %(message)s"
)
handler.setFormatter(formatter)


logger = logging.getLogger()
logger.addHandler(handler)
logger.setLevel(logging.INFO)
