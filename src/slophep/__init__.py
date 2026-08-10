import slophep.Core as Core
import slophep.FormFactors as FormFactors
import slophep.Observables as Observables
import slophep.Tools as Tools


# Setup logging
import logging

handler = logging.StreamHandler()
handler.setLevel(logging.DEBUG)
formatter = logging.Formatter(
    "%(levelname)s:%(name)s: %(message)s"
)
handler.setFormatter(formatter)


logger = logging.getLogger()
# logger.addHandler(handler)
logger.setLevel(logging.INFO)
