import logging

import lxml.etree as et


class SimComponent:
    """
    Super class for all simulation components. Attaches a logger to each component intances.
    """

    def __init__(self):
        self.logger = logging.getLogger(self.__module__ + "." + self.__class__.__name__)
        # debug("creating logger for "+str(self.__class__.__name__),self.logger)

    def to_xml(self) -> et.Element:
        pass

    def from_xml(self, e: et.Element) -> None:
        pass
