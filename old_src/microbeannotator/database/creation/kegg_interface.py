from logging import Logger
from typing import Optional

from microbeannotator.shared.logger import MicrobeAnnotatorLogger
from microbeannotator.shared.constants import MODULE_JSON_URL

import json

import requests

module_text = requests.get(MODULE_JSON_URL)
module_json = json.loads(module_text.text)
print(module_json)


class KeggInterfacer:
    def __init__(self, logger: Optional[Logger] = None) -> None:
        self._logger = logger if logger else MicrobeAnnotatorLogger(self.__class__.__name__)
    
    def _get_module_list(self, )


