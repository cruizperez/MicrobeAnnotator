import json

import requests

from microbeannotator.shared.constants import MODULE_JSON_URL

module_text = requests.get(MODULE_JSON_URL)
module_json = json.loads(module_text.text)
print(module_json)
