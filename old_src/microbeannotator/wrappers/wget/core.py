# from dataclasses import dataclass

import requests
from bs4 import BeautifulSoup
import re
# from microbeannotator.wrappers.common import (ExecutableWrapper,
#                                               WrapperAlgorithmBase,
#                                               get_cmd_with_flags)

# class WgetAlgorithm(WrapperAlgorithmBase):
#     input_url: str


api_url = "http://rest.kegg.jp/list/M00009"
api_url = "https://www.kegg.jp/module/M00001"
response = requests.get(api_url)

# Parse the HTML using BeautifulSoup
if response.status_code == 200:
    html = response.content  # Get the HTML content from the response
    soup = BeautifulSoup(html, 'html.parser')

# Find the map with name "module1"
map_element = soup.find('map', {'name': 'module1'})
print(map_element)
# # Find all area tags within the map
area_tags = map_element.find_all('area')

module_structure = {}
for area in area_tags:
    kegg_id = area['title'].split()
    if len(kegg_id) > 1:
        protein = kegg_id[1]
    kegg_id = kegg_id[0]
    protein = protein.replace("(", "").replace(")", "")
    coords = area['coords'].split(',')
    up_loc = coords[1]
    step_id = area['id']
    level, step_num = step_id.split("_")
    step_num = int(step_num)
    level = int(level.replace('node', ''))
    if level not in module_structure:
        module_structure[level] = {up_loc: [[kegg_id, protein, coords]]}
    elif up_loc not in module_structure[level]:
        module_structure[level][up_loc] = [[kegg_id, protein, coords]]
    else:
        module_structure[level][up_loc].append([kegg_id, protein, coords])

print(module_structure)
