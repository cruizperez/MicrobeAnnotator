from selenium import webdriver
from bs4 import BeautifulSoup
import pandas as pd
import re
import pickle
import ast


""" Script to download and parse KEGG information and store it in data """

def download_kegg_modules(module_name_file, chrome_driver):
    module_ids =[]
    module_names = {}
    module_components_raw = {}
    # Parse module names
    with open(module_name_file) as module_input:
        for line in module_input:
            line = line.strip().split("\t")
            module_ids.append(line[0])
            module_names[line[0]] = line[1]
    driver = webdriver.Chrome(chrome_driver)
    # Access KEGG and download module information
    for identifier in module_ids:
        url = "https://www.kegg.jp/kegg-bin/show_module?" + identifier
        driver.get(url)
        content = driver.page_source
        soup = BeautifulSoup(content, features="lxml")
        # title = soup.findAll('b')[0].get_text()
        formula = soup.findAll('td')[7].get_text().strip()
        module_components_raw[identifier] = formula
    driver.close()
    return module_components_raw


def parse_regular_module_dictionary(bifurcating_list_file, structural_list_file, module_components_raw):
    bifurcating_list = []
    structural_list = []
    # Populate bifurcating and structural lists
    with open(bifurcating_list_file, 'r') as bif_list:
        for line in bif_list:
            bifurcating_list.append(line.strip())
    with open(structural_list_file, 'r') as bif_list:
        for line in bif_list:
            structural_list.append(line.strip())
    # Parse raw module information
    module_steps_parsed = {}
    for key, values in module_components_raw.items():
        values = values.replace(" --", "")
        values = values.replace("-- ", "")
        if key in bifurcating_list or key in structural_list:
            continue
        else:
            module = []
            parenthesis_count = 0
            for character in values:
                if character == "(":
                    parenthesis_count += 1
                    module.append(character)
                elif character == " ":
                    if parenthesis_count == 0:
                        module.append(character)
                    else:
                        module.append("_")
                elif character == ")":
                    parenthesis_count -= 1
                    module.append(character)
                else:
                    module.append(character)
            steps = ''.join(module).split()
            module_steps_parsed[key] = steps
    # Remove modules depending on other modules
    temporal_dictionary = module_steps_parsed.copy()
    for key, values in temporal_dictionary.items():
        for value in values:
            if re.search(r'M[0-9]{5}', value) is not None:
                del module_steps_parsed[key]
                break
    return module_steps_parsed


def create_final_regular_dictionary(module_steps_parsed, module_components_raw, outfile):
    final_regular_dict = {}
    # Parse module steps and export them into a text file
    with open(outfile, 'w') as output:
        for key, value in module_steps_parsed.items():
            output.write("{}\n".format(key))
            output.write("{}\n".format(module_components_raw[key]))
            output.write("{}\n".format(value))
            output.write("{}\n".format("=="))
            final_regular_dict[key] = {}
            step_number = 0
            for step in value:
                step_number += 1
                count = 0
                options = 0
                temp_string = ""
                for char in step:
                    if char == "(":
                        count += 1
                        options += 1
                        if len(temp_string) > 1 and temp_string[-1] == "-":
                            temp_string += "%"
                    elif char == ")":
                        count -= 1
                        if count >= 1:
                            temp_string += char
                        else:
                            continue
                    elif char == ",":
                        if count >= 2:
                            temp_string += char
                        else:
                            temp_string += " "
                    else:
                        temp_string += char
                if options >= 2:
                    temp_string = temp_string.replace(")_", "_")
                    if re.search('%.*\)', temp_string) is None:
                        temp_string = temp_string.replace(")", "")
                    temp_string = "".join(temp_string.rsplit("__", 1))
                    temp_string = temp_string.split()
                if isinstance(temp_string, str):
                    temp_string = temp_string.split()
                temp_string = sorted(temp_string, key=len)
                final_regular_dict[key][step_number] = temp_string
                output.write("{}\n".format(temp_string))
            output.write("{}\n".format("++++++++++++++++++"))
    return final_regular_dict


def export_module_dictionary(dictionary, location):
    pickle_out = open(location,"wb")
    pickle.dump(dictionary, pickle_out)
    pickle_out.close()



def transform_module_dictionaries(bifurcating_data, structural_data, output_bifur, output_struct):
    bifurcating_dictionary = ast.literal_eval(open(bifurcating_data).read())
    export_module_dictionary(bifurcating_dictionary, output_bifur)
    structural_dictionary = ast.literal_eval(open(structural_data).read())
    export_module_dictionary(structural_dictionary, output_struct)


# Execute parsing functions

module_components_raw = download_kegg_modules("00.Module_Names.txt", 'chromedriver')
module_steps_parsed = parse_regular_module_dictionary("01.Bifurcating_List.txt", 
                                "02.Structural_List.txt", module_components_raw)
final_regular_dict = create_final_regular_dictionary(module_steps_parsed, module_components_raw, "05.Modules_Parsed.txt")


export_module_dictionary(final_regular_dict, "../01.KEGG_Regular_Module_Information.pickle")
transform_module_dictionaries("03.Bifurcating_Modules.dict", 
                              "04.Structural_Modules.dict", 
                              "../02.KEGG_Bifurcating_Module_Information.pickle", 
                              "../03.KEGG_Structural_Module_Information.pickle")