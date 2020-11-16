#!/usr/bin/env python

"""
########################################################################
# Author:       Carlos A. Ruiz Perez
# Email:        cruizperez3@gatech.edu
# Intitution:   Georgia Institute of Technology
# Version:      1.0.0
# Date:         Nov 13, 2020

# Description: Maps protein KO information with their respective modules
# and calculates the completeness percentage of each module present.
########################################################################
"""

################################################################################
"""---0.0 Import Modules---"""
import pandas as pd
import re
import pickle
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
from collections import OrderedDict
from sys import exit
import matplotlib

################################################################################
"""---1.0 Define Functions---"""

def ko_match(string):
    """ Looks if string has the form K[0-9]{5}
    
    Arguments:
        string {string} -- String to test
    
    Returns:
        [bool] -- String has form or not
    """
    if re.search(r'^K[0-9]{5}$', string) is not None:
        return 1
    else:
        return 0

def split_compound(string, comp_type):
    """[summary]
    
    Arguments:
        string {[type]} -- [description]
        comp_type {[type]} -- [description]
    
    Returns:
        [type] -- [description]
    """
    if comp_type == 'compound':
        return re.split('[-+]', string)
    elif comp_type == 'and_comp':
        return string.split('_')
    elif comp_type == 'or_comp':
        return string.split(',')

def process_compounds(string, comp_type, genome_annotation):
    string = split_compound(string, comp_type)
    proteins_required = len(string)
    proteins_present = 0
    for option in string:
        if option == '':
            proteins_required -= 1
        elif '+' in option or '-' in option:
            compound = split_compound(option, 'compound')
            proteins_in_compound = len(compound)
            present_in_compound = 0
            for sub_option in compound:
                if ko_match(sub_option) > 0 and sub_option in genome_annotation:
                    present_in_compound += 1
            proteins_present += present_in_compound/proteins_in_compound
        else:
            if ko_match(option) > 0 and option in genome_annotation:
                proteins_present += 1
    return proteins_present, proteins_required

def import_module_dict(dictionary_location):
    pickle_in = open(dictionary_location,"rb")
    dictionary = pickle.load(pickle_in)
    return dictionary

def module_information_importer(input_files):
    print("Importing data... ", end="")
    script_path = Path(__file__)
    script_dir = script_path.parent
    data_folder = Path(script_dir).parent / "data"
    # Import all modules from dictionaries
    regular_modules = import_module_dict(data_folder / "01.KEGG_Regular_Module_Information.pickle")
    bifurcation_modules = import_module_dict(data_folder / "02.KEGG_Bifurcating_Module_Information.pickle")
    structural_modules = import_module_dict(data_folder / "03.KEGG_Structural_Module_Information.pickle")
    # Get genome and module ids
    module_information = {}
    module_ids = []
    genome_names = []
    group_names = []
    # Fill genome names list
    for genome_file in input_files:
        genome_names.append(str(Path(genome_file).name))
    # Get module correspondence
    with open(data_folder / "01.KEGG_DB/00.Module_Names.txt") as correspondence:
        for line in correspondence:
            line = line.strip().split("\t")
            module_information[line[0]] = [line[1], line[2], line[3]]
            # module_id_name[line[0]] = line[1]
            module_ids.append(line[0])
            # module_groups[line[0]] = line[2]
            if line[2] not in group_names:
                group_names.append(line[2])
            # module_colors[line[1]] = [line[2],line[3]]
    # Create module pathway group matrix
    module_group_matrix = pd.DataFrame(0, index=group_names, columns=genome_names)
    # Create module matrix
    metabolism_matrix = pd.DataFrame(0, index=module_ids, columns=genome_names)
    print("Done!")
    return regular_modules, bifurcation_modules, structural_modules, \
           module_information, metabolism_matrix, module_group_matrix

def annotation_file_parser(genome_annotation_file):
    genome_annotation = []
    with open(genome_annotation_file) as annot_input:
        for line in annot_input:
            genome_annotation.append(line.strip())
    return genome_annotation

def regular_module_mapper(module_dictionary, genome_annotation_file):
    genome_annotation = annotation_file_parser(genome_annotation_file)
    regular_module_completenes = []
    for module, final_steps in module_dictionary.items():
        complete_steps = 0
        total_module_steps = len(final_steps)
        for proteins in final_steps.values():
            score_for_step = 0
            steps_per_option = 100
            match = False
            for option in proteins:
                match = False
                # Search for single gene option
                if ko_match(option) > 0:
                    if option in genome_annotation:
                        score_for_step = 1
                        match = True
                    if match == True:
                        steps_per_option = 1
                    elif 1 < steps_per_option:
                        steps_per_option = 1
                elif '%' in option:
                    option = option.replace(')', '')
                    option = option.split('-%')
                    mandatory = split_compound(option[0], 'compound')
                    proteins_required = len(mandatory) + 1
                    proteins_present = 0
                    for prot in mandatory:
                        if ko_match(prot) > 0 and prot in genome_annotation:
                            proteins_present += 1
                    for prot in option[1].split(','):
                        if ko_match(prot) > 0 and prot in genome_annotation:
                            proteins_present += 1
                            break
                    if proteins_present/proteins_required > score_for_step:
                        score_for_step = proteins_present/proteins_required
                        match = True
                    if match == True:
                        steps_per_option = 1
                    elif 1 < steps_per_option:
                        steps_per_option = 1
                elif '_' in option and ',' in option:
                    option = sorted(split_compound(option, 'and_comp'), key=len)
                    highest_score = 0
                    for element in option:
                        if ko_match(element) > 0 and element in genome_annotation:
                            highest_score += 1
                        elif ',' in element:
                            element = sorted(element.split(","), key=len)
                            for sub_element in element:
                                if ko_match(sub_element) > 0 and sub_element in genome_annotation:
                                    highest_score += 1
                                    break
                                elif '+' in sub_element:
                                    proteins_present, proteins_required = process_compounds(sub_element, 
                                    'compound', genome_annotation)
                                    highest_score += proteins_present/proteins_required
                    if highest_score > score_for_step:
                        score_for_step = highest_score
                        match = True
                    if match == True:
                        steps_per_option = len(option)
                    elif len(option) < steps_per_option:
                        steps_per_option = len(option)
                elif len(option.split(",")) > 1:
                    match = False
                    option = sorted(option.split(","), key=len)
                    highest_score = 0
                    steps_to_add = 0
                    for sub_option in option:
                        if ko_match(sub_option) > 0 and sub_option in genome_annotation:
                            if 1 > highest_score:    
                                highest_score = 1
                                steps_to_add = 1
                        elif '+' in sub_option or '-' in sub_option:
                            proteins_present, proteins_required = process_compounds(sub_option, 'compound', genome_annotation)
                            if proteins_present/proteins_required > highest_score:
                                highest_score = proteins_present/proteins_required
                                steps_to_add = 1
                    if highest_score > score_for_step:
                        score_for_step = highest_score
                        match = True
                    if match == True:
                        steps_per_option = steps_to_add
                    elif steps_to_add< steps_per_option:
                        steps_per_option = steps_to_add
                elif '_' in option:
                    match = False
                    proteins_present, proteins_required = process_compounds(option, 'and_comp', genome_annotation)
                    if proteins_present/proteins_required > score_for_step:
                        score_for_step = proteins_present
                        match = True
                    if match == True:
                        steps_per_option = proteins_required
                    elif proteins_required < steps_per_option:
                        steps_per_option = proteins_required
                elif "+" in option or "-" in option:
                    match = False
                    highest_score = 0
                    proteins_present, proteins_required = process_compounds(option, 'compound', genome_annotation)
                    if proteins_present/proteins_required > score_for_step:
                        score_for_step = proteins_present/proteins_required
                        match = True
                    if match == True:
                        steps_per_option = 1
                    elif proteins_required < steps_per_option:
                        steps_per_option = 1
                else:
                    print("Unreccognized module {}. Check your database.".format(option))
            complete_steps += score_for_step
            if steps_per_option > 50:
                steps_per_option = 1
            if steps_per_option > 1:
                total_module_steps += steps_per_option - 1
        regular_module_completenes.append((module, round((complete_steps/total_module_steps)*100, 2)))
    return regular_module_completenes

def bifurcating_module_mapper(module_dictionary, genome_annotation_file):
    genome_annotation = annotation_file_parser(genome_annotation_file)
    bifurcating_module_completenes = []
    for module, versions in module_dictionary.items():
        module_highest = 0
        for version, total_steps in versions.items():
            completed_steps = 0
            total_version_steps = len(total_steps)
            for proteins in total_steps.values():
                score_for_step = 0
                steps_per_option = 100
                match = False
                if isinstance(proteins, (list)):
                    protein = sorted(proteins, key=len)
                    for option in protein:
                        match = False
                        if ko_match(option) > 0:
                            if option in genome_annotation:
                                score_for_step = 1
                                match = True
                            if match == True:
                                steps_per_option = 1
                            elif 1 < steps_per_option:
                                steps_per_option = 1
                        elif '_' in option and ',' in option:
                            option = sorted(split_compound(option, 'and_comp'), key=len)
                            highest_score = 0
                            for element in option:
                                if ko_match(element) > 0 and element in genome_annotation:
                                    highest_score += 1
                                elif ',' in element:
                                    element = sorted(element.split(","), key=len)
                                    score_sub_option = 0
                                    for sub_element in element:
                                        if ko_match(sub_element) > 0 and sub_element in genome_annotation:
                                            highest_score += 1
                                            break
                                        elif '+' in sub_element:
                                            proteins_present, proteins_required = process_compounds(sub_element, 
                                            'compound', genome_annotation)
                                            highest_score += proteins_present/proteins_required
                            if highest_score > score_for_step:
                                score_for_step = highest_score
                                match = True
                            if match == True:
                                steps_per_option = len(option)
                        elif '_' in option:
                            match = False
                            proteins_present, proteins_required = process_compounds(option, 'and_comp', genome_annotation)
                            if proteins_present/proteins_required > score_for_step:
                                score_for_step = proteins_present
                                match = True
                            if match == True:
                                steps_per_option = proteins_required
                            elif proteins_required < steps_per_option:
                                steps_per_option = proteins_required
                        elif '+' in option or '-' in option:
                            match = False
                            highest_score = 0
                            proteins_present, proteins_required = process_compounds(option, 'compound', genome_annotation)
                            if proteins_present/proteins_required > score_for_step:
                                score_for_step = proteins_present/proteins_required
                                match = True
                            if match == True:
                                steps_per_option = 1
                            elif proteins_required < steps_per_option:
                                steps_per_option = 1
                elif ko_match(proteins) > 0:
                    match = False
                    if proteins in genome_annotation:
                        score_for_step = 1
                        match = True
                    if match == True:
                        steps_per_option = 1
                    elif 1 < steps_per_option:
                        steps_per_option = 1
                elif ',' in proteins:
                    options = split_compound(proteins, 'or_comp')
                    for option in options:
                        if ko_match(option) > 0 and option in genome_annotation:
                            score_for_step = 1
                            match = True
                        if match == True:
                            steps_per_option = 1
                        elif 1 < steps_per_option:
                            steps_per_option = 1
                elif '+' in proteins or '-' in proteins:
                    match = False
                    highest_score = 0
                    proteins_present, proteins_required = process_compounds(proteins, 'compound', genome_annotation)
                    if proteins_present/proteins_required > score_for_step:
                        score_for_step = proteins_present/proteins_required
                        match = True
                    if match == True:
                        steps_per_option = 1
                    elif proteins_required < steps_per_option:
                        steps_per_option = 1
                else:
                    print("Unreccognized module {}. Check your database.".format(proteins))
                completed_steps += score_for_step
                if steps_per_option > 50:
                    steps_per_option = 1
                if steps_per_option > 1:
                    total_version_steps += steps_per_option - 1
            if completed_steps/total_version_steps > module_highest:
                module_highest = completed_steps/total_version_steps
        bifurcating_module_completenes.append((module, round(module_highest*100, 2)))
    return bifurcating_module_completenes

def structural_module_mapper(module_dictionary, genome_annotation_file):
    genome_annotation = annotation_file_parser(genome_annotation_file)
    structural_module_completenes = []
    for module, components in module_dictionary.items():
        score_for_components = 0
        module_proteins_present = 0
        module_proteins_required = 0
        for proteins in components:
            if isinstance(proteins, (list)):
                highest_score = 0
                proteins_to_add = 0
                steps_to_add = 100
                for option in proteins:
                    if '_' in option and ',' in option:
                        option = sorted(split_compound(option, 'and_comp'), key=len)
                        proteins_present_option = 0
                        proteins_required_option = 0
                        for element in option:
                            if ko_match(element) > 0:
                                if element in genome_annotation:
                                    proteins_present_option += 1
                                    proteins_required_option += 1
                                else:
                                    proteins_required_option += 1
                            elif ',' in element:
                                element = sorted(element.split(","), key=len)
                                score_sub_element = 0
                                proteins_present_sub_element = 0
                                proteins_required_sub_element = 0
                                for sub_element in element:
                                    if ko_match(sub_element) > 0:
                                        if sub_element in genome_annotation:
                                            score_sub_element = 1
                                            proteins_present_sub_element += 1
                                            proteins_required_sub_element += 1
                                        elif 1 < proteins_required_sub_element and score_sub_element == 0:
                                            proteins_required_sub_element = 1
                                    elif '+' in sub_element or '-' in sub_element:
                                        proteins_present, proteins_required = process_compounds(sub_element, 
                                        'compound', genome_annotation)
                                        if proteins_present/proteins_required > score_sub_element:
                                            score_sub_element = proteins_present/proteins_required
                                            proteins_present_sub_element = proteins_present
                                            proteins_required_sub_element = proteins_required
                                        elif proteins_required < proteins_required_sub_element and score_sub_element == 0:
                                            proteins_required_sub_element = proteins_required
                                proteins_present_option += proteins_present_sub_element
                                proteins_required_option += proteins_required_sub_element
                        if proteins_present_option/proteins_required_option > highest_score:
                            highest_score = proteins_present_option/proteins_required_option
                            proteins_to_add = proteins_present_option
                            steps_to_add = proteins_required_option
                        elif proteins_required_option < steps_to_add and highest_score == 0:
                            steps_to_add = proteins_required_option
                    elif '+' in option or '-' in option:
                        proteins_present, proteins_required = process_compounds(option, 'compound', genome_annotation)
                        if proteins_present/proteins_required > highest_score:
                            highest_score = proteins_present/proteins_required
                            steps_to_add = proteins_required
                            proteins_to_add = proteins_present
                        elif proteins_required < steps_to_add and highest_score == 0:
                            steps_to_add = proteins_required
                module_proteins_present += proteins_to_add
                module_proteins_required += steps_to_add
            elif ko_match(proteins) > 0:
                if proteins in genome_annotation:
                        module_proteins_present += 1
                        module_proteins_required += 1
                else:
                    module_proteins_required += 1
            elif ',' in proteins:
                proteins = sorted(proteins.split(","), key=len)
                score = 0
                proteins_present_option = 0
                proteins_required_option = 100
                for element in proteins:
                    if ko_match(element) > 0:
                        if element in genome_annotation:
                            proteins_required_option = 1
                            proteins_present_option = 1
                            break
                        elif 1 < proteins_required_option and score == 0:
                            proteins_required_option = 1
                    elif '+' in element or '-' in element:
                        proteins_present, proteins_required = process_compounds(element, 'compound', genome_annotation)
                        if proteins_present/proteins_required > score:
                            score = proteins_present/proteins_required
                            proteins_present_option = proteins_present
                            proteins_required_option = proteins_required
                        elif proteins_required < proteins_required_option and score == 0:
                            proteins_required_option = proteins_required
                module_proteins_present += proteins_present_option
                module_proteins_required += proteins_required_option
            elif '+' in proteins or '-' in proteins:
                proteins_present, proteins_required = process_compounds(proteins, 'compound', genome_annotation)
                module_proteins_present += proteins_present
                module_proteins_required += proteins_required  
            else:
                print("Unreccognized module {}. Check your database.".format(proteins))
        score_for_components = module_proteins_present/module_proteins_required
        structural_module_completenes.append((module, round(score_for_components*100, 2)))
    return structural_module_completenes

def global_mapper(regular_modules, bifurcating_modules, structural_modules, annotation_files):
    print("Mapping entries... ", end="")
    full_metabolic_completeness = {}
    
    for file in annotation_files:
        genome_name = str(Path(file).name)
        regular_completeness = regular_module_mapper(regular_modules, file)
        bifurcating_completeness = bifurcating_module_mapper(bifurcating_modules, file)
        structural_completeness = structural_module_mapper(structural_modules, file)
        final_completeness = regular_completeness + bifurcating_completeness + structural_completeness
        full_metabolic_completeness[genome_name] = final_completeness
    print("Done!\n")
    return full_metabolic_completeness

def plot_function_barplots(module_colors, module_group_matrix, metabolism_matrix_dropped_relabel, prefix):
    matplotlib.rcParams['pdf.fonttype'] = 42
    print("Grouping by metabolism type and plotting... ")
    for module in list(metabolism_matrix_dropped_relabel.index):
        for genome in list(metabolism_matrix_dropped_relabel.columns):
            if metabolism_matrix_dropped_relabel.loc[module,genome] >= 80:
                module_group_matrix.loc[module_colors[module][0],genome] += 1
    module_group_matrix_transp = module_group_matrix.T
    emptyness = (module_group_matrix_transp == 0).all()
    if emptyness.all() == True:
        print("There are no modules above 80% completeness. No barplot will be generated.")
    else:
        module_group_matrix_transp = module_group_matrix_transp.loc[(module_group_matrix_transp >= 1).any(1)]
        color_dict = {}
        color_list = []
        for pair in module_colors.values():
            if pair[0] not in color_dict:
                color_dict[pair[0]] = pair[1]
        for pathway in list(module_group_matrix_transp.columns):
            color_list.append(color_dict[pathway])
        Figure, Axis = plt.subplots(1,1, figsize=(25,15))
        module_group_matrix_transp.plot.bar(ax=Axis, stacked=True, color=color_list, figsize=(25,15), legend=False)
        Axis.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), fontsize='medium', markerscale=0.3)
        Axis.set_xlabel("Genomes", fontsize=30)
        Axis.set_ylabel('Number of Modules (>=80% complete)', fontsize=30)
        Axis.tick_params(axis='x', labelsize=15)
        Axis.tick_params(axis='y', labelsize=15)
        Figure.suptitle('Metabolism Module Category per Genome', fontsize=40)
        Figure.subplots_adjust()
        Figure.savefig(prefix + "_barplot.pdf", bbox_inches="tight")
    print("Done!")

def create_output_files(metabolic_annotation, metabolism_matrix, module_information, cluster, prefix):
    matplotlib.rcParams['pdf.fonttype'] = 42
    print("Creating output_file matrix and heatmap... ")
    # Check clustering asked
    if cluster is None:
        row_cluster = False
        col_cluster = False
    elif cluster == 'rows':
        row_cluster = True
        col_cluster = False
    elif cluster == 'cols':
        row_cluster = False
        col_cluster = True
    elif cluster == 'both':
        row_cluster = True
        col_cluster = True
    else:
        exit("Wrong clustering option, select between 'both', 'cols', 'rows' or don't pass this option.")
    # Extract data from metabolism matrix and color dictionary
    for genome, modules in metabolic_annotation.items():
        # Fill metabolism matrix with completeness from the metabolic annotation dictionary
        for module in modules:
            metabolism_matrix.loc[module[0],genome] = module[1]
    # Get module id and its corresponding module name and the module and its corresponding color
    module_id_name = {}
    module_colors = {}
    for module, information in module_information.items():
        module_id_name[module] = information[0]
        module_colors[information[0]] = (information[1],information[2])
    # Get modules that are above 50% complete in at least one genome
    ylabel_text = 'Modules (at least 50% complete in at least one genome)'
    metabolism_matrix_retained = metabolism_matrix.loc[(metabolism_matrix >= 50).any(1)]
    table_empty = (metabolism_matrix_retained == 0).all()
    if table_empty.all() == True:
        metabolism_matrix_retained = metabolism_matrix.loc[(metabolism_matrix >= 0).any(1)]
        ylabel_text = 'Modules'
    colors_for_ticks = []
    colors_for_legend = {}
    metabolism_matrix_retained_relabel = metabolism_matrix_retained.rename(index=module_id_name)
    # Save original annotation table but add module names and pathway
    module_names = []
    module_pathways = []
    for module in metabolism_matrix.index:
        module_names.append(module_information[module][0])
        module_pathways.append(module_information[module][1])
    metabolism_matrix.insert(loc=0, column="name", value=module_names)
    metabolism_matrix.insert(loc=1, column="pathway group", value=module_pathways)
    metabolism_matrix.index.name = 'module'
    metabolism_matrix.to_csv(prefix + "_module_completeness.tab", sep="\t")
    
    # Plot heatmap with clustering
    # Check if clustering is possible with the number of genomes and modules
    table_dim = metabolism_matrix_retained_relabel.shape
    if row_cluster == True and table_dim[0] < 3:
        row_cluster = False
    if col_cluster == True and table_dim[1] < 3:
        col_cluster = False
    sns.set(font_scale=1.6)
    heatmap = sns.clustermap(metabolism_matrix_retained_relabel, xticklabels=True, yticklabels=True, 
                            figsize=(30,25), row_cluster=row_cluster, col_cluster=col_cluster, cbar_pos=(0.01, 1, 0.2, 0.04),
                            cbar_kws= {'orientation': 'horizontal', 'label': 'Module Completeness (%)'}, dendrogram_ratio=0.1)
    # cbar = heatmap.ax_cbar
# here set the labelsize by 20
    # cbar.set_yticklabels(labels=cbar.get_ylabel, fontdict={"fontsize":20})
    for module in heatmap.ax_heatmap.yaxis.get_ticklabels():
        colors_for_ticks.append(module_colors[module.get_text()][1])
        if module_colors[module.get_text()][0] not in colors_for_legend:
            colors_for_legend[module_colors[module.get_text()][0]] = module_colors[module.get_text()][1]
    [mod.set_color(col) for (col,mod) in zip(colors_for_ticks,heatmap.ax_heatmap.yaxis.get_ticklabels())]
    for pathway, color in colors_for_legend.items():
        heatmap.ax_heatmap.scatter(1,1, color=color, label=pathway, zorder=0, s=45)
    heatmap.ax_heatmap.set_xticklabels(heatmap.ax_heatmap.get_xmajorticklabels(), fontsize = 15)
    heatmap.ax_heatmap.set_yticklabels(heatmap.ax_heatmap.get_ymajorticklabels(), fontsize = 15)
    heatmap.ax_heatmap.set_xlabel("Genomes", fontsize=20)
    heatmap.ax_heatmap.set_ylabel(ylabel_text, fontsize=20)
    plt.figlegend(loc="upper right", ncol=2, fontsize=12)
    heatmap.savefig(prefix + "_heatmap.pdf")
    print("Done!")
    return metabolism_matrix_retained_relabel, module_colors

################################################################################
"""---2.0 Main Function---"""

def main():
    import argparse, sys
    # Setup parser for arguments.
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
            description='''This maps KEGG KO identifiers into modules and summarizes the metabolic potential\n'''
            '''Usage: ''' + sys.argv[0] + ''' -i [Input Files] -p [Prefix Output]\n'''
            '''Global mandatory parameters: -i [Input Files] -p [Prefix Output]\n'''
            '''Optional Database Parameters: See ''' + sys.argv[0] + ' -h')
    parser.add_argument('-i', '--input_files', dest='input_files', nargs='+', action='store', required=True,
                        help='Space-separated list of files to parse.')
    parser.add_argument('-p', '--prefix', dest='prefix', action='store', required=True,
                        help='Prefix for the output_file files.')
    parser.add_argument('--cluster', dest='cluster', action='store', required=False,
                        help='Cluster genomes or modules. Select "cols" for genomes, "rows" for modules, or "both".\
                            By default no clustering.')
    args = parser.parse_args()

    input_files = args.input_files
    prefix = args.prefix
    cluster = args.cluster

    # ----------------------------
    # Import data
    regular_modules, bifurcation_modules, structural_modules,  \
    module_information, metabolism_matrix, module_group_matrix = module_information_importer(input_files)
    # Map and annotate genomes
    metabolic_annotation = global_mapper(regular_modules, bifurcation_modules, structural_modules, input_files)
    metabolism_matrix_dropped_relabel, module_colors = create_output_files(metabolic_annotation, metabolism_matrix, module_information, cluster, prefix)
    plot_function_barplots(module_colors, module_group_matrix, metabolism_matrix_dropped_relabel, prefix)
    # ----------------------------

if __name__ == "__main__":
    main()
