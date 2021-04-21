from pathlib import Path
from microbeannotator.database.kofamscan_data import parse_profiles
from microbeannotator.database.kofamscan_data import kofamscan_downloader
from microbeannotator.database.kofamscan_data import kofamscan_data_formatter
from microbeannotator.pipeline.hmmsearch import hmmsearch_launcher
from microbeannotator.pipeline.hmmsearch import protein_to_model_mapper
search = 'prokaryote'


# Download and parse databases
# model_info_path, model_dir = kofamscan_downloader(
#     Path('/mnt/c/Users/Cruiz/Desktop/00.Proteins/02.tests'))
# kofamscan_data_formatter(model_dir)

model_dir = Path("/mnt/c/Users/Cruiz/Desktop/00.Proteins/02.tests/kofam_data/profiles")
model_info_path = Path("/mnt/c/Users/Cruiz/Desktop/00.Proteins/02.tests/kofam_data/ko_list")
# Get hmm model ids depending on subset
model_list = []
with open(f"{model_dir / 'common.list'}", 'r') as infile:
    for line in infile:
        model_list.append(model_dir / line.strip())
with open(f"{model_dir / 'independent.list'}", 'r') as infile:
    for line in infile:
        model_list.append(model_dir / line.strip())
if search == 'prokaryote':
    with open(f"{model_dir / 'prokaryote.list'}", 'r') as infile:
        for line in infile:
            model_list.append(model_dir / line.strip())
elif search == 'eukaryote':
    with open(f"{model_dir / 'eukaryote.list'}", 'r') as infile:
        for line in infile:
            model_list.append(model_dir / line.strip())

# Get protein files to search
protein_list = [
    Path('/mnt/c/Users/Cruiz/Desktop/00.Proteins/01.proteins/GCA_000005845.2.faa'),
    Path('/mnt/c/Users/Cruiz/Desktop/00.Proteins/01.proteins/GCA_000006175.2.faa'),
    Path('/mnt/c/Users/Cruiz/Desktop/00.Proteins/01.proteins/GCA_002763915.1.faa'),
    Path('/mnt/c/Users/Cruiz/Desktop/00.Proteins/01.proteins/GCA_004154025.1.faa')]

output_dir = Path('/mnt/c/Users/Cruiz/Desktop/00.Proteins/02.tests')
# Combine the lists of protein and model paths
for protein_file in protein_list:
    list_proteins_models = protein_to_model_mapper([protein_file], model_list)
    output = output_dir / (Path(protein_file.name).with_suffix('.kofam'))

    combined_hmmsearch_results = hmmsearch_launcher(
        protein_model_list=list_proteins_models,
        profile_info_path=model_info_path,
        output_path=output,
        processes=2, threads=8)
