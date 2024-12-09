#!/usr/bin/env python3

import argparse
import os
import sys
import shutil
import string
import random
import subprocess

##### INPUTS
software_description = '''
                        Generate a multiple sequence alignment (MSA) from query and sequence collection.
                        Wrapper for RepGenR Phylo
                        '''

### Parse input arguments
## setup
parser = argparse.ArgumentParser(description=software_description)

parser.add_argument('-q','--query',required=True,help='Path to query file')
parser.add_argument('-c','--collection',required=True,help='Path to sequence collection')
parser.add_argument('-o','--output',required=True,help='Path to output MSA file (e.g. msa.fasta)')
parser.add_argument('--progressivemauve_single',action='store_true',help='If specified, will not parallelize progressivemauve. Sometimes this leads to a yet unknown error. Running one process at a time has solved the issue.')
parser.add_argument('--query_as_backbone',action='store_true',help='If specified, will use the query as backbone in MSA. Only apply this option if you have a single-contig genome (default: use a reference-sequence as backbone)')


##/
## parse input
if 1 and 'run':
    args = parser.parse_args()
else:
    print('IDE MODE',flush=True)
    if 1:
        args = parser.parse_args(['--query','genomes/Sample01_idbaud.fasta',
                                  '--collection','collection',
                                  '--output','msa.fasta'])

query_file_path = args.query
collection_path = args.collection
output_file_path = args.output
recombination_dataset = args.query

progressivemauve_single = args.progressivemauve_single
query_as_backbone = args.query_as_backbone
##/
## format
query_file_path = os.getcwd() + '/' + query_file_path
collection_path = os.getcwd() + '/' + collection_path
##/
#####/

### Determine query basename
query_basename = os.path.basename(query_file_path)
query_dirname = os.path.dirname(query_file_path)
query_basename_noExt = os.path.splitext(query_basename)[0]
###/

### init tmpdir
tmpdir_path = None
while tmpdir_path == None:
    # generate random name
    random_chars_to_use = string.ascii_letters + string.digits # setup method to generate (characters + digits)
    random_name = ''.join(random.choice(random_chars_to_use) for x in range(15)) # generate random string of X length
    #/
    # compile a candidate path for tmpdir
    candidate_tmpdir_path = 'tmp_recomfi_msa_' + random_name
    #/
    # if this path does not exist then we can use it as tmpdir
    if not os.path.exists(candidate_tmpdir_path):
        tmpdir_path = os.getcwd() + '/' + candidate_tmpdir_path
    #/

if not os.path.exists(tmpdir_path):          os.makedirs(tmpdir_path)
###/

### Setup "genomes" folder inside tmpdir, link-in genomes (uncompressed) or gzip-in genomes (compressed); progressivemauve does not support gzipped files
# init genomes dir
genomes_dir = tmpdir_path+'/'+'genomes'
os.makedirs(genomes_dir)
#/

## link-in or gzip-decompress-in files
# function
def link_or_decompress_file_into_genomes_dir(file_=None,file_source_dir=None,file_target_dir=None):
    if file_.endswith('.gz'):
        print(f'Decompress-in file to tmpdir: {file_}')
        file_basename = os.path.splitext(file_.replace('.gz',''))[0]
        gzip_cmd = ['gzip',
                    '-cd',
                    file_source_dir+'/'+file_,
                    '>',file_target_dir+'/'+file_basename+'.fasta']
        subprocess.call(' '.join(map(str,gzip_cmd)),shell=True)
    else:
        print(f'Link-in file to tmpdir: {file_}')
        file_basename = os.path.splitext(file_)[0]
        os.symlink(file_source_dir+'/'+file_,file_target_dir+'/'+file_basename+'.fasta')
#/
# link-in or gzip-in collection
for file_ in os.listdir(collection_path):
    link_or_decompress_file_into_genomes_dir(file_=file_,file_source_dir=collection_path,file_target_dir=genomes_dir)
#/
# link-in or gzip-in query
link_or_decompress_file_into_genomes_dir(file_=query_basename,file_source_dir=query_dirname,file_target_dir=genomes_dir)
#/
##/
###/

### Determine which genome to use as reference in MSA build
pmauve_ref = None
for file_ in os.listdir(genomes_dir):
    # unless query should be the backbone, use a "random" reference in the collection as backbone
    if not query_as_backbone and file_.find(query_basename_noExt) == -1:
        pmauve_ref = file_
        break # break on first
    #/
    # else, use the query
    else:
        if file_.find(query_basename_noExt) != -1:
            pmauve_ref = file_
    #/

if pmauve_ref == None:
    print('WARNING: Could not determine a reference for progressivemauve. Please ensure that your collection consist of sequence files.')
    print('Terminating!')
    sys.exit()

print(f'Assigned {pmauve_ref} as reference in progressivemauve MSA build')
###/

### Call RepGenR phylo
## Compile CMD
repgenr_phylo_cmd = ['phylo.py',
                     '--mode','accurate',
                     '--all_genomes',
                     '--no_outgroup',
                     '--workdir',tmpdir_path,
                     '--keep_msa',
                     '--halt_after_msa',
                     '--progressivemauve_ref',pmauve_ref]

# check if user specified to single-process progressivemauve
if progressivemauve_single:
    repgenr_phylo_cmd.append('--progressivemauve_single')
#/
##/
## Call CMD
subprocess.call(' '.join(map(str,repgenr_phylo_cmd)),shell=True)
##/
## Validate output exists
msa_output_expected_path = tmpdir_path+'/'+'msa_all.fasta'
if not (os.path.exists(msa_output_expected_path) and os.path.getsize(msa_output_expected_path) > 0):
    print('WARNING: MSA not produced as expected. Please verify that all dependencies have been installed, such as RepGenR.')
    print('Terminating!')
    sys.exit()
##/
###/

### Move-out MSA file and remove tmpdir
os.rename(tmpdir_path+'/'+'msa_all.fasta',output_file_path)
shutil.rmtree(tmpdir_path)
###/
