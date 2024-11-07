from time import sleep

import json
import requests
import traceback
import os
import random
import sys
import shutil

import subprocess as sp
import argparse

chopper = 'scripts/astrochop.py'
BFVD_DB = # path to bfvd

previous_job = ''

# HELPING ROUTINES

# -1. Load inputs

def load_inputs():

    # LOAD INPUTS

    parser = argparse.ArgumentParser(prog = 'FoldSeek_screen', usage = 'foldseek_screen -in <infile> [options]', 
                                         description = 'FoldSeek_screen is a python-based, local tool to screen input protein structures against the alphafold and ESMfold databases using foldseek.',
                                         epilog = 'Example: foldseek_screen -in <infile> ')

    requiredNamed = parser.add_argument_group('Required arguments')
    optionalNamed = parser.add_argument_group('Options with defaults')

    # required inputs
    requiredNamed.add_argument('-in', dest='targets', nargs='+', type=str, required=True, help='Input file(s) and or folder for all-against-all similarity searches')
    # optional inputs
    optionalNamed.add_argument('-n_iterations', dest='n_iterations', type=int, default = 1, help='Number of iterations to carry out (default: 1)')
    optionalNamed.add_argument('-iter_prob', dest='iter_prob', default = 0, type=float, help='Minimum probability of a match to be included for searches in the next iteration (default: 0)')
    optionalNamed.add_argument('-db', dest='db', default = ['afdb50', 'mgnify_esm30'], nargs='+', type=str, help='Databases to search over. (default: afdb50, mgnify_esm30')
    # optionalNamed.add_argument('-bfvd_str', dest='bfvd_str', default = None, type=str, help='Location of the BVFD structures. (default: None')

    # Define inputs
    return parser.parse_args()

def StartFoldSeek(inpdb, iterations = 1, curr_round = 1, databases = ['afdb50', 'mgnify_esm30'], min_prob_for_iteration = 0):
    
    label = inpdb
    if '-' in label:
        label = label.split('-')[1]
        
    outjson = '{}__foldseek_matches_{}_{}.json'.format(label, curr_round, '_'.join(databases))
    
    if not os.path.isfile(outjson):
        inpdb = parse_pdb(inpdb)
        results = {}
        
        for db in databases:

            url = 'https://search.foldseek.com/api/ticket'
            data = {'q': inpdb, 'mode': '3diaa', 'database[]': db}

            response = requests.post(url, data=data)

            if response.status_code != 200:
                return False

            foldseek_json = response.json()
            foldseek_id = foldseek_json['id']

            print(' ... Job submitted for {} in {} database (round {})'.format(label, db, curr_round))
            curr_results = GetFoldseekResults(foldseek_id)
            for key in curr_results:
                results[key] = curr_results[key]
        
        json.dump(results, open(outjson, 'w'), indent = 4)

    else:
        print(' ... Job was already ran for {}. Will load results only.'.format(label))
        results = json.load(open(outjson, 'r'))
        print(' ... ... Found {} hits'.format(len(results)))
        
    if curr_round < iterations:
        matched = list([f for key in results.keys() for f in results[key]['pdb'] if results[key]['prob'] >= min_prob_for_iteration])
        for i, hit in enumerate(matched):
            if hit != label:
                sleep(10)
                it_results = StartFoldSeek(hit, iterations = iterations, curr_round = curr_round + 1, databases = databases, min_prob_for_iteration = min_prob_for_iteration)

                if it_results:
                    for ac in it_results:
                        if ac not in results:
                            results[ac] = it_results[ac]
                else:
                    print('ERROR:', hit)
    
                print(' ... ... {} Total number of hits collected: {}'.format(i, len(results)))
 
    return results

def GetFoldseekResults(foldseek_id):
    
    url = f"https://search.foldseek.com/api/ticket/{foldseek_id}"
        
    repeat = True
    print(' ... ... Job is running')
    while repeat:
        response = requests.get(url).json()
        
        if response['status'] == "ERROR":
            return response['status']
        else:
            sleep(1)
            repeat = response['status'] != "COMPLETE"
    
    print(' ... ... Job is done')
    unique_hits = {}
    npdbs = 0
    
    print(' ... ... Parsing')
    try:
        # Its complete, get the results
        response = requests.get(f"https://search.foldseek.com/api/result/{foldseek_id}/0")
        results = response.json()
        alns = results['results'][0]['alignments'][0]

        if alns:
            
            for aln in alns:
                if aln['eval'] > .1:
                    continue

                if '-assembly' in aln['target']:
                    key = aln['target'].split('-assembly')[0]
                    chain = aln['target'].split('_')[1].split(' ')[0]
                    key = '{}_{}'.format(key, chain)
                else:
                    try:
                        key = aln['target'].split('-')[1]
                    except:
                        key = aln['target'].split('.')[0]
                    
                if key not in unique_hits:
                    hit = {
                        'hit_seq': aln['dbAln'],
                        'query_seq': aln['qAln'],
                        'hit_seq_length': aln['dbLen'],
                        'hit_offset' : aln['dbStartPos']-1,
                        'query_offset' : aln['qStartPos']-1,
                        'seqid' : aln['seqId'],
                        'score' : aln['score'],
                        'eval' : aln['eval'],
                        'prob' : aln['prob'],
                        'coverage' : round(float(aln['qEndPos']-aln['qStartPos']+1)/aln['qLen'], 2),
                        'chains': []
                    }
                    if len(aln['query'].split('_')) > 1:
                        hit['chains'] = [aln['query'].split('_')[1]]
                    try:
                        hit['pdb'] = GetMatchedDomains(key, aln['dbStartPos'],aln['dbEndPos'])
                        npdbs += len(hit['pdb'])
                        unique_hits[key] = hit
                    except Exception as ex:
                        print(ex)
                        print(traceback.format_exc())
                                    
                else:
                    if len(aln['query'].split('_')) > 1:
                        unique_hits[key]['chains'].append(aln['query'].split('_')[1])
                    
    except Exception as ex:
        print(ex)
        print(traceback.format_exc())          
    
    print(' ... ... Found {} hits over {} unique domains'.format(len(unique_hits), npdbs))
    
    return unique_hits

def GetMatchedDomains(code, starti, endi, min_len = 40, min_cov = 0.1):

    pdb = download_pdbs(code)
    domains = chop_domains(pdb)
    
    target_domains = select_domains(domains, starti, endi, min_len = min_len, min_cov = min_cov)
    target_pdbs = save_regions_to_pdb(pdb, target_domains)
            
    return target_pdbs
    

def parse_pdb(inpdb, source_db = 'AFDB'):
    
    if os.path.isfile(inpdb):
        outp = ''
        with open(inpdb, 'r') as inp:
            for line in inp:
                outp += line
        outp = outp.encode('utf-8')
                
    elif source_db == 'AFDB':
        inpdb = download_AF2_pdbs(inpdb)
        outp  = parse_pdb(inpdb)        
        
    return outp

def download_pdbs(code, outfolder = '{}AF2_models'.format(previous_job), bfvd_str = BFVD_DB):
    
    def download_AF2(url: str, dest_folder: str):

        filename = url.split('/')[-1].replace(" ", "_")  # be careful with file names
        file_path = os.path.join(dest_folder, filename)
        
        r = requests.get(url, stream=True)
        if r.ok:
#             print(" ... ... Saving", os.path.abspath(file_path))
            with open(file_path, 'wb') as f:
                for chunk in r.iter_content(chunk_size=1024 * 8):
                    if chunk:
                        f.write(chunk)
                        f.flush()
                        os.fsync(f.fileno())
        else:  # HTTP status code 4XX/5XX
#             print(" ... ... Download failed ({}): status code {}\n{}".format(url, r.status_code, r.text)) 
            pass
        
        return r.ok
    
    def download_ESM(mgy_link: str, dest_folder: str, outmodel: str):
        
        download_structue = sp.Popen(['curl', '-o', outmodel, mgy_link, '-k'], stdout=sp.PIPE, stderr=sp.PIPE, stdin=sp.PIPE)
        stdout, stderr = download_structue.communicate()
                            
        return stderr

    def download_PDB_chain(pdb_link: str, dest_folder: str, outmodel: str, chain: str):
        
        download_structue = sp.Popen(['curl', pdb_link, '-k'], stdout=sp.PIPE, stderr=sp.PIPE, stdin=sp.PIPE)
        stdout, stderr = download_structue.communicate()

        stdout = stdout.decode().split('\n')
        if len(stdout) > 1:
            with open(outmodel,'w') as outp:
                for line in stdout:
                    if line.startswith('ATOM ') and (chain is None or line[20:23].strip() == chain):
                        if chain is not None and chain != 'A':
                            # change the chain id to A because chainsaw only knows how to deal with chains A
                            line = list(line)
                            line[21] = 'A'
                            line = ''.join(line)
                        outp.write('{}\n'.format(line))
                            
        return stderr

    #---------------
        
    if not os.path.isdir(outfolder):
        os.mkdir(outfolder)
    
    if 'MGYP' in code:
        outmodel = '{}/{}.pdb'.format(outfolder, code)
        mgy_link = 'https://api.esmatlas.com/fetchPredictedStructure/{}'.format(code)
        
        if not os.path.isfile(outmodel):
            worked = download_ESM(mgy_link, dest_folder=outfolder, outmodel=outmodel)

    elif 'unrelaxed' in code:
        outmodel = '{}/{}.pdb'.format(outfolder, code)
        if bfvd_str is not None:
            worked = '{}/{}.pdb'.format(bfvd_str, code) 
            if os.path.isfile(worked):
                with open(worked, 'r') as inp:
                    with open(outmodel, 'w') as outp:
                        for line in inp:
                            if line.startswith('ATOM '):
                                outp.write(line)

    elif '_' in code:
        outmodel = '{}/{}.pdb'.format(outfolder, code)
        pdb_link = 'https://files.rcsb.org/download/{}.pdb'.format(code.split('_')[0])
        if not os.path.isfile(outmodel):
            worked = download_PDB_chain(pdb_link, dest_folder=outfolder, outmodel=outmodel, chain=code.split('_')[1])
            
    else:   
        outmodel = '{}/AF-{}-F1-model_v4.pdb'.format(outfolder, code)
        af2_link = 'https://alphafold.ebi.ac.uk/files/AF-{}-F1-model_v4.pdb'.format(code)
        
        if not os.path.isfile(outmodel):
            worked = download_AF2(af2_link, dest_folder=outfolder)
    
    if not os.path.isfile(outmodel):
        return None
    
    return outmodel

def chop_domains(pdb, chopper = chopper):
    
    def parse_domains(outfile):
        with open(outfile, 'r') as chopper_outp:
            for line in chopper_outp:
                if line.startswith('RESULT'):
                    return line.strip().split('\t')[-1].split(',')        
    #--------------
    
    outfile = '{}.astrochop'.format(pdb)
    print(' ... ... ... Running astrochop for {}'.format(pdb))
    
    if not os.path.isfile(outfile):
        chop = sp.Popen(['python3', chopper, pdb], stdout=sp.PIPE, stderr=sp.PIPE, stdin=sp.PIPE)
        stdout, stderr = chop.communicate()

        print(stderr.decode())
        
        with open(outfile, 'w') as outp:
            outp.write(stdout.decode())
    
    outdomains = parse_domains(outfile)
    
    return outdomains

def select_domains(domains, starti, endi, min_len = 40, min_cov = 0.1):

    match_range = set(range(starti, endi))
    
    selected_domains = {}
    for i, domain in enumerate(domains):
        parts = domain.split('_')
        for part in parts:
            part = [int(x) for x in part.split('-')]
            
            if len(part)>1:
                part_range = set(range(part[0], part[1]))

                overlap = match_range.intersection(part_range)
                coverage = len(overlap)/len(match_range)

                if len(part_range) >= min_len and coverage > min_cov:
                    selected_domains[i] = domain
    
    return selected_domains

def save_regions_to_pdb(pdb, domains, outdir='AF2_models_selected_regions'):
    
    def parse_pdb(pdb):
        coordinates = {}
        with open(pdb, 'r') as inpdb:
            for line in inpdb:
                if line.startswith('ATOM '):
                    residue_num = int(line[22:27].strip())
                    if residue_num not in coordinates:
                        coordinates[residue_num] = []
                    coordinates[residue_num].append(line)
        return coordinates

    #--------------------------
    
    coordinates = parse_pdb(pdb)
    saved_pdbs = []
    
    if not os.path.isdir(outdir):
        os.mkdir(outdir)
        
    for domain in domains:
        outfile = '{}/{}.domain{}.pdb'.format(outdir, pdb.split('/')[-1].replace('.pdb', ''), domain)
        interval = domains[domain]
        saved_pdbs.append(outfile)
        
        if not os.path.isfile(outfile):
            with open(outfile, 'w') as outpdb:
                for part in interval.split('_'):
                    part = [int(i) for i in part.split('-')]
                    for j in range(part[0], part[1]+1):
                        if j in coordinates:
                            for atom in coordinates[j]:
                                outpdb.write(atom)
                        
    return saved_pdbs


if __name__ == '__main__':

    args = load_inputs()

    for target in args.targets:
        if os.path.isdir(target):
            for curr_target in os.listdir(target):
                f = '{}/{}'.format(target,curr_target)
                all_hits = StartFoldSeek(f, iterations = args.n_iterations, databases = args.db, min_prob_for_iteration = args.iter_prob)
                json.dump(all_hits, open('all_hits_{}_{}.json'.format(curr_target, '_'.join(args.db)), 'w'), indent=4)
        else:
            all_hits = StartFoldSeek(target, iterations = args.n_iterations, databases = args.db, min_prob_for_iteration = args.iter_probZ)
            json.dump(all_hits, open('all_hits_{}_{}.json'.format(target, '_'.join(args.db)), 'w'), indent=4)