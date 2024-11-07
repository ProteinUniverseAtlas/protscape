import subprocess as sp
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx

from tqdm import tqdm
from scipy import stats

import sys
import os
import argparse
import time
import json
import random
import pacmap
import numba
import matplotlib
import shutil

# # to setup the environment for protscape.py
# source "${HOME}/conda/etc/profile.d/conda.sh"
# source "${HOME}/conda/etc/profile.d/mamba.sh"
# mamba create -n protscape -c bioconda -c salilab -c conda-forge -c anaconda python=3.11 foldseek dssp mmseqs2 libboost=1.73.0 numpy matplotlib networkx tqdm pacmap
# mamba activate protscape
# cd code/neighborhood_correlation/
# python setup.py install

# # to initiate the environment for protscape.py
# mamba activate protscape

# HELPING ROUTINES

# -1. Load inputs

def load_inputs():

    # LOAD INPUTS

    parser = argparse.ArgumentParser(prog = 'ProtScape', usage = 'protscape -infile <infile> [options]', 
                                         description = 'ProtScape is a python-based, local tool to generate protein sequence or structire similarity networks, based on MMseqs2 or FoldSeek.',
                                         epilog = 'Example: protscape -infile <infile>')

    requiredNamed = parser.add_argument_group('Required arguments')
    optionalNamed = parser.add_argument_group('Options with defaults')

    # required inputs
    requiredNamed.add_argument('-infile', dest='infile', type=str, required=True, help='Input file or folder for all-against-all similarity searches')
    # optional inputs
    optionalNamed.add_argument('-max_evalue', dest='max_evalue', type=float, default = 1e-4, help='E-value threshold for match collection (default: 1e-4)')
    optionalNamed.add_argument('-find_clusters', dest='find_clusters', type=bool, default = False, help='Boolean statement to find clusters automatically using networkx (deafult: False)')
    optionalNamed.add_argument('-n_iterations', dest='n_iterations', type=int, default = 5, help='Maximum number of iterations to find clusters iteratively with networkx (default: 5)')
    optionalNamed.add_argument('-evalue_clustering', dest='evalue_clustering', type=float, default = 1e-4, help='Initial E-value threshold for finding clusters with networkx (default: 1e-4)')
    optionalNamed.add_argument('-evalue_sigma', dest='evalue_sigma', type=float, default = 0, help='Sigma coefficient for iteratively setting the E-value threshold for finding clusters with networkx (default: 0)')
    optionalNamed.add_argument('-cluster_size', dest='cluster_size', type=int, default = 10, help='Minimum size of clusters to keep (default: 10)')
    optionalNamed.add_argument('-min_nc', dest='min_nc', type=float, default = 0.05, help='Minimum neighborhood correlation to be considerd (default: 0.05)')
    optionalNamed.add_argument('-tmp_folder', dest='tmp_folder', default = '/tmp',type=str, help='The temporary folder (default: /tmp)')
    optionalNamed.add_argument('-mode', dest='mode', default = 'mmseqs',type=str, help='The method to use (default: mmseqs). Allowed: mmseqs, foldseek, mmseqs_NC, foldseek_NC')
    optionalNamed.add_argument('-plddt_weighted', dest='plddt_weighted', default=False,type=bool, help='Flag to weight scores by the pLDDT of the matched proteins')
    optionalNamed.add_argument('-infasta', dest='infasta', default = None,type=str, help='An input fasta file if no sequences from structure are to be used')
    optionalNamed.add_argument('-inseqcon', dest='inseqcon', default = None,type=str, help='An input all-against-all mmseqs-like file to replace the all-against-all step. Inner structure should be "query,target,evalue"')

    # Define inputs
    return parser.parse_args()

# 0. General

def sequence_index(infasta, ref_fasta = None):

    seq_index = {}
    count = -1

    if os.path.isfile(infasta):
        with open(infasta, 'r') as inf:
            for line in inf:
                if line.startswith('>'):
                    if count > -1:
                        seq_index[curr_code] = {'description': curr_header.strip(), 'sequence': curr_seq, 'index': count}

                    curr_code = line.strip('>').split()[0]
                    if 'tr|' in curr_code or 'sp|' in curr_code:
                        curr_code = curr_code.split('|')[1]
                    elif '|' in curr_code:
                        curr_code = curr_code.split('|')[0].strip()

                    curr_header = line.strip('>')
                    curr_seq = ''
                    count += 1
                else:
                    curr_seq += line.strip()

            seq_index[curr_code] = {'description': curr_header, 'sequence': curr_seq, 'index': count}

        if count >= len(seq_index):
            print('\nWARNING: The number of sequences read does not match the total number of sequences in the file.\nMaybe there are duplicates in the file. Check!\n')

    elif os.path.isdir(infasta):
        if ref_fasta is None:
            outfasta = '{}.fasta'.format(infasta)
            if not os.path.isfile(outfasta):
                with open(outfasta, 'w') as outp:
                    dir_content = [file for file in os.listdir(infasta) if file.endswith('.pdb')]
                    for file in tqdm(dir_content, total=len(dir_content), desc=' ... Collecting sequences',  bar_format="{l_bar}{bar}{r_bar} [ time left: {remaining}, time spent: {elapsed}]"):
                        sequence = get_sequence_structure('{}/{}'.format(infasta,file))
                        outp.write('>{}\n{}\n'.format(file, sequence))

            seq_index, infasta = sequence_index(outfasta)

        else:
            seq_index, infasta = sequence_index(ref_fasta)

    return seq_index, infasta

def get_header_map(infile, seq_db):
    
    header_map = dict()

    dir_content = [file for file in os.listdir(infile) if file.endswith('.pdb')]
    for file in tqdm(dir_content, total=len(dir_content), desc=' ... Creating header map',  bar_format="{l_bar}{bar}{r_bar} [ time left: {remaining}, time spent: {elapsed}]"):

        try:    
            curr_map = [i for i in seq_db if file.split('.')[0] in i.split('.')[0]][0]
            header_map[file] = curr_map
        except:
            pass

    return header_map

def get_sequence_structure(pdb_file):

    out_file = "{}_dssp.out".format(pdb_file)

    if not os.path.isfile(out_file):
        
        run_dssp = sp.Popen(['mkdssp', '-i', pdb_file, '-o', out_file], stdout=sp.PIPE, stderr=sp.PIPE, stdin=sp.PIPE)
        stdout, stderr = run_dssp.communicate() 
        if len(stderr) > 0:
            run_dssp = sp.Popen(['dssp', '-i', pdb_file, '-o', out_file], stdout=sp.PIPE, stderr=sp.PIPE, stdin=sp.PIPE)
            stdout, stderr = run_dssp.communicate() 
            if len(stderr) > 0:
                print(stderr)
                exit()
                   
    dssp_aa = parse_DSSPout(out_file)['AA']
    sequence = ''
    
    for i in range(len(dssp_aa)):
        sequence += dssp_aa[i]

    os.remove(out_file)
    
    return sequence

def parse_DSSPout(dsspout_file):

    found_table = False
    data_dict = {'ResNum':[], 'AA': [], 'SecStr': [], 'Phi': [], 'Psi': [], 'Chains': []}
    
    with open(dsspout_file, "r") as dsspout:
        for line in dsspout:
            if "#  RESIDUE AA STRUCTURE BP1" in line:
                found_table = True
            elif found_table:
                resnum = line[5:10].strip()
                secstr = line[16]
                secstr_compl = line[23]
                
                if secstr in [' ', 'S', 'B', 'T']:
                    secstr = 'L'
                elif secstr in ['G', 'I']:
                    secstr = 'H'

                if secstr == 'L':
                    secstr = '-'
                    
                res = line[13]
                phi = line[103:109].strip()
                psi = line[110:116].strip()
                chain = line[10:13].strip()

                if res != '!':
                    data_dict['ResNum'].append(resnum)
                    data_dict['AA'].append(res)
                    data_dict['SecStr'].append(secstr)
                    data_dict['Phi'].append(phi)
                    data_dict['Psi'].append(psi)
                    data_dict['Chains'].append(chain)

    for i in range(len(data_dict['SecStr'])):
        if i > 0 and i < len(data_dict['SecStr'])-1:
            if data_dict['SecStr'][i] == '-':
                if data_dict['SecStr'][i-1] != '-' and data_dict['SecStr'][i+1] != '-':
                    if data_dict['SecStr'][i-1] == data_dict['SecStr'][i+1]:
                        data_dict['SecStr'][i] = data_dict['SecStr'][i-1] 
    return data_dict

def replace_tabs(intable):

    lines = [line for line in open(intable,'r')]

    with open(intable, 'w') as outp:
        for line in lines:
            line = line.split()
            line = ' '.join(line)
            outp.write('{}\n'.format(line))

    return intable

def extract_plddts_from_structures(pdbs):

    print(' ... ... Extracting pLDDTs from {} structures'.format(len(pdbs)))

    plddts = dict()

    for f in tqdm(pdbs, total=len(pdbs), desc=' ... ... Extracting pLDDTs',  bar_format="{l_bar}{bar}{r_bar} [ time left: {remaining}, time spent: {elapsed}]"):
        key = f.split('/')[-1]
        plddts[key] = []

        with open(f, 'r') as inf:
            for line in inf:
                if line.startswith('ATOM '):
                    plddts[key].append(float(line[60:67]))
        plddts[key] = np.mean(plddts[key])

    return plddts


# 1. For link finding 

def all_against_all(infile, mode = 'mmseqs', tmp_folder = '/tmp/', n_iterations = 1, 
                    max_evalue = 1e-4, sensitivity = 7.5, coverage = 0, n_cpus = 2, infasta = None):

    print('\nRunning all-against-all')

    outpath = '{}_results'.format(mode)
    if not os.path.isdir(outpath):
        os.mkdir(outpath)

    searchresults = '{}/{}.{}'.format(outpath, infile.split('/')[-1].split('.')[0], mode)

    if '_' in mode:
        mode = mode.split('_')[0]
        outp_format = 'query,target,bits'
    elif 'foldseek' in mode:
        outp_format = 'query,target,evalue,qtmscore,ttmscore,alntmscore'
    else:
        outp_format = 'query,target,evalue'

    if mode == 'mmseqs' and os.path.isdir(infile):
        _, infile = sequence_index(infile)

    if mode == 'foldseek' and os.path.isfile(infasta):
        
        print('\n ... Generating targets folder')

        targets = '{}/{}_targets'.format(tmp_folder, infile.split('/')[-1].split('.')[0])
        if not os.path.isdir(targets):
            os.mkdir(targets)

        with open(infasta, 'r') as inf:
            for line in inf:
                if line.startswith('>'):
                    pdb = line.strip().split('>')[1]
                    if not pdb.endswith('.pdb'):
                        pdb += '.pdb'

                    shutil.copyfile('{}/{}'.format(infile, pdb), '{}/{}'.format(targets, pdb))

        infile = targets

    if not os.path.isfile(searchresults):
        if mode == 'mmseqs':
            sp.run(['mmseqs', 'easy-search', infile, infile, searchresults, tmp_folder, 
                '-e', str(max_evalue), '-s', str(sensitivity), '-c', str(coverage), 
                '--num-iterations', str(n_iterations), '--threads', str(n_cpus), 
                '--format-output', outp_format])

        elif mode == 'foldseek':
            sp.run(['foldseek', 'easy-search', infile, infile, searchresults, tmp_folder, 
                '-e', str(max_evalue), '-s', str(sensitivity), '-c', str(coverage), 
                '--threads', str(n_cpus), 
                '--format-output', outp_format])

    return searchresults

def neighborhood_correlation(intable, min_nc = 0.05):

    intable = replace_tabs(intable)

    searchresults = '{}.ncorr'.format(intable)

    if not os.path.isfile(searchresults):
        sp.run(['NC_standalone', '-f', intable, '--nc_thresh', str(min_nc), '-o', searchresults])

    # convert NC to e-value-like
    lines = [line for line in open(searchresults, 'r')]
    with open(searchresults, 'w') as outp:
        for line in lines:
            line = line.split()
            line[-1] = str(1-float(line[-1]))
            # line[-1] = str((2/(1+np.exp(-5*float(line[-1]))))-1)
            line = '\t'.join(line)
            outp.write('{}\n'.format(line))

    return searchresults

def weight_by_plddt(infile, seq_con):

    print(' ... Weighting by plddt')

    if os.path.isdir(infile):

        outfile = '{}.plddt_weighted'.format(seq_con)

        if not os.path.isfile(outfile):

            plddts = extract_plddts_from_structures(['{}/{}'.format(infile, f) for f in os.listdir(infile) if f.endswith('.pdb')])
            json.dump(plddts, open('{}_plddts.json'.format(infile), 'w'), indent = 4)

            print(' ... ... Weighting')
            with open(outfile, 'w') as outf:
                with open(seq_con, 'r') as inf:
                    for line in inf:
                        i, j, distance = line.strip().split()

                        weight = np.mean([plddts[i], plddts[j]])
                        weighted_distance = float(distance)/weight

                        outf.write('{}\t{}\t{}\n'.format(i, j, weighted_distance))

        seq_con =  outfile
    
    return seq_con

# 2. For writting CLANS-like output

# 2.1. To define clusters

def get_edges(edge_data, evalue_clustering, name_prefix = 0):

    # edges are saved as a dictionary whose value is proportional to the evalue
    edges = {}

    # if edge_data is a file, it is an mmseqs/foldseek result file that needs to be parsed
    try:
        if os.path.isfile(edge_data):
            mmseqs_results = edge_data
            with open(mmseqs_results, 'r') as curr_thread:
                for line in curr_thread:
                    line = line.split()
                    if len(line) >= 0:
                        in_node = line[0]
                        out_node = line[1]
                        dist = float(line[2])
                        if dist > 0 and dist <= evalue_clustering:
                            edges[(in_node,out_node)] = -np.log10(dist)

    # if edge_data is not a file, then it is a networks graph object with weighted edges
    except:
        for in_node,out_node,data in edge_data.edges(data=True):
            if data['weight'] >= evalue_clustering:
                edges[(in_node,out_node)] = data['weight']

    # weights = list(edges.values())
    # weight_threshold = np.median(weights) + stats.median_abs_deviation(weights)

    # plt.clf()
    # plt.hist(weights, bins = round(len(set(weights))/2), range=(min(weights), max(weights)))
    # plt.vlines(x = weight_threshold, ymin=0, ymax=200, color='red')
    # plt.savefig('weights_hist_{}.png'.format(name_prefix))

    return edges

def get_connected_components(edges, min_size):

    graph = nx.from_edgelist(list(edges.keys()))
    nx.set_edge_attributes(graph, values = edges, name = 'weight')

    connected_components = [graph.subgraph(c).copy() for c in nx.connected_components(graph) if len(c) >= min_size]

    return connected_components

def extract_clusters(mmseqs_results, evalue_clustering = 1e-10, min_size = 10, sigma = 0, name_prefix = '', n_iterations = 5, curr_iteration = 1):
    
    print('\nExtracting clusters')
    print('... Iteration {} for name "{}" (Threshold: {}) ...'.format(curr_iteration, name_prefix, evalue_clustering))

    edges = get_edges(edge_data = mmseqs_results, evalue_clustering = evalue_clustering, name_prefix = name_prefix)
    clusters = get_connected_components(edges, min_size)

    if curr_iteration <= n_iterations:
        all_clusters = []
        for i, cluster in enumerate(clusters):
            if len(cluster) > min_size*2:

                name = name_prefix + '{}.'.format(i)

                cluster_weights = [data['weight'] for node1,node2,data in cluster.edges(data=True)]
                evalue_clustering = np.median(cluster_weights) + stats.median_abs_deviation(cluster_weights)*sigma
                # evalue_clustering = np.mean(cluster_weights) + np.std(cluster_weights)

                subclusters = extract_clusters(mmseqs_results = cluster, evalue_clustering = evalue_clustering, min_size = min_size, sigma = sigma, name_prefix = name, n_iterations = n_iterations, curr_iteration = curr_iteration+1)

                if len(subclusters) <= 1:
                    subclusters = [cluster]

            else:
                subclusters = [cluster]

            all_clusters += subclusters

    else:
        all_clusters = clusters

    return all_clusters

# 2.2. To write the file

def get_clusters_colors(clusters, cmap='Spectral'):

   if type(cmap) == str:
      cmap = matplotlib.cm.get_cmap(cmap)
      norm = matplotlib.colors.Normalize(vmin=0, vmax=len(clusters))

      colours = [cmap(norm(i)) for i in range(len(clusters))]

   else:
      colours = [matplotlib.colors.to_rgb(cmap[c])+(1,) for c in cmap]

   colours = [[int(255*j) for j in colours[i]] for i in range(len(colours))]
   
   return colours

def write_clans_file(infile, mmseqs_results, infasta = None, mode = 'mmseqs', outpath = None, weighted=False, evalue_clustering = 1e-4, clusters = None, cmap = 'Spectral'):

    print('\nWritting clans map')

    seq_db, infasta = sequence_index(infile, infasta)

    if os.path.isdir(infile):
        header_map = get_header_map(infile, seq_db)
    else:
        header_map = {key: key for key in seq_db}

    print()

    if outpath is not None:
        if not os.path.isdir(outpath):
            os.mkdir(outpath)
        outfile = '{}/{}_{}.clans'.format(outpath, infasta.split('/')[-1].split('.')[0], mode)
    else:
        outfile = '{}_{}.clans'.format(infasta.split('/')[-1].split('.')[0], mode)

    if weighted:
        outfile = outfile.replace('.clans', '_plddt_weighted.clans')

    no_seq = len(seq_db)
    count = 0
    with open(outfile, 'w') as outp:
        
        # write parameters
        print(" ... Writing parameters")
        outp.write('sequences={}\n'.format(no_seq))
        if clusters is not None:
            outp.write('<param>\ndotsize=0\npval={}\ncluster2d=true\n</param>\n'.format(evalue_clustering))
        else:
            outp.write('<param>\ndotsize=5\npval={}\ncluster2d=true\n</param>\n'.format(evalue_clustering))
            
        # write sequences
        print(" ... Writing sequences")
        print(" ... ... There are {} sequences".format(no_seq))
        outp.write('<seq>\n')
        for i in seq_db:
            outp.write('>{}\n{}\n'.format(seq_db[i]['description'], seq_db[i]['sequence']))
        outp.write('</seq>\n')

        # write random positions for each sequence
        print(" ... Writing coordinates")
        outp.write('<pos>\n')
        for i in seq_db:
            x,y,z = round(random.uniform(-1, 1), 3), round(random.uniform(-1, 1), 3), round(random.uniform(1, 1), 3)
            outp.write('{} {} {} {}\n'.format(seq_db[i]['index'], round(x, 3), round(y, 3), round(z, 3)))
        outp.write('</pos>\n')

        # write hsps
        print(" ... Writing HSPs")
        outp.write('<hsp>\n')
        with open(mmseqs_results, 'r') as curr_thread:
            for line in curr_thread:
                line = line.split()
                if len(line) > 0:
                    in_node = seq_db[header_map[line[0]]]['index']
                    out_node = seq_db[header_map[line[1]]]['index']

                    outp.write('{} {}:{}\n'.format(in_node, out_node, line[2].replace('E-', 'e-')))
        outp.write('</hsp>\n')

        if clusters is not None:
            print(" ... Writing clusters")
            no_nodes = len([i for cluster in clusters for i in cluster])
            print(' ... ... There are {} connected components/clusters counting for {} nodes ({}% of all nodes)'.format(len(clusters), no_nodes, round(no_nodes*100/no_seq, 1)))

            colors = get_clusters_colors(clusters, cmap = cmap)
            outp.write('<seqgroups>\n')
            all_members = ''
            for i, cluster in enumerate(clusters):
                name = 'cluster_{}'.format(i)
                color = ";".join([str(col) for col in colors[i]])
                members = ';'.join([str(seq_db[node]['index']) for node in cluster])
                members = members + ";"

                outp.write('name={}\n'.format(name))
                outp.write('type=0\n')
                outp.write('size=9\n')
                outp.write('color={}\n'.format(color))
                outp.write('numbers={}\n'.format(members))

                all_members += members
            
            # write the 'all' group
            color = "0;0;0;255"
            outp.write('name=all\n')
            outp.write('type=0\n')
            outp.write('size=12\n')
            outp.write('color={}\n'.format(color))
            outp.write('numbers={}\n'.format(all_members))

            # write the singletons group
            color = "177;177;177;255"
            outp.write('name=all\n')
            outp.write('type=0\n')
            outp.write('size=6\n')
            outp.write('color={}\n'.format(color))
            outp.write('numbers={}\n'.format(';'.join([str(i) for i in range(0, no_seq)])+';'))

            outp.write('</seqgroups>\n')


# MAIN CODE

if __name__ == '__main__':

    args = load_inputs()

    # RUN MAIN PIPELINE

    main_start = time.time()
    
    if args.inseqcon is None:
        seq_con = all_against_all(args.infile, mode = args.mode, tmp_folder = args.tmp_folder, max_evalue = args.max_evalue, infasta = args.infasta)
    else:
        seq_con = args.inseqcon

    if 'NC' in args.mode:
        seq_con = neighborhood_correlation(seq_con, min_nc = args.min_nc)

    if args.plddt_weighted:
        seq_con = weight_by_plddt(args.infile, seq_con)

    if args.find_clusters:
        clusters = extract_clusters(seq_con, evalue_clustering = args.evalue_clustering, min_size = args.cluster_size, sigma = args.evalue_sigma, n_iterations = args.n_iterations)
    else:
        clusters = None

    out_net  = write_clans_file(args.infile, seq_con, infasta = args.infasta, mode = args.mode, weighted=args.plddt_weighted, evalue_clustering=args.max_evalue, clusters = clusters)

    main_end = time.time()
    main_numb_seconds = main_end - main_start
    print("\n#### Finished job after: {} \n".format(time.strftime('%H hours %M min %S sec', time.gmtime(main_numb_seconds))))

