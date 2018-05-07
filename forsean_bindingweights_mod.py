#!/usr/bin/python

"""
Shilpa Kobren (primary author) and Adam Berman
Princeton University
7 May 2018

Find binding sites in human proteins from weight vector files
"""

import os
import gzip
import pickle
import json


def save_obj(obj, name):
    with open('obj/'+ name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def load_obj(name):
    with open('obj/' + name + '.pkl', 'rb') as f:
        return pickle.load(f)


def interaction_track(track_name):
  """
  :param track_name: string corresponding to the name of a track (in my weight vector files)
  :return: boolean indicating whether the track corresponds to binding site information (True) or not (False)
  """

  if 'Homology_binding' not in track_name:
    current_tracks = []
    for bd_track in track_name.split(','):
      if not bd_track.startswith('PF') or bd_track.endswith(':complete'):
        continue
      current_tracks.append(bd_track)
    if len(current_tracks) < 1:
      return False  # consider only interaction-based or within-domain weight based tracks
  return True


def get_indices_from_interval(intervals_str):
  """
  :param intervals_str: comma-separated intervals (e.g., "2-10,20-45,90-100")
  :return: set of all indices as specified by the intervals string
  """

  current_indices = set()

  for current_interval in intervals_str.split(','):
    start_index, end_index = map(int, current_interval.split('-')[:2])
    for current_index in xrange(start_index, end_index + 1):
      current_indices.add(current_index)

  return current_indices


genes_and_bindingweights = {}

# Find binding sites in human proteins from weight vector files
if __name__ == "__main__":
  
  path_to_bindingfiles = '/Users/adamberman/metabolite-scorer/Aggregate/'

  # Iterate through file hierarchy to get all binding weight data for all genes in all chromosomes
  for chromosome in os.listdir(path_to_bindingfiles):
    if not chromosome.startswith('.') and os.path.isdir(path_to_bindingfiles+chromosome+'/'):
      print "CHROMOSOME: " + chromosome ####
      for gene_id in os.listdir(path_to_bindingfiles+chromosome+'/'):
        if not gene_id.startswith('.') and os.path.isdir(path_to_bindingfiles+chromosome+'/'+gene_id+'/'):
          print "GENE ID: " + gene_id ####

          # keep track of the max functional weight / score at each 0-index
          max_func_weights = {}

          for example_weightvector in os.listdir(path_to_bindingfiles+chromosome+'/'+gene_id+'/'):
            if not example_weightvector.startswith('.'):
              print "WEIGHT VECTOR: " + example_weightvector
              
              # Open file
              wv_handle = gzip.open(path_to_bindingfiles+chromosome+'/'+gene_id+'/'+example_weightvector)

              for wv_line in wv_handle:
                if wv_line.startswith('#') or len(wv_line.split('\t')) < 5:
                  continue

                # Collect data
                prot_id, _, track_name, interval, func_weights = wv_line[:-1].split('\t')[:5]
                
                # Store maximum binding weight at each position across all tracks
                if interaction_track(track_name):
                  # which protein indices were we even able to model
                  prot_positions = get_indices_from_interval(interval)  
                  # dictionary of index -> binding value
                  func_weights = {int(v.split(':')[0]): float(v.split(':')[1]) for v in func_weights.split(',')} 
                  for index, fweight in func_weights.items():
                    if index in max_func_weights:
                      max_func_weights[index] = max(fweight, max_func_weights[index])
                    else:
                      max_func_weights[index] = fweight

          print "MAX WEIGHTS:"
          print max_func_weights

          # Store genes as corresponding maximum binding weight by position dictionary in overall dictionary
          genes_and_bindingweights[gene_id] = max_func_weights

          """
          NOTE: track_name will tell you if the binding scores came from homology modeling (Homology_binding) or a domain (otherwise). 
                The span of the homology/domain match is included, and the ligand being bound to is also specified. 
          """

# Store genes_and_bindingweights dictionary in a JSON file for later use
with open('genes_and_bindingweights.json', 'w') as fp:
    json.dump(genes_and_bindingweights, fp)
