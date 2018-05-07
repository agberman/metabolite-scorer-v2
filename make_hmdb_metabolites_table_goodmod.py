# Adam Berman
# Princeton University
# 5 May 2017


import os
import xmltodict
import xml.etree.ElementTree as ET
import itertools


# Parse metabolites
all_metabolites = []

# Iterate over all metabolites
# Folder containing all metabolite's xml files (named "hmdb_metabolites") can be found 
# at the following link: http://www.hmdb.ca/system/downloads/current/hmdb_metabolites.zip
for filename in os.listdir('/Users/adamberman/Independent_Work_v2/hmdb_metabolites'):

	if filename.startswith('HMDB'):

		# Parse the metabolite
		root = ET.parse('/Users/adamberman/Independent_Work_v2/hmdb_metabolites/' + filename).getroot()
		origin = root.find('ontology').find('origins').find('origin')

		# Keep only endogenous metabolites
		if (origin != None) and (origin.text == 'Endogenous'):

			# Collect all genes associated with the given metabolite
			gene_names = []
			for protein in root.find('protein_associations').findall('protein'):
				gene_name = protein.find('gene_name').text
				if gene_name != None:
					gene_names.append(gene_name.strip())

			# Only continue if the metabolite has at least three associated genes
			#if (len(gene_names) >= 3):

			# Only continue if the metabolite has at least one associated gene:
			if gene_names:

				# Collect the name and KEGG ID of the metabolite
				name = root.find('name').text.strip()
				kegg = root.find('kegg_id').text
				if kegg == None:
					kegg = "None"
				else:
					kegg = kegg.strip()

				# Collect the synonyms of the metabolite
				synonyms = []
				for synonym in root.find('synonyms').findall('synonym'):
					synonyms.append(synonym.text.strip())

				# Collect the SMILES identifier of the metabolite
				smiles = root.find('smiles').text
				if smiles != None:
					smiles = smiles.strip()
				
				# InChI instead of SMILES
				#if smiles == None:
					
				inchi = root.find('inchi').text
				if inchi != None:
					inchi = inchi.strip()
					inchi_key = root.find('inchikey').text
					if inchi_key != None:
						inchi_key = inchi_key.strip()

				# Add new metabolite to all_metabolites
				metabolite = {'name': name, 'synonyms': synonyms, 'kegg': kegg, 'gene_names': gene_names, 'smiles': smiles, 'inchi': inchi, 'inchi_key': inchi_key, 'filename': filename}
				all_metabolites.append(metabolite)


print 'Number of endogenous metabolites: '
print len(all_metabolites)
#20,440 out of 61,388 metabolites that are endogenous and have at least one associated gene
#'''

# Simple test version of all_metabolites
'''
all_metabolites = [{'name': 'dog', 'synonyms': ['canine', 'pupper', 'poochy'], 'kegg': 'kegg1', 'gene_names': ['1AA', 'B1B'], 'smiles': 'IOU1K9_1', 'inchi': 'inchi1', 'inchi_key': 'inchi_key1', 'filename': 'Dogger_1.xml'}, 
					{'name': 'hound', 'synonyms': ['dog', 'pupper', 'doggo'], 'kegg': 'kegg2', 'gene_names': ['B2B', 'BBB'], 'smiles': 'IOU1K9_2', 'inchi': 'inchi1', 'inchi_key': 'inchi_key2', 'filename': 'Dogger_2.xml'}, 
					{'name': 'houndy', 'synonyms': ['new', 'newer', 'newest'], 'kegg': 'kegg3', 'gene_names': ['EEE', 'FFF'], 'smiles': 'IOU1K9_3', 'inchi': 'inchi3', 'inchi_key': 'inchi_key3', 'filename': 'Dogger_3.xml'},
					{'name': 'New', 'synonyms': ['Title', 'hound', 'poochy'], 'kegg': 'kegg4', 'gene_names': ['GGG', 'HHH'], 'smiles': None, 'inchi': 'inchi4', 'inchi_key': 'inchi_key4', 'filename': 'Dogger_4.xml'},
					{'name': 'doggy', 'synonyms': ['canine', 'pupper', 'poochy'], 'kegg': 'kegg1', 'gene_names': ['YYY', 'ZZZ'], 'smiles': 'IOU1K9_1', 'inchi': 'inchi10', 'inchi_key': 'inchi_key1', 'filename': 'Dogger_1.xml'}]
'''




# Remove synonymous metabolites from all_metabolites
smiles_dictionary = {}
for metabolite in all_metabolites:

	genes = metabolite.get('gene_names')
	smiles = metabolite.get('smiles')
	kgg = metabolite.get('kegg')
	syn = metabolite.get('synonyms')
	name = metabolite.get('name')
	inchi = metabolite.get('inchi')
	inchi_key = metabolite.get('inchi_key')
	name_and_syn = syn + [name]

	merged = False
	for key, value in smiles_dictionary.iteritems():

		# Merge metabolite with existing smiles_dictionary entry if both metabolites have the same set of associated genes
		if (set(genes) == set(value.get('gene_names'))):

			# Perform merge operation
			new_syn = list(set(name_and_syn + value.get('synonyms')))
			if key in new_syn:
				new_syn.remove(key)

			new_gn = list(set(metabolite.get('gene_names') + value.get('gene_names')))
			new_gn.sort()

			# Create new entry and store it into the dictionary
			new_value = {'synonyms': new_syn, 'kegg': value.get('kegg'), 'gene_names': new_gn, 
							'smiles': value.get('smiles'), 'inchi': value.get('inchi'), 
							'inchi_key': value.get('inchi_key'), 'filename': value.get('filename')}

			smiles_dictionary[key] = new_value

			# Denote that merge occured
			merged = True
			break

	# If the metabolite was not merged, create new smiles_dictionary entry for metabolite
	if not merged:
		smiles_dictionary[name] = {'synonyms': syn, 'kegg': metabolite.get('kegg'), 
									'gene_names': metabolite.get('gene_names'), 
									'smiles': metabolite.get('smiles'), 
									'inchi': metabolite.get('inchi'), 
									'inchi_key': metabolite.get('inchi_key'), 
									'filename': metabolite.get('filename')}

print 'Number of endogenous metabolites after merging synonyms: '
print len(smiles_dictionary)


# Print SMILES metabolites into table
print 'Metabolite' + '\t' + 'SMILES' + '\t' + 'KEGG' + '\t' + 'Gene Name'
for key, value in smiles_dictionary.iteritems():
	n = key
	if n == None:
		n = 'None'
	s = value.get('smiles')
	if s == None:
		s = 'None'
	k = value.get('kegg')
	if k == None:
		k = 'None'
	gns = value.get('gene_names')
	if gns == None:
		gns = 'None'

	for gn in gns:
		print n + '\t' + s + '\t' + k + '\t' + gn







# TEST CODE

'''
# Simple test version of all_metabolites
all_metabolites = [{'name': 'dog', 'synonyms': ['canine', 'pupper', 'poochy'], 'gene_names': ['1AA', 'B1B'], 'smiles': 'IOU1K9_1', 'filename': 'Dogger_1.xml'}, 
					{'name': 'hound', 'synonyms': ['dog', 'pupper', 'doggo'], 'gene_names': ['B2B', 'BBB'], 'smiles': 'IOU1K9_2', 'filename': 'Dogger_2.xml'}, 
					{'name': 'hound', 'synonyms': ['new', 'newer', 'newest'], 'gene_names': ['EEE', 'FFF'], 'smiles': 'IOU1K9_3', 'filename': 'Dogger_3.xml'},
					{'name': 'New', 'synonyms': ['Title', 'hound', 'poochy'], 'gene_names': ['GGG', 'HHH'], 'smiles': 'IOU1K9_4', 'filename': 'Dogger_4.xml'}]
'''

'''
# Simple test version of all_inchi_metabolites
all_inchi_metabolites = [{'name': 'dog', 'synonyms': ['canine', 'pupper', 'poochy'], 'gene_names': ['1AA', 'B1B'], 'inchi': 'inchi_1', 'inchi_key': 'inchi_key_1', 'filename': 'Dogger_1.xml'}, 
					{'name': 'hound', 'synonyms': ['dog', 'pupper', 'doggo'], 'gene_names': ['B2B', 'BBB'], 'inchi': 'inchi_2', 'inchi_key': 'inchi_key_2', 'filename': 'Dogger_2.xml'}, 
					{'name': 'hound', 'synonyms': ['new', 'newer', 'newest'], 'gene_names': ['EEE', 'FFF'], 'inchi': 'inchi_3', 'inchi_key': 'inchi_key_3', 'filename': 'Dogger_3.xml'},
					{'name': 'New', 'synonyms': ['Title', 'hound', 'poochy'], 'gene_names': ['GGG', 'HHH'], 'inchi': 'inchi_4', 'inchi_key': 'inchi_key_4', 'filename': 'Dogger_4.xml'}]
'''

'''
all_metabolites = [{'name': 'dog', 'synonyms': ['canine', 'pupper', 'poochy'], 'kegg': 'dog_kegg', 'gene_names': ['1AA', 'B1B'], 'smiles': 'IOU1K9_1', 'filename': 'Dogger_1.xml'}, 
					{'name': 'hound', 'synonyms': ['dog', 'pupper', 'doggo'], 'kegg': 'hound_kegg', 'gene_names': ['B2B', 'BBB'], 'smiles': 'IOU1K9_2', 'filename': 'Dogger_2.xml'}, 
					{'name': 'Title', 'synonyms': ['new', 'newer', 'newest'], 'kegg': 'Title_kegg', 'gene_names': ['EEE', 'FFF'], 'smiles': 'IOU1K9_3', 'filename': 'Dogger_3.xml'},
					{'name': 'New', 'synonyms': ['Title', 'hound', 'poochy'], 'kegg': 'New_kegg', 'gene_names': ['GGG', 'HHH'], 'smiles': 'IOU1K9_4', 'filename': 'Dogger_4.xml'}]
'''


'''
all_metabolites = [{'name': 'test1', 'synonyms': ['a', 'b', 'c'], 'kegg': 'same_kegg', 'gene_names': ['1AA', 'B1B'], 'smiles': 'IOU1K9_1', 'filename': 'Dogger_1.xml'}, 
					{'name': 'test2', 'synonyms': ['d', 'e', 'f'], 'kegg': 'diff_kegg', 'gene_names': ['B2B', 'BBB'], 'smiles': 'IOU1K9_2', 'filename': 'Dogger_2.xml'}, 
					{'name': 'test3', 'synonyms': ['g', 'h', 'i'], 'kegg': 'same_kegg', 'gene_names': ['EEE', 'FFF'], 'smiles': 'IOU1K9_3', 'filename': 'Dogger_3.xml'},
					{'name': 'test4', 'synonyms': ['j', 'k', 'l'], 'kegg': 'same_kegg', 'gene_names': ['GGG', 'HHH'], 'smiles': 'IOU1K9_4', 'filename': 'Dogger_4.xml'}]
'''

# GRAVEYARD
# Merge metabolite with existing smiles_dictionary entry if both metabolites have the same set of associated genes, or either their SMILES or INCHI strings are the same
		#if (set(genes) == set(value.get('gene_names'))) or (((inchi != None) and (inchi == value.get('inchi'))) or ((smiles != None) and (smiles == value.get('smiles')))):
		#if ((inchi == value.get('inchi')) or (smiles == value.get('smiles'))):

		#if ((key in name_and_syn) 
		#	or (not set(value.get('synonyms')).isdisjoint(name_and_syn)) 
		#	or ((kgg != "None") and (kgg == value.get('kegg'))) 
		#	or (smiles == value.get('smiles'))):

#print '{0:40s} {1:500s} {2:30s}'.format(n, s, gn)
