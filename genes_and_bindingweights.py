# Adam Berman
# Princeton University
# 7 May 2018

import json
import csv

# Import Shilpa's BioLiP genes and binding weights
with open('genes_and_bindingweights.json', 'r') as fp:
    genes_and_bindingweights = json.load(fp)

# Clean binding weight data
genes_and_bindingweights = {str(key):value for key,value in genes_and_bindingweights.items()}
for gene, sub_dict in genes_and_bindingweights.iteritems():
	genes_and_bindingweights[gene] = {int(key):value for key,value in genes_and_bindingweights[gene].items()}
genes_and_bindingweights = {key.split('.', 1)[0]:value for key,value in genes_and_bindingweights.items()}



# Create mutations and positions dictionary
mutations_and_positions = {}

ensemble = ''
start_position = ''
end_position = ''

# Import list of 62,315 missense mutations with positional information
with open('mutations_with_positions.csv', 'rb') as csvfile:
	reader = csv.reader(csvfile, delimiter = '	')
	next(reader)
	for row in reader:
		for element in row:
			element_split = element.replace('"', '').split(',')
			ensemble = element_split[1]
			start_position = int(element_split[2])-1 # Correct for 1-indexed TCGA positions
			end_position = int(element_split[3])-1 # Correct for 1-indexed TCGA positions

		if ensemble not in mutations_and_positions:
			mutations_and_positions[ensemble] = [[start_position, end_position]]
		else:
			mutations_and_positions[ensemble].append([start_position, end_position])



# Look up each mutation in each gene in mutations_and_positions in Shilpa's data (genes_and_bindingweights), 
# annotating mutations whose positional range includes at least one binding weight greater than or equal to 
# 0.1 with a 1 (a positionally "significant" mutation), and all other mutations with a 0.

# Initialize dictionary that counts the number of positionally significant mutations per gene
genes_with_num_sig_mutations = {}
for gene in mutations_and_positions.iterkeys():
	# Remove genes that do not appear in Shilpa's positional data at all
	if gene in genes_and_bindingweights:
		genes_with_num_sig_mutations[gene] = 0

# 16,097 of the 16,133 genes in TCGA's mutations_and_positions appear in Shilpa/Biolip's genes_and_bindingweights

# Check each position within between the start position and end position of each mutation;
# if any position has a binding weight score >= 0.1, 
for gene, mutation_positions in mutations_and_positions.iteritems():
	if gene in genes_and_bindingweights:
		for mutation_position in mutation_positions:
			max_bindingweight = 0
			for pos in range(mutation_position[0], mutation_position[1]+1):
				if pos in genes_and_bindingweights[gene]:
					if genes_and_bindingweights[gene][pos] > max_bindingweight:
						max_bindingweight = genes_and_bindingweights[gene][pos]
			if max_bindingweight >= 0.1:
				genes_with_num_sig_mutations[gene] += 1


# Print out genes_with_num_sig_mutations to a CSV file to be read into R
with open('genes_with_num_sig_mutations.csv', 'wb') as csv_file:
    writer = csv.writer(csv_file)
    for key, value in genes_with_num_sig_mutations.iteritems():
       writer.writerow([key, value])
