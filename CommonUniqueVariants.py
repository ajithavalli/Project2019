import vcf
import argparse
from collections import defaultdict
import pandas as pd

# example run:
# python CommonUniqueVariants.py -v vcf1.vcf -v vcf2.vcf -v vcf3.vcf -v vcf4.vcf -c cosmic.vcf 

vcf_variants = defaultdict(set)
vcf_positions = defaultdict(set)
variant_features = defaultdict(dict)

def parseVCF(vcffile):
	vcfreader = vcf.Reader(open(vcffile))
	samples = vcfreader.samples
	sampleName = None
	features = ['AF', 'AC', 'DP', 'ExcessHet','FS','MQ','QD',
	'GQ', 'DP', 'QUAL']

	if len(samples) == 0:
		print("No sample found, consolidating it as known_mutations")
		sampleName = vcffile.split("/")[-1].replace(".vcf", "") + "_known_mutations"

	for record in vcfreader:
		alts = record.ALT
		chrm = record.CHROM.replace('chr', '')
		chrm = chrm.replace('Chr', '')
		if sampleName and "_known_mutations" in sampleName:
				vcf_positions[sampleName].add(chrm + "_" + str(record.POS))
				alts = record.ALT
				for alt in alts:
					keyh = "_".join([chrm, str(record.POS), record.REF, str(alt)])
					vcf_variants[sampleName].add(keyh)
					for field in features:
						if field in ('AC', 'AF'):
							variant_features[sampleName][keyh][field] = record.INFO[field][0]
						elif field in ('DP', 'ExcessHet','FS','MQ','QD'):
						    variant_features[sampleName][keyh][field] = record.INFO[field]
						elif field in ('QUAL'):
							variant_features[sampleName][keyh][field] = record.__getattribute__(field)
						else:
							print(field + " not found for " + keyh)
		else:
			for sample in samples:
				sampleName = vcffile.split("/")[-1].replace(".vcf", "") + "_" + sample
				vcf_positions[sampleName].add(chrm + "_" + str(record.POS))
				try:
				    alts = record.genotype(sample).gt_bases.split("|")
				except:
					continue
				if len(alts) == 1:
					alts = record.genotype(sample).gt_bases.split("/")
				for alt in alts:
					keyh = "_".join([chrm, str(record.POS), record.REF, str(alt)])
					vcf_variants[sampleName].add(keyh)
					variant_features[sampleName][keyh] = {}
					for field in features:
						if field in ('AF', 'AC'):
							variant_features[sampleName][keyh][field] = record.INFO[field][0]
						elif field in ('ExcessHet','FS','MQ','QD') and field in record.INFO:
						    variant_features[sampleName][keyh][field] = record.INFO[field]
						elif field in ('GQ', 'DP'):
						    variant_features[sampleName][keyh][field] = record.genotype(sample)[field]
						elif field in ('QUAL'):
							variant_features[sampleName][keyh][field] = record.__getattribute__(field)
						else:
							print(field + " not found for " + keyh)


def parseCosmic(cosmicVCF):
	vcfreader = vcf.Reader(open(cosmicVCF))
	vcfreader._parse_info
	cosmicMutations = set()
	cosmicMutationScores = dict()
	for record in vcfreader: 
		for alt in record.ALT:
			chrm = record.CHROM.replace('chr', '')
			chrm = chrm.replace('Chr', '')
			identifier = "_".join([chrm, str(record.POS), record.REF, str(alt)])
			cosmicMutations.add(identifier)
			# Fill in fathmm score consolidation here
			# cosmicMutationScores[identifier] = vcfreader._parse_info(record.INFO) 
	return cosmicMutations, cosmicMutationScores


def common_mutations(vcf_positions, vcf_variants, cosmicMutations):
	samples = vcf_positions.keys()
	print("SAMPLES being compared:", samples)
	
	allCommon = set()
	commonVariants = {}
	cosmicOverlap = {}

	## Common across all samples:
	for sample in samples:
		if allCommon == set():
			allCommon = vcf_variants[sample]
		else:
			allCommon = allCommon.intersection(vcf_variants[sample])
	print("Common Across all samples:", allCommon)

	## Common variants between every 2 samples and with cosmic:
	for sample1 in samples:
		commonVariants[sample1] = {}
		for sample2 in samples:
			if sample1 == sample2:
				continue
			else:
				commonVariants[sample1][sample2] = len(vcf_variants[sample1].intersection(vcf_variants[sample2]))
		cosmicOverlap[sample1] = len(vcf_variants[sample1].intersection(cosmicMutations))

	print("COMMON VARIANTS matix:")
	print(pd.DataFrame(commonVariants).to_string())
	print("COSMIC OVERLAP:", cosmicOverlap)

	fout = open("AllCommon_PCA.txt", 'w')
	line = ["Sample"]
	features = ['AF', 'AC', 'DP', 'ExcessHet','FS','MQ','QD',
	'GQ', 'DP', 'QUAL']
	for variant in allCommon: 
		for field in features:
			line.append(variant + "_" + field)
	fout.writelines("\t".join(line) + "\n")


	for sample in samples:
		line = [sample]
		for variant in allCommon:
			for field in features:
				if field in variant_features[sample][variant]:
					line.append(str(variant_features[sample][variant][field]))
				else:
					line.append('0')
		fout.writelines("\t".join(line) + "\n")
	fout.close()


	return (commonVariants, cosmicOverlap)

def unique_mutations(vcf_variants):
	samples = vcf_variants.keys()
	vcf_exclusives = {}

	## Unique variants in every sample
	for sample1 in samples:
		vcf_exclusives[sample1] = vcf_variants[sample1]
		for sample2 in samples:
			if sample1 == sample2:
				continue
			else:
				vcf_exclusives[sample1] = vcf_exclusives[sample1] - vcf_variants[sample2]

		fout = open(sample1 + "_uniqueVariants.txt", 'w')
		for variant in vcf_exclusives[sample1]:
			fout.writelines(variant + "\n")
		fout.close()

	return vcf_exclusives

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-v", "--vcfs", action="append", help="list vcfs")
    parser.add_argument("-c", "--cosmic", help="cosmic vcf")
    parser.add_argument("-o", "--output_path", help="output path")

    args = parser.parse_args()
    for eachVcf in args.vcfs:
    	print("PARSING", eachVcf)
    	parseVCF(eachVcf)

    cosmicMutations, cosmicMutationScores = parseCosmic(args.cosmic)
    (commonVariants, cosmicOverlap) = common_mutations(vcf_positions, vcf_variants, cosmicMutations)
    #vcf_exclusives = unique_mutations(vcf_variants)

if __name__ == "__main__":
    main()
