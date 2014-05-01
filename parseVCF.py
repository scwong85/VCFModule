import sys

#usage python parseVCF.py vcffile outputfile
#example python parseVCF.py myvcf.txt output.txt
#output explanatin:
#column delimited file
#SCAFFOLD_ID	POSITION	BASE*	DEPTH	BASECOUNT(A C G T)***
#*Top 2 most common bases will be encapsulated in []
#**Repetition of A, C, G, T base count for N individuals from pileup file

class ParseVCF:
	
	vcf_file = ''
	output_file = ''
	chisq_dict = {1: 3.841, 
				  2: 5.991,
				  3: 7.815,
				  4: 9.488,
				  5: 11.07,
				  6: 12.59,
				  7: 14.07,
				  8: 15.51,
				  9: 16.92}
	
	#constructor
	def __init__(self, vcf_file, output_file):
		self.vcf_file = vcf_file
		self.output_file = output_file

	#factorial
	def factorial(self, n):
		if n == 1:
			return n
		else:
			return n * self.factorial(n-1)
	
	#return total genotype combination
	def combinationWithRepeat(self, alleleCount):
		return self.factorial(alleleCount + 2 - 1) / (self.factorial(2)*self.factorial(alleleCount -1))
		
	#get unique value of a list	
	def uniqueListLength(self, mylist):
		tmp = []
		for item in mylist:
			if item not in tmp:
				tmp.append(item)
		return tmp, len(tmp)
	
	#return maf from an allele dict
	def minorAlleleFrequency(self, mydict):	
		maf = 1
		allele = ''
		for key in mydict.keys():
			item = mydict[key]
			if float(item) <= maf:
				maf = float(item)
				allele = key
		return allele, maf
				
	#count number of allele in genotype
	def countAllele(self, allele, genotype):
		count = 0
		for i in genotype:
			if str(allele) == str(i):
				count += 1
		return count	
	
	#return genotype set
	def genoSet(self, alleleCount):
		dict = {}
		for i in range(0, alleleCount):
			for j in range(0, alleleCount):
				key = '%s%s' %(i, j)
				dict[key] = 0
		return dict	
		
	#return genotype position in the genotype likelihood list
	def getGenotypePosition(self, genocode):
		a1 = int(genocode[0])	#allele 1
		a2 = int(genocode[-1])	#allele 2
		return (a2*(a2+1)/2)+a1	#position for a1/a2
			
	#read and process vcf file, and write to output file	
	def parseVcf(self, line):
		line = line.strip().split('\t')
		allele = {}
		chrom = line[0]
		pos = line[1]
		ref = line[3]				
		alt = line[4]				
		qual = line[5]
		genoInfo = line[8]				#genotype info GT:PL:GQ, order might varied
		ind = line[9:]					#individual genotype call
		alt = alt.split(',')			#alternative allele			
		total_alt = len(alt)			#total alternative allele
		
		#fill allele dict
		allele[0] = ref
		for i in range(0, total_alt):
			allele[i+1] = alt[i] 
		###
		
		#determining order of genotype info
		genoInfo = genoInfo.split(':')
		GTi = -1
		PLi = -1
		GQi = -1
		for i in range(len(genoInfo)):
			if genoInfo[i] == 'GT':
				GTi = i
			elif genoInfo[i] == 'PL':
				PLi = i
			elif genoInfo[i] == 'GQ':
				GQi = i
		###
				
		geno_dict = self.genoSet(total_alt + 1)		#return dictionary with total combination of genotype		
		total = 0									#total reliable genotype
	
		#calculate genotype information for each sample		
		for sample in ind:
			sample = sample.split(':')
			geno = sample[GTi]
			geno = geno.replace('/','').replace('|','')
			
			#extracting likelihood information
			likelihood = sample[PLi]
			geno_qual = sample[GQi]
			geno_L = likelihood.split(',')			#likelihood of each genotype
			total_genotype = len(geno_L)			#total genotypes
			l_list = []								#list to store likelihood
			
			#get likelihood for each genotype
			#lh_geno = {}
			#for genotype in geno_dict.keys():		#key would be 00, 01, 10, 11, 12, etc... 0 = ref, 1,2,3... = alt
				#likelihood_pos = self.getGenotypePosition(genotype)
				#lh_geno[genotype] = geno_L[likelihood_pos]
			
			for l_value in geno_L:
				l_list.append(int(l_value))			
			l_list.sort()
			###
						
			#check likelihood difference between the most likely and the second most likely
			min_diff = l_list[1] - l_list[0]
			###
			
			#if there is likelihood and the difference > 10
			unique_likelihood, likelihood_length = self.uniqueListLength(geno_L)
			#if len(geno_L) > 3:
				#print len(geno_L), geno_L, likelihood_length, min_diff, l_list
			
			if likelihood_length > 1  and min_diff > 10:
				geno_dict[geno] = geno_dict[geno] + 1
				total += 1
				#if likelihood_length > 3:
				#	print geno_dict
			###
			###end of processing one individual
		###end of processing one line
			
		#calculate allele frequency
		afreq = {}
		for i in allele.keys():
			afreq[allele[i]] = 0
			allele_count = 0
			total_count = 0
			for j in geno_dict.keys():		#j = 00, 01, 10, 11, 02, etc...
				if str(i) in j:				#if allele i was found in genoype j
					count_i = self.countAllele(i, j)
					allele_count += count_i*(int(geno_dict[j]))
				total_count += int(geno_dict[j])
			if total_count > 0:
				afreq[allele[i]] = float(allele_count) / (2*total_count) 
		#print afreq
		###
				
		#calculate HWE exact test
		#formula = p2 + q2 + r2 + pq + qp + pr + rp + qr + rq...
		chisq = 0
		if total > 0:
			used = []
			for genotype in geno_dict.keys():	#00, 01, 10, 11, etc...
				a1 = int(genotype[0])
				a2 = int(genotype[1])
				if afreq[allele[a1]] != 0 and afreq[allele[a2]] != 0:
					if a1 == a2:					#homo
						oh = geno_dict[genotype]	#observed homo 
						eh = afreq[allele[a1]]*afreq[allele[a1]]*total	#expected = p*p*n
						chisq += (((oh - eh)*(oh - eh))/(eh))
					else:							#hetero
						if genotype not in used:
							g1 = '%s%s' %(a1, a2)
							g2 = '%s%s' %(a2, a1)
							oh = geno_dict[g1] + geno_dict[g2]
							eh = 2*afreq[allele[a1]]*afreq[allele[a2]]*total
							chisq += (((oh - eh)*(oh - eh))/(eh))
							used.append(g1)
							used.append(g2)
				else:
					chisq = 'NA'
					break
		else:
			chisq = 'NA'
		###
		
		if chisq != 'NA' and total_alt <= 10:
			threshold = self.chisq_dict[total_alt]
			if chisq <= threshold:
				inHWE = 'YES'
			else:
				inHWE = 'NO'
		else:
			inHWE = 'NA'
		
			
		#transforming geno_dict key
		geno_dictA = {}	
		for key in geno_dict.keys():
			c1 = int(key[0])
			c2 = int(key[1])
			try:		
				a1 = allele[c1]
				a2 = allele[c2]
				akey = '%s%s' %(a1,a2)
				geno_dictA[akey] = geno_dict[key]
			except KeyError:
				pass
			
		#get minor allele frequency
		minorAllele, maf = self.minorAlleleFrequency(afreq)
		mafString = '%s: %s' %(minorAllele, maf)		
				
		string = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' %(chrom, pos, ref, ', '.join(alt), qual, geno_dictA, len(ind), float(total)/len(ind), mafString, chisq, inHWE)
		return string
	
	#read vcf file
	def read_vcf(self):
		f = open(self.vcf_file, 'r')
		w = open(self.output_file, 'w')
		header = 'Chromosome\tPosition\tRef\tAlt\tQuality\tGenotypeCount\tSampleSize\tCallRate\tMAF\tHWEChiSq\tinHWE?\n'
		w.write(header)
		for line in f:
			if '#' not in line:
				results = 1
				p_string = self.parseVcf(line)
				#write output
				if results*1 == 1:	
					w.write(p_string)
		f.closed
		w.flush()
		w.closed
				
	#run	
	def run(self):
		self.read_vcf()

	
p = ParseVCF(sys.argv[1], sys.argv[2])
p.run()