configfile: "config.json"
FOLDER=config['folder_H'] # this option is only for testing on different machines, originaly it can be only one run_folder in config file

#check if paired end
reads=['1']
paired_end=config['paired_end'] #Miseq test data is paired end!! others single read
if paired_end=='True':
	reads=['1','2']
print('paired end: ',paired_end)
print(config['machine'])
SAMPLES=config['samples_H'] #CHANGE SAMPLE SHEET IF CHANGING CUSTOM <->ILLUMINA BARCODES for testing
n_lanes=config['n_lanes']
custom_lanes=config['custom_lanes']
no_align=config['no_align']



#make list of lists as lanes with samples in lanes
presamps = [[] for x in range(len(n_lanes))]

for n in n_lanes:
	for sample in SAMPLES: 
		if sample[4] == n:
			presamps[n_lanes.index(n)].append(sample)



#generating output of bcl2fastq, samples full names in list of lists
samps_in_lanes=[]
for l in presamps:
	if config['machine']=='M' or config['machine']=='H':			
		l=[(sample+'_S'+str(n)+'_L00'+sample[4]+'_R') for n,sample in enumerate(l,1)]
		l=[sample+read+'_001.fastq.gz' for sample in l for read in reads]
		samps_in_lanes.append(l)
	elif config['machine']=='N':
		l=[(sample+'_S'+str(n)+'_L001') for n,sample in enumerate(l,1)]
		l=[sample+'_R'+read+'_001.fastq.gz' for sample in l for read in reads]
		samps_in_lanes.append(l)


'''
#Dictioary with lanes and samples where for lane with custom barcode there's only one sample. 
#In case je demultiplex will be in seperate rule this code can be useful. 

samps_in_lanes_D={}
for l in n_lanes:	
		samps_custom_bcl2fastq=[]
		if l not in custom_lanes:
			samps_in_lanes_D[l]=samps_in_lanes[n_lanes.index(l)]
			#print(l)
		else:			
			for r in reads:				
				if config['machine']=='N':				
					samp='lane'+l+'_S1'+'_R'+r+'_001.fastq.gz'
				else:
					samp='lane'+l+'_S1'+'_L00'+l+'_R'+r+'_001.fastq.gz'
				samps_custom_bcl2fastq.append(samp)				
				samps_in_lanes_D[l]=samps_custom_bcl2fastq
print(samps_in_lanes_D)
'''

#Dictionary with all samples assigned to particular lane, for agregate_input checkpoint.bcl2fastq.get 
samps_in_lanes_D_all={}
for l in n_lanes:	
	samps_in_lanes_D_all[l]=samps_in_lanes[n_lanes.index(l)]


#list of samples in Aligned_lane* folders
aligned_samps_in_lanes=[]
for l_aligned in presamps:	
	l_aligned=[sample+'.bam' for sample in l_aligned]
	aligned_samps_in_lanes.append(l_aligned)


samps_aligned_D={}
for l in n_lanes:	
	if l not in no_align:
		samps_aligned_D[l]=aligned_samps_in_lanes[n_lanes.index(l)]



all_samples = [FOLDER+'/Aligned_lane'+lane+'/'+sample for lane in samps_aligned_D for sample in samps_aligned_D[lane]] 


#dictionary of samples which shouldn't be aligned, only bcl2fastq should be done for them
not_aligned_samps_D={}
for l in n_lanes:
	if l in no_align:		
		not_aligned_samps_D[l]=aligned_samps_in_lanes[n_lanes.index(l)]
#below 3 lines add samples to input in 'all' rule, for which only bcl2fastq should be done. Their name is consistent with get_sample_wc rule output. 		
for lane in not_aligned_samps_D:
	for sample in not_aligned_samps_D[lane]:
		all_samples.append(FOLDER+'/fake/get_WC_fake/Aligned_lane'+lane+'/'+sample)



#below function produces files with barcodes and names for je demultiplex script in case of custom lane 
def barcodes_to_file_for_dpx():
	BarCdtxt = FOLDER+'/SampleSheetOriginal.csv'
	for lane in custom_lanes:
		t_list=[]
		with open(BarCdtxt, 'r') as f:
			for line in f:			
				if line.split(',')[2][4] == lane:
					s=line.split(',')[2]
					b=line.split(',')[4]	
					if len(reads)==2:
						for s_n in samps_in_lanes[int(lane)-1]:
							if s_n.split('_')[0]==s and s_n.split('_')[3]=='R1':
								n_r1=s_n
							elif s_n.split('_')[0]==s and s_n.split('_')[3]=='R2':
								n_r2=s_n	
						t_list.append((s,b[:int(len(b)/2)]+':'+b[int(len(b)/2):],n_r1,n_r2))				
					else:
						for s_n in samps_in_lanes[int(lane)-1]:
							if s_n.split('_')[0]==s:
								n=s_n
						t_list.append((s,b,n))
		filename=FOLDER+'/log_pipeline/lane'+lane+'_barcodes.txt'		
		if not os.path.exists(os.path.dirname(filename)):
			os.makedirs(os.path.dirname(filename))
		with open(filename, 'w') as fp:
			if len(reads)==2: 
				fp.write('\n'.join('{} {} {} {}'.format(x[0],x[1],x[2],x[3]) for x in t_list))
			else:
				fp.write('\n'.join('{} {} {}'.format(x[0],x[1],x[2]) for x in t_list))

barcodes_to_file_for_dpx()


def bcl2fastq_threads_setting(wildcards):
    """Set up threads for the four bcl2fastq stages
    Illumina defaults:
    - 4 threads for reading the data:    -r, --loading-thread
    - 4 threads for writing the data:    -w, --writing-threads
    - 20% for demultiplexing data:       -d, --demultiplexing-threads
    - 100% for processing demultiplexed: -p, --processing-threads
    
    Percent here given as percent of total CPU on system which in our
    case should be the number of threads.
    Processes seem to be mutually exclusive (hard to tell though) and
    IO bound
    """
    num_threads=5
    r = min([4, num_threads])
    w = min([4, num_threads])
    p = num_threads
    return "-r {} -w {} -p {}".format(r, w, p)


def get_use_base(wildcards):	
	file_SScfg=FOLDER+'/SampleSheet_lane{}.cfg'.format(wildcards.lane)
	open(file_SScfg, 'a').close()	
	with open(file_SScfg,'r') as f:		
		for line in f:				
			if line.startswith('use_bases'):				
				use_bases=line.split('=')[1].strip()
				#print(use_bases)		
				return(use_bases)

def get_barcode_mode(wildcards):
	file_SScfg=FOLDER+'/SampleSheet_lane{}.cfg'.format(wildcards.lane)
	open(file_SScfg, 'a').close()
	with open(file_SScfg,'r') as f:		
		for line in f:						
			if line.startswith('barcode_mode'):				
				barcode_mode=line.split('=')[1].strip()
				print(barcode_mode)				
				return(barcode_mode)


def get_ref_genome(wildcards):
	file_SScfg=FOLDER+'/SampleSheet_lane{}.cfg'.format(wildcards.lane)
	open(file_SScfg, 'a').close()
	with open(file_SScfg,'r') as f:		
		for line in f:						
			if line.split('=')[0].strip()=='genome':				
				refgen=line.split('=')[1].strip()
				refgen_path='/g/solexa/bin/genomesNew/{}/{}.fa'.format(refgen,refgen)							
				return(refgen_path)



wildcard_constraints:
    lane="\d+"

rule all:
	input:
		all_samples


rule SS_convert:
	input: 
		folder=FOLDER,
		conf='/g/solexa/home/zawada/Desktop/SNAKE_MAKE/sampleSheet_convert/SampleSheetConverter.cfg'
	output:
		SS = expand(FOLDER+'/SampleSheet_lane{{lane}}.{filetype}', filetype = ["csv", "cfg"]),
	conda:
		"/g/solexa/home/zawada/Desktop/SNAKE_MAKE/environment.yaml"
	message:'Running rule SS_convert for lane {wildcards.lane}'	
	
	shell:
		'python2.7 sampleSheet_convert/SampleSheetConverter.py --run {input.folder} --lane {wildcards.lane} --config {input.conf} >> SS_output.txt 2>&1'



		
checkpoint Bcl2fastq:
	input:
		folder=FOLDER,
		SS_csv = rules.SS_convert.output.SS[0],
		bcl2fastq_path='/usr/local/bin/bcl2fastq' #use different path in case of runung on server		
	output:				
		samps_out=directory(FOLDER+'/Unaligned_lane{lane}')		
	log:
		FOLDER+'/Unaligned_lane{lane}/nohup.out'		
	params:
		use_bases=get_use_base,
		barcode_mode=get_barcode_mode,
		threads_n = bcl2fastq_threads_setting
		#use_bases='Y50,IIIIIIn'
	#resources:
		#mem_mb=100
	conda: #doesn't work with "run"
		"/g/solexa/home/zawada/Desktop/SNAKE_MAKE/environment.yaml"
	message:'Use Bases: {params.use_bases}, lane: {wildcards.lane}, barcode_mode: {params.barcode_mode}' 	
	run:
		if config['machine']=='M' or config['machine']=='H':			
			shell('{input.bcl2fastq_path} --tiles s_{wildcards.lane} --use-bases-mask {params.use_bases} {params.threads_n} --barcode-mismatches 0 --mask-short-adapter-reads 0 --minimum-trimmed-read-length 0 --runfolder-dir {input.folder} --output-dir {output.samps_out} --sample-sheet {input.SS_csv} > {log} 2>&1')
			if wildcards.lane in config['custom_lanes']:
				if len(reads)==2:					
					print(wildcards.lane+' is a custom lane')						
					shell('je demultiplex F1={input.folder}/Unaligned_lane{wildcards.lane}/lane{wildcards.lane}_S1_L00{wildcards.lane}_R1_001.fastq.gz F2={input.folder}/Unaligned_lane{wildcards.lane}/lane{wildcards.lane}_S1_L00{wildcards.lane}_R2_001.fastq.gz BF={input.folder}/log_pipeline/lane{wildcards.lane}_barcodes.txt BPOS=BOTH BM=BOTH BRED=false O={input.folder}/Unaligned_lane{wildcards.lane}')
				else:					
					print(wildcards.lane+' is a custom lane')
					shell('je demultiplex F1={input.folder}/Unaligned_lane{wildcards.lane}/lane{wildcards.lane}_S1_L00{wildcards.lane}_R1_001.fastq.gz BF={input.folder}/log_pipeline/lane{wildcards.lane}_barcodes.txt O={input.folder}/Unaligned_lane{wildcards.lane}')					

		elif config['machine']=='N':
			shell('{input.bcl2fastq_path} --use-bases-mask {params.use_bases} {params.threads_n} --no-lane-splitting --barcode-mismatches 0 --mask-short-adapter-reads 0 --minimum-trimmed-read-length 0 --runfolder-dir {input.folder} --output-dir {output.samps_out} --sample-sheet {input.SS_csv}  > {log} 2>&1')
			if wildcards.lane in config['custom_lanes']:
				if len(reads)==2:				
					print(wildcards.lane+' is a custom lane')					
					shell('je demultiplex F1={input.folder}/Unaligned_lane{wildcards.lane}/lane{wildcards.lane}_S1_R1_001.fastq.gz F2={input.folder}/Unaligned_lane{wildcards.lane}/lane{wildcards.lane}_S1_R2_001.fastq.gz BF={input.folder}/log_pipeline/lane{wildcard.lane}_barcodes.txt BPOS=BOTH BM=BOTH BRED=false O={input.folder}/Unaligned_lane{wildcards.lane}')					
				else:			
					print(wildcards.lane+' is a custom lane')
					shell('je demultiplex F1={input.folder}/Unaligned_lane{wildcards.lane}/lane{wildcards.lane}_S1_R1_001.fastq.gz BF={input.folder}/log_pipeline/lane{wildcards.lane}_barcodes.txt O={input.folder}/Unaligned_lane{wildcards.lane}')		
			### BELOW LINE is renaming files from nextseq
			shell('for f in {input.folder}/Unaligned_lane{wildcards.lane}/lane*; do mv "$f" ${{f//_R/_L001_R}}; done') # risky way of renaming, '_R' can be somewhere in the path
			shell('ls {input.folder}/Unaligned_lane{wildcards.lane}')



def agregate_input(wildcards):
	print(checkpoints.Bcl2fastq.get(lane=wildcards.lane))
	samples = samps_in_lanes_D_all[wildcards.lane]
	return(samples)

checkpoint move:
	input:
		rules.Bcl2fastq.output,
		folder=FOLDER,
		processed =  lambda wc: expand(FOLDER+'/Unaligned_lane{{lane}}/{samples}', samples = agregate_input(wc))	
	output:
		#to_align=expand(FOLDER+'/Aligned_lane{{lane}}/{{sample}}_{strand}_sequence.txt.gz', strand=reads) - this way an error associated with "the touch of death" stack overflow topic
		fake=FOLDER+'/fake/works{lane}.txt'
	shell:
		"""	
		for f in {input.folder}/Unaligned_lane{wildcards.lane}/lane*; do mv "$f" `echo "${{f//_R/_}}_sequence.txt.gz" | rev | cut -d"_" -f1,3,6- | rev` ; done	
		mv {input.folder}/Unaligned_lane{wildcards.lane}/lane* {input.folder}/Aligned_lane{wildcards.lane}
		touch {output.fake}
		"""

def get_unaligned(wildcards):
	checkpoint_output = checkpoints.move.get(**wildcards).output[0]
	return expand("{dir}/Aligned_lane{lane}/{sample}_{strand}_sequence.txt.gz",
		dir = FOLDER,
		lane = wildcards.lane,
		sample=wildcards.sample,
		strand = reads) 


#Below rule doesn't do anything, it's producing fake output. 
#I couldn't use there "run" instead of "shell", I think because of input with checkpoint,
#what is necessary for move rule ("touch of death" error) 
rule get_sample_wc: 
	input:
		to_align=get_unaligned,
		folder=FOLDER
	output:
		second_fake=FOLDER+'/fake/get_WC_fake/Aligned_lane{lane}/{sample}.bam'
		#bam=FOLDER+'/Aligned_lane{lane}/{sample}.bam'
	params:
		refgen=get_ref_genome
	message:'runing get_sample_wc for sample {sample}'
	shell:
		'touch {output.second_fake}'
		

rule align:
	input:
		rules.get_sample_wc.output.second_fake,
		folder=FOLDER		
	output:
		bam=FOLDER+'/Aligned_lane{lane}/{sample}.bam'
	params:
		rg=r"@RG\tID:{sample}\tSM:{sample}\tPL:illumina\tCN:EMBL",#"@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}\\tPL:illumina\\tCN:EMBL" in commandline works as well
		refgen=get_ref_genome,
		refgen_STAR=config['ref_gen_STAR'], # can be from config or function like for bwa refgen (above)
		prefix=FOLDER+'/Aligned_lane{lane}/{sample}' # name of output file will have to be changed to ...{sample}.bam so it matches output.bam
	message:'reference genome : {params.refgen}, lane: {wildcards.lane}'
	run:
		if wildcards.lane in config['STAR_lanes']:
			shell('STAR --runThreadN 1 --genomeDir {params.refgen_STAR} --outSAMtype BAM SortedByCoordinate --outFileNamePrefix {params.prefix} --readFilesCommand zcat --readFilesIn {params.prefix}_1_sequence.txt.gz')
			shell('for f in {params.prefix}*sortedByCoord.out.bam; do mv "$f" `echo "${{f//Aligned.sortedByCoord.out/}}"` ; done')
		else:
			if len(reads)==1:
				shell('bwa mem -R "{params.rg}" -t 1 {params.refgen} {params.prefix}_1_sequence.txt.gz 2>> {input.folder}/log_pipeline/gc_pipeline.lane{wildcards.lane}.bwa.err | samtools view -bT {params.refgen} - > {output.bam} 2>> {input.folder}/log_pipeline/gc_pipeline.lane{wildcards.lane}.samtools.err')
				
			else:		
				shell('bwa mem -R "{params.rg}" -t 1 {params.refgen} {params.prefix}_1_sequence.txt.gz {params.prefix}_2_sequence.txt.gz 2>> {input.folder}/log_pipeline/gc_pipeline.lane{wildcards.lane}.bwa.err | samtools view -bT {params.refgen} - > {output.bam} 2>> {input.folder}/log_pipeline/gc_pipeline.lane{wildcards.lane}.samtools.err')
		
		#below command will remove "_1" from fastq file names in case it's a single read. use wilcards in double {{}}, {{wildcards.lane}}/{{wildcards}}
		#if config['machine']=="N":
			#shell('for f in *txt.gz; do mv "$f" `echo "$f" | rev | cut -d"_" -f1,3- | rev` ; done') # remove "_1" from fastq file names in case it's a single read. use wilcards in double {{}}, {{wildcards.lane}}/{{wildcards}}

