# General
<p>This pipeline was written in response to the pandemic outbreak of the novel coronavirus known as COVID-19. The inability to test citizens at the rate required to stop or even slow the spread of this virus is very alarming to both scientists and citizens. The idea behind this pipeline is to improve the testing speed via the use of liquid handling robots combined with amplicon sequencing to massively parallelize sample processing...</p>

#### Basic COVID-19 testing facts:
<p>Currently, COVID-19 testing is conducted through a variety of molecular methods, including RT-qPCR, RT-PCR followed by standard PCR and gel electrophoresis, immunofluorescence microscopy, immunoblotting for viral protein, and CT scans of infected lung tissue. COVID-19 is an +ssRNA virus, therefore analysis of its genetic material requires reverse-transcription of the RNA --> complementary DNA (cDNA). From a molecular diagnostic prospective, the RT-qPCR is the most definitively accurate and quick method to perform. RT-qPCR is simply reverse transcription of the RNA into cDNA followed by standard PCR amplification of the cDNA products at a specific known viral genome location.
</p>

#### Ideas behind amplicon sequencing in COVID-19:
<p>Amplicon sequencing is simply target enriched sequencing. Instead of amplifying and quantifying one viral target in one patient per well (RT-qPCR), we would reverse transcribe and amplify multiple viral genome targets in each patient sample. Then, the samples will be combined into the same sample tube by patient and amplicon library assembly is performed using specific barcodes for each patient sample. To make the process safe and fast, we would need automated liquid handling robots to perform the heavy lifting/pipetting. After sequencing, we can demultiplex those sample reads to each individual, align them to the viral genome, and perform standard +/- virus analysis in addition to SNP calling and a myriad of other analyses.</p>  

#### Why amplicon sequencing instead of RT-qPCR:
I realized that using even the least capable next-generation sequencer  available (like an Illumina MiSeq, 15Gb w/25x10<sup>6</sup> reads) could not only drastically increase testing capacity for simple +/- diagnosis, but also increase its accuracy by using multiple viral amplicons. In addition, we may be to learn about the COVID-19 virus as it spreads and evolves via SNP analysis of our amplicons, instead of just getting +/- answers. To do this, we need to perform parallelized sequencing of multiple viral cDNAs in many patients simultaneously.</p>
<p>This pipeline will:</p>
<ul>
<li> Automatically demultiplex patient amplicon sequencing data </li>
<li> Align the sequencing data to the COVID-19 genome </li>
<li> Generate a custom report for each sample/individual describing the analysis results</li>
</ul>


#### Assumptions about sample throughput to feed into the analysis pipeline:
<p> I envision an automated process using liquid handling robots to extract viral RNA from thousands of patient samples, perform RT-PCR on them, amplify the viral amplicon targets (now cDNA) with PCR, build libraries and sequence.</p>

### Goals of this pipeline:
<ol>
<li>Reproducible</li>
<li>Scalable </li>
<li>Simple</li>
<li>Fast</li>
</ol>

### Reproducible & Scalable Analysis:
<p>The relatively simple operations of this pipeline are written in the <a href="https://www.nextflow.io">nextflow</a> pipeline development software to be run in a pre-built singularity container which is also freely available <a>here</a>. To further facilitate analysis in areas without strong computational infrastructure, I've developed a publicly available AWS image (AMI ID: ami-0681e8be831a1a855, AMI Name: COVID19 Targeted Sequencing Analysis Pipeline) which already has the necessary software to run the analysis built into the included singularity container.</p>


### Pipeline Minimal Compute Requirements:
<ul>
<li>CPU architecture: >= 2 </li>
<li>RAM: >= 8 Gb </li>
</ul>

# Installation:

<p>This pipeline is written using multiple layers of abstraction to provide high reproduciblity and portability. A user can choose from the following options to configure and run this pipeline: 
<ol>
<li>Run in the cloud with an <a href="https://aws.amazon.com/free/?trk=ps_a131L0000085EJuQAM&trkCampaign=acq_paid_search_brand&sc_channel=ps&sc_campaign=acquisition_US&sc_publisher=google&sc_category=core-main&sc_country=US&sc_geo=NAMER&sc_outcome=acq&sc_detail=amazon%20web%20services&sc_content=Brand_amazon_web_services_e&sc_segment=423740514695&sc_medium=ACQ-P|PS-GO|Brand|Desktop|SU|AWS|Core|US|EN|Text&s_kwcid=AL!4422!3!423740514695!e!!g!!amazon%20web%20services&ef_id=Cj0KCQjwx7zzBRCcARIsABPRscODB5HYuzBwvlVnnA5ob9O5LMgOlsdQer9H-vadHQlijFuRmHFPYXUaAtysEALw_wcB:G:s&all-free-tier.sort-by=item.additionalFields.SortRank&all-free-tier.sort-order=asc">Amazon Web Services</a> image</li>
<li>Run using <a href="https://sylabs.io/docs/">Singularity</a> container environments</li>
<li>Run after building the included <a href="https://docs.conda.io/en/latest/">conda</a> environment ./condaEnv/create_env.sh</li>
</ol>

### 1) Pipeline Minimal Software Requirements for running on AWS cloud resources:
<p>                            *****<strong>Nothing</strong>*****                                   </p>
<p>Simply use the publicly available <a href="https://aws.amazon.com/free/?trk=ps_a131L0000085EJuQAM&trkCampaign=acq_paid_search_brand&sc_channel=ps&sc_campaign=acquisition_US&sc_publisher=google&sc_category=core-main&sc_country=US&sc_geo=NAMER&sc_outcome=acq&sc_detail=amazon%20web%20services&sc_content=Brand_amazon_web_services_e&sc_segment=423740514695&sc_medium=ACQ-P|PS-GO|Brand|Desktop|SU|AWS|Core|US|EN|Text&s_kwcid=AL!4422!3!423740514695!e!!g!!amazon%20web%20services&ef_id=Cj0KCQjwx7zzBRCcARIsABPRscODB5HYuzBwvlVnnA5ob9O5LMgOlsdQer9H-vadHQlijFuRmHFPYXUaAtysEALw_wcB:G:s&all-free-tier.sort-by=item.additionalFields.SortRank&all-free-tier.sort-order=asc">Amazon Web Services</a> image AMI ID: ami-0681e8be831a1a855 with the AMI Name: COVID19 Targeted Sequencing Analysis Pipeline. It contains singularity and the singularity image file COVID19_Analysis.sif. Get your fastq data into the instance with the sabre compatible key (described below) and run using the: 'Via the prebuilt Singularity Container' command below.</p>


### 2) Pipeline Minimal Software Requirements for running locally with a singularity container:
<ul>
<li><a href="https://sylabs.io/docs/">Singularity</a> (version:3.5.3) and the prebuilt COVID19_Analysis.sif container</li>
</ul>
<p>To get the public prebuilt singularity container (COVID19_Analysis.sif) execute the below wget command:</p>
<p><code>wget https://covid19-amplicon-analysis-singularity-image.s3.us-east-2.amazonaws.com/COVID19_Analysis.sif</code></p>


### 3) Pipeline Minimal Software Requirements to build and run from a conda environment:
<p>Install the following functional programs in your $PATH:</p>
<ul>
<li><a href="https://www.nextflow.io">Nextflow</a> version: 20.01.0.5264</li>
<li><a href="http://bowtie-bio.sourceforge.net/bowtie2/index.shtml">Bowtie2</a> version: 2.3.5.1</li>
<li><a href="http://www.htslib.org">Samtools</a> version: 1.9</li>
<li><a href="https://github.com/najoshi/sabre">Sabre Demultiplexer</a> version: 1.0</li>
<li><a href="https://www.python.org/downloads/">Python</a> version: 3.7.6 with packages: pandas(v1.0.2), matplotlib (v3.1.3), numpy (v1.18.1), argparse (included with python dist) </li>
</ul>

<p><strong>If you plan to install everything instead of using the singularity instance I recommend using <a href="https://docs.conda.io/en/latest/">conda</a> environments as below</strong></p>

<p>First retrieve this git repo:</p>
<p><code>git clone https://github.com/Sandman2127/COVID19-Amplicon-Sequencing-Analysis-Pipeline.git</code></p> 
<p><code>cd COVID19-Amplicon-Sequencing-Analysis-Pipeline</code></p>
<p>Once conda is functional run the following command:</p>
<p><code>sh ./condaEnv/create_env.sh</code></p>
<p></p>
<p>This will create the conda environment COVID19-Amplicon-Seq. To activate it type:</p>
<p><code>source activate COVID19-Amplicon-Seq </code></p> 
<p></p>
<p>To test functionality type:</p>  
<p><code>bowtie2 --help</code></p> 

# Running Amplicon Sequencing Analysis:

#### Sabre barcode file:
<p>Should be a tab delimited text file with the barcode in the first column and the sampleName.fq (.fq is required after sample name):</p>
<p>CTCTCCAG <\t> RED_CLOUD_001.fq</p>
<p>TAATTG <\t> RED_CLOUD_002.fq</p>
<p></p>

#### FastQ input file:
<p>Any standard fastq file coming from an Illumina MiSeq, HiSeq, or NovaSeq should be able to run successfully in this setup</p>

#### Amplicon Sites Bedfile:
<p>This bed file provides the analysis program analyzeSam.py with the expected sites in the COVID19 genome where your amplicons are. It is in very basic bed format as below:</p>
<p>lcl|NC_045512.2_gene_1 <\t> 3921 <\t> 4071</p>
<p>lcl|NC_045512.2_gene_10 <\t> 421 <\t> 571</p>
<p>lcl|NC_045512.2_gene_3 <\t> 281 <\t> 431</p>
<p></p>
<p>If you intend to use the test dataset, you can find a test AmpliconSites.bed in the /lib/reportScripts/exampleBed/AmpliconSites.bed or you can simply copy it from <a href="https://github.com/Sandman2127/COVID19-Amplicon-Sequencing-Analysis-Pipeline/blob/master/lib/reportScripts/exampleData/test.bed">here</a>.</p>

## Getting and running the test data
<p>First retrieve this git repo:</p>
<p><code>git clone https://github.com/Sandman2127/COVID19-Amplicon-Sequencing-Analysis-Pipeline.git</code></p>
<p><code>cd COVID19-Amplicon-Sequencing-Analysis-Pipeline</code></p>
<p>to get the fastq test data:</p>
<p><code>wget https://covid-19-amplicon-test-data.s3.us-east-2.amazonaws.com/COVID19_simulated_SE_reads.fastq</code></p>
<p>to get the barcode data:</p>
<p><code>wget https://covid-19-amplicon-test-data.s3.us-east-2.amazonaws.com/COVID19_sample_barcode_file.txt</code></p>
<p>Jump to Run Nextflow Pipeline Command and run depending on your configuration</p>

## Run Nextflow Pipeline Command:

#### Via the prebuilt singularity Container
<p><code>singularity run /path/to/COVID19_Analysis.sif nextflow run /path/to/AnalyzeMultiplexedSamples.nf --genome ./COVID-19/genome --barcodes ./COVID19_sample_barcode_file.txt --bed /path/to/AmpliconSites.bed --inputF ./COVID19_simulated_SE_reads.fastq</code></p>

#### Without singularity, assuming all required programs are in the $PATH
<p><code>nextflow run /path/to/AnalyzeMultiplexedSamples.nf --genome ./COVID-19/genome --barcodes ./COVID19_sample_barcode_file.txt --bed /path/to/AmpliconSites.bed --inputF ./COVID19_simulated_SE_reads.fastq</code></p>

## Patient Data Outputs:


### Example Report Output:
https://sandman2127.github.io/COVID19-Amplicon-Sequencing-Analysis-Pipeline/

#### Main Diagnosis Result:
<img src="https://github.com/Sandman2127/COVID19-Amplicon-Sequencing-Analysis-Pipeline/blob/master/lib/exampleImages/MAINSCREEN.png" alt="Main Diagnostic Screen">

#### Read alignment quality:
<img src="https://github.com/Sandman2127/COVID19-Amplicon-Sequencing-Analysis-Pipeline/blob/master/lib/exampleImages/MAQ.png" alt="COVID19 amplicon alignment alignment quality image">

#### COVID19 amplicon alignment genome position:
<img src="https://github.com/Sandman2127/COVID19-Amplicon-Sequencing-Analysis-Pipeline/blob/master/lib/exampleImages/POS.png" alt="COVID19 amplicon alignment genome position">


### Conclusions
<p>Since I started with 25x10<sup>6</sup> simulated reads split across 5 amplicons in 96 separate samples, we expect the following: (25x10<sup>6</sup>/5)/96) == 52,631 reads per sample. Clearly the data above indicate perfect alignment of all reads in the sample. This will not be the case in real world samples, but I am confident we can determine a false positive rate for each real world COVID-19 primer set. If a sample has reads above this false positive rate at multiple loci we would consider a sample positive for COVID-19

## Performance:
<p>On a prebuilt dataset with 25x10<sup>6</sup> reads spread across 5 COVID19 amplicons from 96 samples with a 4 CPU 16 Gb (AWS m5a.xlarge) the analysis completes in 26 minutes @ a cost of $0.17/hr. I can easily see this scaling into 10s of thousands of samples processed per hour for less than $3 per hour.</p>

## How was the test data built:
##### Simulating COVID-19 amplicon sequencing

<p>I used <a href="https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm"> ART</a> to simulate 4.16x10<sup>6</sup> reads of 150 bp in length using Miseq V3 error profiles over 5 randomly chosen COVID amplicons using the below commands.</p>

##### Commands for ART data simulation

<p><code>art_illumina -ss MSv3 --samout -amp -na -i ./amplicons.fa -l 150 --fcov 5000000 -o single_ended_COVID_amplicon --rndSeed 127</code><p>

<p>I then barcoded each read using a custom script and a key file I had from previous GBS analysis. Those barcoded reads are found in our test data and can be effectively demultiplexed by sabre.</p>

