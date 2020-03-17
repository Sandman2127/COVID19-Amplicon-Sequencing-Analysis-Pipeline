# General
<p>This pipeline was written in response to the pandemic outbreak of the novel coronavirus known as COVID-19. The inabilty to test citizens at the rate required to stop or even slow the spread of this virus is alarming and my idea to improve testing speed is simple. Use amplicon sequencing and combine it with liquid handling robots to massively parallelize sample processing...</p>

#### Basic COVID-19 testing facts:
<p>Currently COVID-19 testing can be done via a variety of molecular methods including: immunofluorescence microscopy, RT-PCR followed by standard PCR and gel electrophoresis, western blots for viral protein and CT scans of infected lung tissue. From a molecular diagnostic prespective the RT-PCR method is the most accurate and quick to perform. COVID-19 is an RNA virus, and as such its detection requires reverse transcription of the RNA into DNA before PCR or subsequent sequencing.</p>

#### Background ideas:
<p>This pipeline was designed to take the RT-PCR method a step further and analyze sequencing data of the viral cDNA. I realized that using even the least capable next generation sequencer available (like an Illumina MiSeq, 15Gb w/25x10^6 reads) could not only drastically increase testing capacity for simple +/- diagnosis, but we may be also be able to learn about the COVID-19 virus as it spreads and evolves by sequencing, instead of just performing electorphoresis. To do this we need to perform parallelized sequencing of multiple viral cDNAs in many patients simultaneously.</p>
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
<p>The relatively simple operations of this pipeline are written in the <a href="https://www.nextflow.io">nextflow</a> pipeline development software to be run in a pre-built singularity container which is also freely available <a>here</a>. To further facilitate analysis in areas without strong computational infrastructure I've developed a publicly available AWS image (AMI ID: ami-0681e8be831a1a855, AMI Name: COVID19 Targeted Sequencing Analysis Pipeline) which already has the necessary software to run the analysis built into the included singularity container.</p>


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
<li>Run after building the included <a href="https://docs.conda.io/en/latest/">conda</a> environment ./condaEnv/COVID19_Analysis_conda_environment.yml</li>
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
<p>Install the following functional programs in your path</p>
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
<p>...</p>
<p></p>

#### FastQ input file:
<p>Any standard fastq file coming from an Illumina MiSeq, HiSeq or NovaSeq should be able to run successfully in this setup</p>

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
<p><code>singularity run /path/to/COVID19_Analysis.sif nextflow run /path/to/AnalyzeMultiplexedSamples.nf --genome ./COVID-19/genome --barcodes ./COVID19_sample_barcode_file.txt --inputF ./COVID19_simulated_SE_reads.fastq</code></p>

#### Without singularity, assuming all required programs are in the $PATH
<p><code>nextflow run /path/to/AnalyzeMultiplexedSamples.nf --genome ./COVID-19/genome --barcodes ./COVID19_sample_barcode_file.txt --inputF ./COVID19_simulated_SE_reads.fastq</code></p>

## Patient Data Outputs:

#### <em>Read alignment by contig and bowtie2 mapping quality score (MAQ)</em>
<img src="https://github.com/Sandman2127/COVID19-Amplicon-Sequencing-Analysis-Pipeline/blob/master/lib/resultImages/Amplicon_readMAQ.png"/>
<p>The mapping quality of reads by contig can tell us how accurately our amplified viral sequences are matching the COVID-19 genome. High numbers of mapping reads on multiple amplicons indicates a positive diagnosis for COVID-19, low mapping at almost all will indicate a negative result (i.e. no COVID-19 detected). This simulated data maps very highly to the genome, <strong>in real life</strong> we will need to empirically determine the false positive alignment rate of each primer set. However, assuming we can generate high specificity primers for multiple viral genome sites, we can likely drastically improve the sensitivity of COVID-19 detection.</p>
<p></p>
<p></p>

#### <em>Histogram showing starting position of mapping reads by contig</em>
<img src="https://github.com/Sandman2127/COVID19-Amplicon-Sequencing-Analysis-Pipeline/blob/master/lib/resultImages/Amplicon_Pos.png"/>
<p>Theoretically in an amplicon sequencing experiment one should see a tight distribution. Ideally a single large peak of read start mapping positions. This plot immediately tells the reviewer what site in the viral genome the reads are mapping to, and how many. With background knowledge of the expected primer anealing sites, expected read depth and false positive rate a knowledgable reviewer would instantly know if the alignments indicate a positive or negative COVID-19 result.</p>

## Performance:
<p>On a prebuilt dataset with 25 million reads spread across 5 COVID19 amplicons from 96 samples with a 4 CPU 16 Gb (AWS m5a.xlarge) the analysis completes in 26 minutes @ a cost of $0.17/hr. I can easily see this scaling into 10s of thousands of samples processed per hour for less than $3 per hour.</p>

## How was the test data built:
##### Simulating COVID-19 amplicon sequencing

<p>I used <a href="https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm"> ART</a> to simulate 4.166 Million reads of 150 bp in length using Miseq V3 error profiles over 5 randomly chosen COVID amplicons using the below commands.</p>

##### Commands for ART data simulation

<p><code>art_illumina -ss MSv3 --samout -amp -na -i ./amplicons-0-5.fa -l 150 --fcov 4166666 -o single_ended_COVID_amplicon --rndSeed 127</code><p>

<p>I then barcoded each read using a custom script and a key file I had from previous GBS analysis. Those barcoded reads are found in our test data and can be effectively demultiplexed by sabre.</p>

