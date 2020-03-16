# General
<p>I wrote this pipeline in early March 2020 as a response to the pandemic outbreak of the novel coronavirus known as COVID-19. The inabilty to test citizens at the rate required to stop or even slow the spread of this virus was alarming and my idea to improve testing speed is simple. Use what we already have...</p>

#### Basic COVID-19 testing facts:
<p>Currently testing can be done via a variety of methods including immunofluorescence microscopy, RT-PCR followed by standard PCR and gel electrophoresis, Western blots for viral protein and CT scans of infected lung tissue. From a molecular diagnostic prespective the RT-PCR method is the most accurate and quick to perform. COVID-19 is an RNA virus, and as such its detection requires reverse transcription of the RNA into DNA before PCR or subsequent sequencing.</p> 

#### My plan:
<p>I designed this pipeline to take the RT-PCR method one step further and perform sequencing of the viral cDNA. I realized that using even the least capable next generation sequencer available (like an Illumina MiSeq, 15Gb w/25x10^6 reads) could drastically increase testing capacity. To do this we need to perform parallelized sequencing of multiple viral cDNAs in many patients simultaneously.</p>
<p>This pipeline will:</p>
<ul>
<li> Automatically demultiplex patient amplicon sequencing data </li>
<li> Align the sequencing data to the COVID-19 genome </li>
<li> Generate a custom report for each sample/individual describing the analysis results</li>
</ul>


#### Assumptions about sample throughput to feed into the analysis pipeline:
<p> I envision an automated process using liquid handling robots to extract viral RNA from thousands of patient samples, perform RT-PCR on them, amplify the viral amplicon targets (now cDNA) with PCR, build libraries and sequence.</p>

#### Goals of this pipeline:
<ol>
<li>Reproducible</li>
<li>Scalable </li>
<li>Simple</li>
<li>Fast</li>
</ol>

#### Reproducible & Scalable Analysis:
<p>The relatively simple operations of this pipeline are written in the <a href="https://www.nextflow.io">nextflow</a> pipeline development software to be run in a pre-built singularity container which is also freely available <a>here</a>. To further facilitate analysis in areas without strong computational infrastructure I've developed a publicly available AWS image (AMI ID: ami-0681e8be831a1a855, AMI Name: COVID19 Targeted Sequencing Analysis Pipeline) which already has the necessary software to run the analysis built into the included singularity container.</p>


#### Pipeline Minimal Compute Requirements:
<ul>
<li>CPU architecture: >= 2 </li>
<li>RAM: >= 8 Gb </li>
</ul>

#### Pipeline Minimal Software Requirements for running on AWS cloud resources:
<p>                            *****<strong>Nothing</strong>*****                                   </p>
<p>Simply use the publicly available image AMI ID: ami-0681e8be831a1a855 with the AMI Name: COVID19 Targeted Sequencing Analysis Pipeline. It contains singularity and the singularity image file COVID19_Analysis.sif. Get your fastq data into the instance with the sabre compatible key (described below) and run using the: 'Via the prebuilt Singularity Container' command below.</p>


#### Pipeline Minimal Software Requirements for running locally w/Singularity:
<ul>
<li><a href="https://sylabs.io/docs/">Singularity</a> (version:3.5.3) and the prebuilt COVID19_Analysis.sif container</li>
</ul>
<p>To get the public prebuilt singularity container (COVID19_Analysis.sif) <strong>execute the below command on a machine with AWScli</strong></p>
<p>aws s3 cp s3://covid19-amplicon-analysis-singularity-image/COVID19_Analysis.sif ~ </p>


#### Pipeline Minimal Software Requirements without singularity image file:
<p>Install the following functional programs in your path</p>
<ul>
<li><a href="https://www.nextflow.io">Nextflow</a> version:20.01.0.5264</li>
<li><a href="http://bowtie-bio.sourceforge.net/bowtie2/index.shtml">Bowtie2</a> version:>=2.3.5.1</li>
<li><a href="http://www.htslib.org">Samtools</a> version:1.9</li>
<li><a href="https://github.com/najoshi/sabre">Sabre Demultiplexer</a> version:1.0</li>
<li><a href="https://www.python.org/downloads/">Python</a>version: 3.7.6 with packages: pandas(v1.0.2), matplotlib (v3.1.3), numpy (v1.18.1), argparse (included with python dist) </li>
</ul>
<p>If you plan to install everything instead of using the singularity instance I recommend using <a href="https://docs.conda.io/en/latest/">conda</a> environments</p>


## Running Amplicon Sequencing Analysis:
#### Sabre barcode file:
<p>Should be a tab delimited text file as below:</p>
<p>CTCTCCAG RED_CLOUD_001.fq</p>
<p>TAATTG   RED_CLOUD_002.fq</p>
<p>...</p>
<p></p>

#### FastQ input file:
<p>Any standard fastq file coming from an Illumina MiSeq, HiSeq or NovaSeq should be able to run successfully in this setup</p>

## Run Nextflow Pipeline Command:

#### Via the prebuilt Singularity Container
<p>Singularity run /path/to/.sif nextflow run /path/to/AnalyzeMultiplexedSamples.nf --genome /path/to/COVID-19/genome --barcodes /path/to/sabre_barcode_file.txt --inputF /path/to/multiplexed.fastq </p>

#### Without singularity, assuming all required programs are in the $PATH
<p>nextflow run /path/to/AnalyzeMultiplexedSamples.nf --genome /path/to/COVID-19/genome --barcodes /path/to/sabre_barcode_file.txt --inputF /path/to/multiplexed.fastq</p>

## Patient Data Outputs:

## Performance:
<p>On a prebuilt dataset with 25 million reads spread across 6 COVID19 amplicons from 96 samples with a 4 CPU 16 Gb (AWS m5a.xlarge) the analysis completes in 26 minutes @ a cost of 0.17/hr. I can easily see this scaling into 10s of thousands of samples processed per hour for less than $3 per hour.</p>

##### Simulating COVID-19 amplicon sequencing

I used <a href="https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm"> ART</a> to simulate 4.166 Million reads of 150 bp in length using Miseq V3 error profiles over 6 COVID amplicons using the below commands.

##### Commands for ART data simulation

art_illumina -ss MSv3 --samout -amp -na -i ./amplicons-0-5.fa -l 150 --fcov 4166666 -o single_ended_COVID_amplicon --rndSeed 127
