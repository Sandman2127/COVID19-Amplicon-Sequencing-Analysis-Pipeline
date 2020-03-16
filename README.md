# General
<p>I wrote this pipeline in early March 2020 as a response to the pandemic outbreak of the novel coronavirus known as COVID-19. The inabilty to test citizens at the rate required to stop or even slow the spread of this virus was alarming and my idea to improve testing speed is simple. Use what we already have...</p>

#### Basic Testing facts:
<p>Currently testing can be done via a variety of methods including immunofluorescence microscopy, RT-PCR followed by standard PCR and gel electrophoresis, Western blots for viral protein and CT scans of infected lung tissue. From a molecular diagnostic prespective the RT-PCR method is the most accurate and quick to perform. COVID-19 is an RNA virus, and as such its detection requires reverse transcription of the RNA into DNA before sequencing.</p> 

#### My plan:
<p>I designed this pipeline to take the RT-PCR method one step further and perform sequencing of the viral cDNA. I realized that using even the least capable next generation sequencer available (like an Illumina MiSeq, 15Gb w/25x10^6 reads) could drastically increase testing capacity. To do this we need to perform parallelized sequencing of multiple viral cDNAs in many patients simultaneously.</p>
<p>This pipeline will:</p>
<ul>
<li> Automatically demultiplex patient amplicon sequencing data </li>
<li> Align the sequencing data to the COVID-19 genome </li>
<li> Generate a custom report for each sample/individual describing the analysis results</li>
</ul>


#### Assumptions about sample throughput to feed into the analysis pipeline:
<p> I envision an automated process using liquid handling robots to extract viral RNA from thousands of patient samples, perform RT-PCR on them, amplify the viral amplicon targets with PCR, build libraries and sequence.</p>


#### In order for this analysis pipeline to benefit everyone it needs to have the following charecteristics:


#### Goals of this pipeline:
<ol>
<li>Reproducible</li>
<li>Scalable </li>
<li>Simple</li>
<li>Fast</li>
</ol>

### Simulating COVID-19 amplicon sequencing

I used ART to simulate 4.166 Million reads of 150 bp in length using Miseq V3 error profiles over 6 COVID amplicons

##### Commands for ART data simulation

art_illumina -ss MSv3 --samout -amp -na -i ./amplicons-0-5.fa -l 150 --fcov 4166666 -o single_ended_COVID_amplicon --rndSeed 127


