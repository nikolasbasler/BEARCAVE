# BEARCAVE

A collective data repository for archiving, trimming, mapping and
maintaining Illumina resequencing data

**Purpose and Philosophy**

Many researchers who work with Illumina resequencing data will write
their personal computational pipelines that can process their data in an
automated manner. However, when working as a group on shared or similar
datasets it may be necessary to agree on a fixed workflow to keep the
results of different projects consistent and comparable. In such cases,
the storage of both raw and processed data files in a centralised shared
data repository accessible by the group provides a more efficient
approach than maintaining duplicated, redundant copies of data files in
each individual’s home directories or local systems.

The BEARCAVE provides an environment for establishing and maintaining
such a shared data repository. The BEARCAVE distribution includes all
necessary scripts and will automatically download and install the
necessary software via the Miniconda package and environment manager
([<span
class="underline">https://conda.io/docs/index.html</span>](https://conda.io/docs/index.html)): Cutadapt \[1\] for adapter trimming and short
read removal, FLASH \[2\] for merging paired-end reads, BWA \[3\] for
mapping to a reference genome and SAMtools \[4\] for sorting, quality
filtering and duplicate removal. It therefore can be easily used on multiple servers or systems, while ensuring
consistent results across systems and also also provides a simple user-friendly
method for inserting and archiving raw data files. The generated results and processed data files are stored within a logical folder structure to provide a framework
for collective research work. Moreover, all datasets
are automatically assigned with a unique identification code such that
the sequence of processing and scripts can be tracked from final bam
file all the way back to the original raw data file(s). A script is also
included to automatically generate a list of all data sets in the
BEARCAVE and their stage of processing. The workflow is not completely
automated and requires user-input at certain key steps. This helps to
avoid accidentally overwriting valuable data and allows for better
troubleshooting, in case it is needed. The scripts are also constructed
in a way to help avoid inadvertent use of incorrect copies or versions
of software programs or data files.

We originally developed the BEARCAVE as a shared resource for multiple
researchers working with large datasets from ancient and modern bears
(*Ursidae*). However, it is applicable to any taxonomic group and is
fully customisable allowing optimisation for different research
requirements and environments. The framework and workflow of the
BEARCAVE has been optimised based on years of experience of the
practicalities of working with such data, and most importantly what can
go wrong and how it can be avoided. Scripts are provided for the
processing of single-end and paired-end sequencing of both
single-stranded and double-stranded libraries, allowing changes in
minimum read length threshold and mapping specificity. We find these
alternatives accommodate all our needs, but scripts can be easily
modified or new scripts added should additional functions be required.

**Important notes:**

-   **How to cite**: If you publish results that have been generated
    using the BEARCAVE, please refer to this article: &lt;not published
    yet&gt;

    -   **Note**, however, that the BEARCAVE does not contain any
        original software. Please also refer to the appropriate
        articles for Cutadapt \[1\], FLASH \[2\], BWA \[3\] and SAMtools
        \[4\] where applicable.

-   All scripts are located in and must be executed from the
    BEARCAVE/scripts/ folder. They are executable Bash scripts. To run
    them, simply navigate to BEARCAVE/scripts and type:

    `./script.sh`

-   View the scripts you want to use. At the top of each script there is
    a section with a short explanation on how to use it.

-   All scripts that include merging of paired-end reads or mapping are
    set to use **10 CPU cores**. You may want to adjust this setting for
    your system by changing the “threads” variable in the top section of
    the appropriate scripts.

-   The BEARCAVE was tested for Scientific Linux 6.10, Fedora 29, Debian
    9 stretch, Ubuntu 16.04 and Ubuntu 18.04.

-   Permissions are automatically set so that only members of the
    installer’s user group can use and access the BEARCAVE.

-   The examples given in this manual refer to a test dataset and
    reference genome that come with the BEARCAVE package.

**Known issues**

-   Having multiple identical and simultaneously running jobs (i.e.
    identical scripts and identical input data) will overwrite each
    other’s output files and therefore break them. This is also true
    if the jobs are issued by different users.

**Installation**

-   Download the BEARCAVE package from GitHub ([<span
    class="underline">https://github.com/nikolasbasler/BEARCAVE</span>](https://github.com/nikolasbasler/BEARCAVE)),
    copy it into a folder of your choice and extract it:

    `tar -zxvf BEARCAVE_v*tar.gz` This will create a new folder called BEARCAVE.

-   Navigate to BEARCAVE/ and execute the installation script:

    `./install.sh`

-   You will have to accept the licence agreement of Miniconda and and
    otherwise confirm the default settings by hitting ENTER.

-   After the installation is done, you will find the software, scripts
    and the following sub-folder structure within BEARCAVE/:

**BEARCAVE/**
- **rawdata/**
- **old\_metadata/**
- **refgenomes/**
- **scripts/**
- **software/**
- - **miniconda3/**
- - - **bin/**
- - - **&lt;several more subfolders&gt;**
- **test\_data/**
- **trimdata/**
- - **trimlogs/**

**Indexing Reference Genomes**

The BEARCAVE can not only be used for bears. You can download any
reference genome you wish to use for mapping or use the polar bear
mitochindrial genome provided in /test\_data/:

-   Create a folder for the reference genome you want to use in
    /refgenomes/, and copy your reference genome into this folder,e.g.

    > **/refgenomes/PolarBear_mt/**

-   Navigate to /scripts/ and execute the index\_ref script with the
    path and file name of your new reference genome as argument, e.g.:

    `./index_ref.sh ../refgenomes/PolarBear_mt/polarbear_mt.fasta.gz`

-   The reference genome will be indexed with BWA and SAMtools and a new
    “mapped” folder will automatically be created:

    > **/mappedPolarBear_mt/**

**Inserting Sequencing Data**

-   Create a new folder in /rawdata/ named like the sample your data
    originates from (**this folder name should be used whenever a script
    asks for "sample"**), e.g.
    >**/rawdata/BearA/**

-   Copy your zipped FASTQ files into the new folder. They might look
    like this (this particular example file can be found in /test\_data/):

    > **/rawdata/BearA/BearA-01_S2_L001_R1_001.fastq.gz**

-   Note that file names must not contain any + (plus) symbols.

-   Create a comma separated table containing metadata about the
    dataset, such as the sample code, library and sequencing run. This
    table must have the following entries as header and should be filled
    out accordingly (view the /scripts/add\_to\_meta.sh for more
    information on the format). Missing data should be inserted as NA.
    For an example metadata file, have a look at /test\_data/BearA.meta.
    <table>
    <tr class="header">
    <td><p><center>SEQ_RUN</p></td>
    <td>An identifier for the original file. <b> SEQ_RUN has to be the beginning of the file name</b> (e.g. BearA-01_S2 for the file in the example above).</td>
    </tr>
    <tbody>
    <tr class="odd">
    <td><p><center>SAMPLE</p></td>
    <td>The sample name. Identical to the folder name in /rawdata/.</td>
    </tr>
    <tr class="odd">
    <td><p><center>TAXON</p></td>
    <td>The samples taxon.</td>
    </tr>
    <tr class="odd">
    <td><center>LOCALITY</p></td>
    <td>Where the sample has been found.</td>
    </tr>
    <tr class="odd">
    <td><p><center>COUNTRY</p></td>
    <td>The country the sample is from.</td>
    </tr>
    <tr class="odd">
    <td><p><center>AGE</center></p></td>
    <td>Modern, historical or ancient.</td>
    </tr>
    <tr class="odd">
    <td><p><center>DATABASE_NO</center></p></td>
    <td>If the data has been published, state the accession number here, if not use "unpublished" instead.</td>
    </tr>
    <tr class="odd">
    <td><p><center>LIBRARY_NO</center></p></td>
    <td>The library number for internal reference.</td>
    </tr>
    <tr class="odd">
    <td><p><center>EXTRACT_METH</center></p></td>
    <td>Method used for DNA extraction (e.g. “Dab” for Dabney [5]).</td>
    </tr>
    <tr class="odd">
    <td><p><center>LIBRARY_METH</p></td>
    <td>Single-stranded (SS) or double-stranded (DS).</td>
    </tr>
    <tr class="odd">
    <td><p><center>PLATFORM</p></td>
    <td>The device used for sequencing.</td>
    </tr><tr class="odd">
    <td><p><center>READ_LENGTH</p></td>
    <td>Maximum read length used during the sequencing run.</td>
    </tr>
    <tr class="odd">
    <td><p><center>MODE</center></p></td>
    <td>Single-end (SE) or paired-end (PE).  </td>
    </tr>
    <tr class="odd">
    <td><p><center>SEQ_PRIMER</center></p></td>
    <td>The primer used for sequencing (e.g. “standard” or CL72 [5]).</td>
    </tr>
    <tr class="odd">
    <td><p><center>RIGHTS</center></p></td>
    <td>Who should be contacted before publishing this data (separate names with an underscore, avoid commas and blank spaces).</td>
    </tr>
    </table>


-   To add the content of your table to the BEARCAVE’s metadata file and
    prepare your raw data for further use, navigate to /scripts/ and
    execute the add\_to\_meta.sh using the newly created table as
    argument, e.g.

    `./add_to_meta.sh ../test_data/BearA.meta`

-   The entries from your table can now be found in the
    /rawdata/metadata.txt file together with a unique three-character
    identifier that will serve as a file name prefix. Consult this file
    for information about all previously entered datasets. If you open
    the file in a spreadsheet viewer (such as LibreOffice Calc), make
    sure to “format quoted fields as text”, otherwise some entries might
    accidentally be recognised as numbers and automatically change.

-   Note that the prefix has automatically been added to your rawdata
    file name, which now might  like this (since the prefix is
    generated randomly, the prefix in your file will very likely be
    different):

    `/rawdata/BearA/pr1+BearA-01_S2_L001_R1_001.fastq.gz`

-   The addition of the prefix will lead to a different order in which
    the files are listed. To view them ordered alphabetically, ignoring
    the prefix, use this command:

    `ll | sort -t '+' -k 2`

-   You can add several datasets at once. Simply include a line for each
    dataset in the table you have created.

**Trimming (and paired-end merging)**

-   All adapter trimming is performed by Cutadapt (version 1.12)
    \[1\] with an adapter overlap value (‑O) of 1. If not specifically
    stated otherwise, a minimum read length (-m) of 30 is used.

-   For merging of paired-end data, FLASH (version 1.2.11)
    \[2\] is used with a maximum overlap value for read pairs (-M) of 75.

-   Note that all paths to files and programs are hard-coded into the
    scripts. This is to reduce the potential for alternative software
    versions or data files located outside the BEARCAVE being used for
    analysis.

-   Navigate to /scripts/ and choose the appropriate trimming script:

<table>
<tr class="header">
<td><p>trim_merge_DS_PE_CL72.sh</p>
<p>Software: Cutadapt [1], FLASH [2]</p></th>
<td>For trimming and merging double-stranded libraries from a paired-end sequencing run that used the CL72 sequencing primer [6]. <strong>Arguments:</strong> SAMPLE, PREFIX, SEQ_RUN. <strong>Output:</strong> 1 log file for trimming, 1 log file for merging, 3 FASTQ files (mappable_R1.fastq and mappable_R2.fastq contain unmerged reads, mappable.fastq contains merged reads). Process output with a PE mapping script. If your data is from an ancient sample, you may want to delete the unmerged files and consider only merged data. Then the output has to be processed with an SE mapping script.</th>
</tr>
<tbody>
<tr class="odd">
<td><p>trim_merge_DS_PE_standard.sh</p>
<p>Software: Cutadapt [1], FLASH [2]</p></td>
<td>For trimming and merging double-stranded libraries from a paired-end sequencing run that used the standard Illumina sequencing primer. <strong>Arguments/Output:</strong> See trim_merge_DS_PE_CL72.sh.</td>
</tr>
<tr class="even">
<td><p>trim_merge_SS_PE_CL72.sh</p>
<p>Software: Cutadapt [1], FLASH [2]</p></td>
<td>For trimming and merging single-stranded libraries from a paired-end sequencing run that used the CL72 sequencing primer [6]. <strong>Arguments:</strong> SAMPLE, PREFIX, SEQ_RUN. <strong>Output:</strong> 1 log file for trimming, 1 log file for merging, 1 FASTQ file (files containing unmerged reads are deleted as they are not used for mapping ancient data). Process output with an SE mapping script (because merged reads are not paired anymore).</td>
</tr>
<tr class="odd">
<td><p>trim_merge_SS_PE_CL72_var_min_length.sh</p>
<p>Software: Cutadapt [1], FLASH [2]</p></td>
<td>For trimming with a custom read length cutoff and merging single-stranded libraries from a paired-end sequencing run that used the CL72 sequencing primer [6]. <strong>Arguments:</strong> SAMPLE, PREFIX, SEQ_RUN, minimum read length. <strong>Output:</strong> See trim_merge_SS_PE_CL72.sh.</td>
</tr>
<tr class="even">
<td><p>trim_SE.sh</p>
<p>Software: Cutadapt [1]</p></td>
<td>For trimming any libraries from a single-end sequencing run. <strong>Arguments:</strong> SAMPLE, PREFIX, SEQ_RUN. <strong>Output:</strong> 1 log file for trimming, 1 FASTQ file. Process output with an SE mapping script.</td>
</tr>
<tr class="odd">
<td><p>trim_SE_var_min_length.sh</p>
<p>Software: Cutadapt [1]</p></td>
<td>For trimming any libraries from single-end sequencing run with a custom minimum read length cutoff. <strong>Arguments:</strong> SAMPLE, PREFIX, SEQ_RUN, minimum read length. <strong>Output:</strong> See trim_SE.sh.</td>
</tr>
</tbody>
</table>

-   Run the appropriate script, e.g. like this:

    `./trim_SE.sh BearA pr1 BearA-01_S2`

-   When the script is done, you will find a new "processing" folder in
    /trimdata/ containing the output files mentioned in the table above.

    > **/trimdata/BearA_processing/**

-   View the log file(s) and confirm the trimming went as expected.

-   Move the log file(s) into /trimdata/trimlogs/ and the zipped FASTQ
    file(s) one folder up into /trimdata/.

-   Delete the „processing" folder, which should now be empty.

**Combining Datasets (optional)**

-   Navigate to /scripts/ and call the combine\_files.sh script using
    the sample name, the prefixes of the files you want to combine and
    either ‘mappable’ (for merged read pair files), ‘mappable\_R1’ (for
    unmerged R1 files) or ‘mappable\_R2’ (for unmerged R2 files) as
    arguments, e.g.:

    `./combine_files.sh BearA mappable pr1 pr2`

-   If one of the files already contains combined datasets, use its
    whole set of prefixes as one argument. E.g., if you want to combine
    the files __pr1+BearA_mappable.fastq.gz__ and __pr2\_pr3+BearA_mappable.fastq.gz__
    call the script like this:

    `./combine_files.sh BearA mappable pr1 pr2_pr3`

-   Note, that if you are using files from a read length cutoff other
    than 30, you will have to extend the sample name argument according
    to the file, e.g. BearA_28bp

-   When the files are combined, the script will ask you if the original
    files should be deleted. This is recommended, as the original files
    are redundant to the combined ones and are therefore unnecessarily
    occupying hard disk space. Keep the cave clean.

-   The combined output file can be found in a "processing" folder in
    /trimdata/.

-   The order of prefixes in the file name reflects the order of the
    sequences in the file.

-   Move the combined file one folder up into /trimdata/ and delete the
    "processing" folder, which should now be empty.

**Mapping**

-   All mapping is performed by BWA (version 0.7.15) \[3\] "aln" and "samse" (for
    single-end data) or "sampe" (for paired-end data). If not
    specifically stated otherwise, a mismatch parameter (-n) of 0.04 is
    used.

-   The alignment is edited by SAMtools (version 1.3.1.) \[4\]: Low quality reads are removed by SAMtools "view"
    with a threshold (-q) of 30, the alignment is sorted with "sort" before duplicate reads are collapsed with "rmdup" and the alignment is indexed with "index".

-   Like in the trimming and merging scripts, all paths to files and
    programs are hard-coded into the scripts. This is to reduce the
    potential for alternative software versions or data files located
    outside the BEARCAVE being used for analysis.

-   Navigate to /scripts/ and chose the appropriate mapping script:

<table>
<thead>
<tr class="header">
<td><p>map_SE_0.01mismatch.sh</p>
<p>Software: BWA [3], SAMtools [4]</p></th>
<td>For mapping single-end data with a BWA mismatch value (-n) of 0.01. <strong>Arguments:</strong> PREFIX, reference genome<strong>*</strong>, 3‑character taxon identifier, SAMPLE.</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td><p>map_SE.sh</p>
<p>Software: BWA [3], SAMtools [4]</p></td>
<td>For mapping single-end data. <strong>Arguments:</strong> PREFIX, reference genome<strong>*</strong>, 3‑character taxon identifier, SAMPLE.</td>
</tr>
<tr class="even">
<td><p>map_modern_PE_0.01mismatch.sh</p>
<p>Software: BWA [3], SAMtools [4]</p></td>
<td>For mapping (modern) paired-end data with a BWA mismatch value (-n) of 0.01. <strong>Arguments:</strong> PREFIX, reference genome<strong>*</strong>, 3‑character taxon identifier, SAMPLE.</td>
</tr>
<tr class="odd">
<td><p>map_modern_PE.sh</p>
<p>Software: BWA [3], SAMtools [4]</p></td>
<td>For mapping modern paired-end data. <strong>Arguments:</strong> PREFIX, reference genome<strong>*</strong>, 3‑character taxon identifier, SAMPLE.</td>
</tr>
</tbody>
</table>

<strong>*</strong> Refers to the folder name in /refgenomes/.

<!-- -->

-   Run the appropriate mapping script e.g. like this:

    `./map_SE.sh pr1 PolarBear_mt urs BearA`

-   If you are mapping a combined data file, you should enter its whole
    set of prefixes as one argument, like this: pr1\_pr2\_pr3.

-   When the script is done, you will find a "processing" folder in the
    appropriate “mapped“ folder in BEARCAVE/, containing the output
    files mentioned in the table above, e.g.:

> **/mappedPolarBear_mt/pr1+BearA_panda_map_processing/**

-   View the log file and confirm the mapping has worked as expected and
    then move it to the "logs" folder of the respective "mapped" folder,
    e.g. like this:

    `mv *.log ../PolarBear_mt_logs/`

-   Also contained in this folder, there is a .bam file (your alignment)
    and a .bam.bai file (the index of the alignment). Their file names
    consist of the prefix of the dataset, the sample name, the taxon
    identifier, the reference used for mapping and the total Gb that
    mapped, e.g. like this:

> **pr1+BearA_urs_PolarBear_mt_.00666.bam**

-   Move these files (both the .bam and the .bam.bai) one folder up into
    its appropriate "mapped" folder and delete the "processing" folder,
    which should now be empty.

-   The file naming system allows for easy parsing of the accumulated
    set of bam files, e.g. this command lists all bam files for the
    taxon “urs”

    `ls -1 *urs*bam`

**State of the Cave (optional)**

-   For a detailed list of all raw data that entered the BEARCAVE, view
    the metadata file /rawdata/metadata.txt.

-   To get an overview of which datasets have been trimmed and mapped
    (and to which references), navigate to /scripts/ and execute the
    cavescan.sh without any arguments.

-   This will generate a comma separated table that lists all samples
    and prefixes associated with them and if they have been trimmed
    and/or mapped to the references or not (based on the existence of
    files in the respective "mapped" folders): **/state\_of\_the\_cave.txt**

-   When opening the table in a spreadsheet viewer (such as LibreOffice
    Calc), make sure to “format quoted fields as text”, otherwise some
    identifiers might accidentally be recognised as numbers and
    automatically change.

**Quick Guide for Repeated Use**

-   **Inserting data:**

    -   If not already existing, create a folder for your sample in
        /rawdata/.

    -   Copy the zipped FASTQ files of your sequencing data there. File
        names must not contain a + (plus) symbol.

    -   Create a comma separated table containing the metadata for your
        dataset. Consult the above paragraph and the
        /scripts/add_to_meta.sh for more information on the required
        format and the /test_data/BearA.meta for an example.

    -   Navigate to /scripts/ and execute the add_to_meta.sh with your
        metadata table as argument. A random file name prefix will
        automatically be added to your raw data.

-   **Trimming (and paired-end merging):**

    -   Navigate to /scripts/ and chose the appropriate trimming (and
        merging) script for your data. Consult the above paragraph and
        the scripts themselves for more information.

    -   Execute the script with the PREFIX of your dataset, the SAMPLE
        name and the SEQRUN as arguments, as stated in the
        metadata.txt (and the read length cutoff, if applicable).

    -   When done, navigate to the automatically created "processing"
        folder in /trimdata/ and view the log file(s) to confirm the
        trimming (and merging) was successful.

    -   Move the log file(s) into /trimdata/trimlogs/ and the FASTQ       
        file(s) one folder up into /trimdata/.

    -   Delete the "processing" folder, which should now be empty.

-   **Combining datasets (optional):**

    -   Navigate to /scripts/ and execute the cobine\_files.sh with the
        SAMPLE name, either ‘mappable’, ‘mappable_R1’ or
        ‘mappable_R2’ and all PREFIXes of the datasets you wish to
        combine as arguments.

    -   Navigate to the automatically created "processing" folder in
        /trimdata/.

    -   The file name of the output file lists all prefixes of files
        which were combined. When using this combined file for mapping
        or further combining, refer to it using the whole set of
        prefixes as one argument (e.g. pr1_pr2_pr3 for the file
        pr1_pr2_pr3+BearA_mappable.fastq.gz).

    -   Move the combined file from the "processing" folder one folder
        up into /trimdata/.

    -   Delete the "processing" folder, which should now be empty.

-   **Mapping**

    -   Navigate to /scripts/ and chose the appropriate mapping script
        for your data. Consult the above paragraph and the scripts
        themselves for more information.

    -   Execute the script. As arguments use the PREFIX of your dataset,
        the reference genome to be used for mapping, a three-character
        taxon identifier of your sample and the SAMPLE name.

    -   When done, navigate to the automatically created "processing"
        folder within the "mapped" folder corresponding to the
        reference used for mapping and view the log file to confirm
        the mapping was successful.

    -   Move the log file into the "logs" folder within the appropriate
        "mapped" folder.

    -   Move the .bam and .bam.bai files one folder up into their
        appropriate "mapped" folder.

    -   Delete the "processing" folder, which should now be empty.

-   **State of the Cave (optional)**

    -   Consult the /rawdata/metadata.txt for a detailed list of all rawdata in the BEARCAVE.

    -   For an overview of which datasets have been trimmed and mapped
        (and to which references), navigate to /scripts/ and execute the
        cavescan.sh without any arguments.

    -   This will generate a comma separated table (BEARCAVE/
        state_of_the_cave.txt) that lists all samples and prefixes
        associated with them and if they have been trimmed and/or mapped
        to the references or not.

**References**

1. Martin M. Cutadapt removes adapter sequences from high-throughput
sequencing reads. EMBnet.journal. 2011;17:10.

2. Magoč T, Salzberg SL. FLASH: fast length adjustment of short reads
to improve genome assemblies. Bioinformatics. academic.oup.com; 2011;27:
2957–2963.

3. Li H, Durbin R. Fast and accurate short read alignment with
Burrows–Wheeler transform. Bioinformatics. Oxford University Press;
2009;25: 1754–1760.

4. Li H, Handsaker B, Wysoker A, Fennell T, Ruan J, Homer N, et al. The
Sequence Alignment/Map format and SAMtools. Bioinformatics. 2009;25:
2078–2079.

5.  Dabney J, Knapp M, Glocke I, Gansauge M-T, Weihmann A, Nickel B, et al. Complete mitochondrial genome sequence of a Middle Pleistocene cave bear reconstructed from ultrashort DNA fragments. Proc Natl Acad Sci U S A. 2013;110: 15758–15763.

6. Gansauge M-T, Meyer M. Single-stranded DNA library preparation for
the sequencing of ancient or damaged DNA. Nat Protoc. 2013;8:
737–748.
