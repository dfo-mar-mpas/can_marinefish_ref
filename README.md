# Environmental DNA reference database for Canadian Marine Fish species using 12s and 16s genes.
<img src="inst/hexlogo-01.png" align="left" width="200px"/>

Here are two custom DNA reference libraries for marine fish (_Actinopterygii_ only) in Canada (Pacific and Atlantic Oceans) for the 12S and 16S amplicons described in [He et al., 2022](https://cdnsciencepub.com/doi/10.1139/cjfas-2021-0215). 

<br clear="left"/>

These are reference libraries specifically constructed for amplicons of the following markers:

1.	12S modified MiFish ([Miya et al., 2015](https://royalsocietypublishing.org/doi/10.1098/rsos.150088)). 
    Forward primer = 5’- **CGTGCCAGCCACCGCGGTT** -3’ 
    Reverse primer = 5’- **CATAGTGGGGTATCTAATCCCAGTTTG** -3’ 

2.	16S Fish ([McInnes et al., 2017](https://www.frontiersin.org/articles/10.3389/fmars.2017.00277/full)). 
    Forward primer = 5’- **AGCGYAATCACTTGTCTYTTAA** -3’ 
    Reverse primer = 5’- **CRBGGTCGCCCCAACCRAA** -3’
    
They are formatted for use with the FuzzyID2 software package for taxonomic assignment ([Shi et al., 2018](https://onlinelibrary.wiley.com/doi/10.1111/1755-0998.12738)).


## File Descriptions:
1. **12S_reference_library_Actinopterygii.fasta**: is a fasta formatted file with headers formatted for use with FuzzyID2 software in the format of GBaccession_Family_Genus_Species. *Note that original Genbank Accession numbers are maintained with the exceptions detailed for shared haplotypes between species.

2. **12S_haplogroups_list.csv**: is a csv formatted file containing the list of haplotypes shared between species and genera. The header is formatted as (1) accession number change to group initials (2), followed by group number and 4 zeros and individual number. The rest of the header contains the Group Name for Family (i.e., Agonidae1) followed by the group code (i.e., AG1) followed by a unique species identifier that is sequential for all groups (i.e., species1). Each haplotype in a group will have the exact same Family, Genus and Species text to avoid any program error. For example, each entry for Agonidae1 group is AG100001_Agonidae1_AG1_species1, AG100002_Agonidae1_AG1_species1, and AG100003_Agonidae1_AG1_species1.

3. **16S_reference_library_Actinopterygii.fasta**: is a fasta formatted file with headers formatted for use with FuzzyID2 software in the format of GBaccession_Family_Genus_Species. *Note that original Genbank Accession numbers are maintained with the exceptions detailed for shared haplotypes between species.

4. **16S_haplogroups_list.csv**: is a csv formatted file containing the list of haplotypes shared between species and genera. The header is formatted as (1) accession number change to group initials (2), followed by group number and 4 zeros and individual number. The rest of the header contains the Group Name for Family (i.e. Agonidae1) followed by the group code (i.e. AG1) followed by a unique species identifier that is sequential for all groups (i.e. species1). Each haplotype in a group will have the exact same Family, Genus and Species text to avoid any program error. For example, each entry for Agonidae1 group is AG100001_Agonidae1_AG1_species1, AG100002_Agonidae1_AG1_species1, and AG100003_Agonidae1_AG1_species1.

5. **Coad_OBIS_fish_list_Canada.csv**: the final species list used as a template to gather region specific reference DNA sequences -> all species that spend at least part of their life cycle in marine or brackish waters in the Atlantic and Pacific Oceans within Canadian waters.


## The general overview of reference library construction

  1. [Determine species list for _Actinopterygii_ for Canadian marine waters in Atlantic and Pacific Oceans.](#step1)
  2. [Gather GenBank entries for genes and species of choice.](#step2)
  3. [Perform in silico PCR on all entries -> reflib1.](#step3)
  4. [Determine which entries failed in silico PCR and manually align + reflib1 -> reflib2.](#step4)
  5. [Identify unique haplotypes and collapse entries + reflib2 -> reflib3.](#step5)
  6. [Calculate 95% confidence intervals of intraspecific distances and generate list of GenBank accession numbers greater than cut-off.](#step6)
  7. [Visually inspect potential GenBank ID errors using phylogenetic trees and remove entries from reference library + reflib3 -> reflib4.](#step7)

Specific methods for reference library construction under numbered headings of general steps above.

**1. Determine species list** <a name="step1"/> 

The final list of marine species in Canada was comprised of the list derived from [OBIS](https://obis.org/) observations from within Canada and supplemented by a list of Canadian marine fish species compiled by Ichthyologist Brian Coad [website](http://www.briancoad.com/main.asp?page=whatsnew.asp). This list totals 1543 species in Actinopterygii and is available as [**Coad_OBIS_fish_list_Canada.csv**](https://github.com/dfo-mar-mpas/can_marinefish_ref/blob/main/data/Coad_OBIS_fish_list_Canada.csv).

**2. Gather Genbank entries using NCBI E-Utilities** <a name="step2"/>

  2.1 _‘esearch’_ and _‘efetch’_ species names and gene names according to list. Include or do separate searches for alternate naming schemes (e.g., small/large ribosomal sub-unit). Specify GenBank format for downloads and do not download RefSeq entries as these are duplicates.(`<br>`)
  2.2 Note missing species by looking at files with zero size.(`<br>`)
  2.3 Concatenate all files into single GenBank formatted file.(`<br>`)


**Perform in _silico_ PCR** <a name="step3"/>
  1. Convert to obitools database
      - obiconvert -d /home/kristen/Documents/ncbi20201008/ --ecopcrdb-output=16S.db /home/kristen/Documents/barcoding_gap/16S_species/REFLIB_PIPELINE/species_gb/all.gb
  2. Run ecoPCR for specific primers with the following flags (e.g. 16S given)
      - ecoPCR -d 16S.db -e 3 -l 100 -L 300 AGCGYAATCACTTGTCTYTTAA CRBGGTCGCCCCAACCRAA > 16S_ecopcr_out_raw.txt
  3. Remove first 13 lines of output file. 
  4. Reformat into FuzzyID2 reference library fasta format -> reflib1

**Determine which entries failed _in silico_ PCR** <a name="step4"/>
  1. work from file in reflib directory.
  2. obtain list of gb accession numbers in reflib1 
      - sed -n '/^>/p' 16S_reflib1.fasta > accession_list_reflib1.txt
      - reformat in excel to just accession numbers, sort in excel.
  3. Obtain list of gb accession numbers in all.gb (raw gb entries). 
      - sed -n '/^ACCESSION/p' all.gb > accession_list_all.txt
      - reformat in excel to just accession numbers, sort in excel.
  4. Determine which accession numbers were cut by in silico PCR
      - comm -23 <(sort accession_list_all.csv | uniq) <(sort accession_list_reflib1.csv | uniq) > accession_list_not_ecopcr.txt
  5. use list of accession numbers in esearch and efetch commands, this time download fasta format.
  6. Make reflib1 into a blast database and blast the extra sequences.
  7. Steps to make reflib1 into database:
      - extract species names from each entry in reflib1 and look up taxids at [website](https://www.ncbi.nlm.nih.gov/Taxonomy/TaxIdentifier/tax_identifier.cgi)
      - truncate headers in reflib1 to 50 characters
      merge the reflib1 truncated header file with the taxid file to create a taxid map for making the local blast db
      - makeblastdb -in 16S_reflib1_headertrunc50.fasta -parse_seqids -blastdb_version 5 -taxid_map 16S_taxid_map.csv -title "16Sreflib1" -dbtype nt
      - blast the remaining sequences to reflib1 to see which sequences are a potential match, output is a list of query accession numbers that had a match. 
      -blastn -db 16S_reflib1_headertrunc50.fasta -query not_ecopcr_oneliner2.fasta -evalue 1e-6 -outfmt '6 qseqid' -max_target_seqs 1 > blast_out.txt.
  8. Take output and acquire the full header and sequence from the original query file not_ecopcr_oneliner2.fasta using awk batch cmd.
  9. Separate output into batches and align manually.
  10. Format trimmed sequences for fuzzyid2:
      - print headers to file
  sed -n '/^>/p' blast_outpu_all_trimmed.fasta > blast_output_all_headers.txt
      - delete headers to make a sequence file
  sed -e '/^>/d' blast_output_all_trimmed.fasta > blast_output_all_trimmed_seq.txt
      - parse headers in excel to fuzzyid2 format by deleting all information except accession number and species name. Add in column for family name and export to text file
      - add in new headers to sequences.
  paste -d \\n blast_output_all_trimmed_formattedheader.txt blast_output_all_trimmed_seq.txt > blast_output_all_formattedfuzzyid2.fasta
      - add in trimmed sequences to reflib1 -> reflib2.
      - remove short sequences, those sequences that are obviously cut-off and not natural length variation. If there are conspecifics with much larger sequences then remove the entry with the sequence that was cut off.


**Identify unique haplotypes and collapse entries** <a name="step5"/>
  1. Align reflib2_Acitnopterygii using mafft and load into R.
  2. Execute steps in reflib2haplos.R **Use genetic distances for discovering unique haplotypes, do not use identities!
  3. Go through list of unique haplos and record changes to reflib2. The rules for collapsing unique haplotypes are: maximum of three entries for each unique haplotype per species. When a haplotype is shared between 2 or more species, record those species and form a group. Name group and change accession number to group initials (2), followed by group number and 4 zeros and individual number. The rest of the header contains the Group Name for Family (i.e., _Agonidae1_) followed by the group code (i.e., _AG1_) followed by a unique species identifier that is sequential for all groups (i.e., _species1_). Each haplotype in a group will have the exact same Family, Genus and Species text to avoid any program error. For example, each entry for Agonidae1 group is AG100001_Agonidae1_AG1_species1, AG100002_Agonidae1_AG1_species1, and AG100003_Agonidae1_AG1_species1.
  4. Integrate new header file into edited reflib2 -> save as reflib3.
  5. Groups with 2 or more species and/or genera are potential Genbank misidentifications. 
  6. Parse all entries per group's Family in reflib2, create phylogeny in R, and inspect. (NB. do not used collapsed library to parse Family entries).

**Calculate 95% cut-off value of intraspecific distances and generate list of GenBank accession numbers greater than cut-off.** <a name="step6"/>
1. Generate list of unique species from reflib3 using sed and excel
2. Parse reflib3 for each species and create one new file for each species using awk.
3. Align each species file separately using mafft.
4. Analyze each aligned species file in R, calculating intraspecific K2P distances min, mean, max, and stdev. Write to file.
5. Calculate the grand mean of means and the average sd. Calculate 95% cut-off as 4.5 average sd* of grand mean. 
(*Chebyshev’s inequality was used to determine the 95% confidence interval as the distribution of average pairwise intraspecific variation was heavily skewed towards zero).
a. Calculate the average sd using the formula: Average S.D. = √ ((n<sub>1</sub>-1)s<sub>1</sub><sup>2</sup> +  (n<sub>2</sub>-1)s<sub>2</sub><sup>k</sup> + … +  (n<sub>k</sub>-1)s<sub>k</sub><sup>2</sup>) /  (n<sub>1</sub>+n<sub>2</sub> + … + n<sub>k</sub> – k) where nk: Sample size for kth group, sk: Standard deviation for kth group, and k: Total number of groups
6. Highlight each species with a max value greater than the 95% cut-off and generate a Family phylogeny.

**Visually inspect potential GenBank ID errors using phylogenetic trees and remove entries from reference library.** <a name="step7"/>

For each species that was identified outside of the range, a NJ tree for all sequences within the Family was generated and individual pairwise distance matrices were examined. Obvious outliers were identified as (1) single Genbank entries that were placed outside of a monophyletic species clade, (2) single Genbank entries that had genetic distances with all other conspecifics above the cutoff. In these cases, the single Genbank entry for that species was removed. Once removed, each species fell within the 95% confidence interval. Less obvious outliers occurred in species with large distributions that may represent phylogeographic variation. The entire species was removed in these cases as we could not reasonably assume the cause of high intraspecific distance for any specific Genbank record. 

Unique haplotypes that were shared among species, genera, and families were also examined. All haplotypes shared among families were examined for outliers, which were obvious in all cases and those entries removed so each unique haplotype was shared at or below the genus level. Most cases where haplotypes were shared among species or genera were treated as groups where species level could not be resolved. This approach was favoured over including all species in an overall genus group as haplotype sharing only applied to a subset of species within a genus, offering the infest level of taxonomic assignment possible. 

