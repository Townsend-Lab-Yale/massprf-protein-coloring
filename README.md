# massprf-protein-coloring
maps the output of massprf into a file that can be used in Chimera to color proteins

## example calls
genMASSPRF_Chimera(pdbFile = "~/../Downloads/SARS-CoV-1_ORF7a_Gene.pdb",MASSPRF_Nuc_Fasta = "orf7a_gene.fasta",MASSPRF_Table = "orf7a_MASS_PRF.tsv",scaling=1)


genMASSPRF_Chimera(pdbFile = "S_23.pdb",MASSPRF_Nuc_Fasta = "s_gene.fasta",onlySig=F,MASSPRF_Table = "~/../Downloads/S1_MASSPRF.tsv",scaling=6,outfile = "S_cc.txt")


## Inputs
pdbFile: The protein structure you are trying to map onto, in pdb format. REQUIRED INPUT


MASSPRF_Nuc_Fasta: the nucleotide sequence representing your MASSPRF input, in open reading frame and in fasta format. REQUIRED INPUT


MASSPRF_Table: The output table from MASSPRF containing the gamma values and CI. REQUIRED INPUT


outfile: where the resultant commands for Chimera to color will be stored. Defaults chimeraColoring.txt

    (hint: execute the following on the Chimera command line): read <outfile>; 
    
    
scaling: the scale factor from MASSPRF--look in the output of MASSPRF to find this. 1 or a multiple of 3. REQUIRED INPUT


sigSetting: what are you using to determine significance? only affects scale factor 1; otherwise sitewise significance is used. Defaults average


onlySig: T or F do you want to only color significant sites. Defaults T


rgb1 and rgb2: vectors of size 3, corresponding to a rgb value for coloring the selection intensity. Defaults red and blue


bins: how many equally spaced apart color categories you want to use for your data. Defaults 100


ehColor: what do you want to color 'eh' residues (missing data, nonsignificant) as a vector size three of rgb. Defaults grey
