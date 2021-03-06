// ThermonucleotideBLAST
// 
// Copyright (c) 2008, Los Alamos National Security, LLC
// All rights reserved.
// 
// Copyright 2007. Los Alamos National Security, LLC. This software was produced under U.S. Government 
// contract DE-AC52-06NA25396 for Los Alamos National Laboratory (LANL), which is operated by Los Alamos 
// National Security, LLC for the U.S. Department of Energy. The U.S. Government has rights to use, 
// reproduce, and distribute this software.  NEITHER THE GOVERNMENT NOR LOS ALAMOS NATIONAL SECURITY, 
// LLC MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.  
// If software is modified to produce derivative works, such modified software should be clearly marked, 
// so as not to confuse it with the version available from LANL.
// 
// Additionally, redistribution and use in source and binary forms, with or without modification, 
// are permitted provided that the following conditions are met:
// 
//      * Redistributions of source code must retain the above copyright notice, this list of conditions 
//        and the following disclaimer.
//      * Redistributions in binary form must reproduce the above copyright notice, this list of conditions 
//        and the following disclaimer in the documentation and/or other materials provided with the distribution.
//      * Neither the name of Los Alamos National Security, LLC, Los Alamos National Laboratory, LANL, 
//        the U.S. Government, nor the names of its contributors may be used to endorse or promote products 
//        derived from this software without specific prior written permission.
// 
// 
// THIS SOFTWARE IS PROVIDED BY LOS ALAMOS NATIONAL SECURITY, LLC AND CONTRIBUTORS "AS IS" AND ANY 
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY 
// AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LOS ALAMOS NATIONAL SECURITY, LLC 
// OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, 
// OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, 
// STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS 
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

enum {
	SEQUENCE_ONTOLOGY,
	REGION,
	INTERIOR_CODING_EXON,
	PCR_PRODUCT,
	READ_PAIR,
	SCRNA,
	MATCH_SET,
	MATCH_PART,
	GENE_PART,
	OPERATOR,
	TRANSPOSABLE_ELEMENT,
	EXPRESSED_SEQUENCE_MATCH,
	CLONE_INSERT_END,
	POLYPEPTIDE,
	SEQUENCE_VARIANT_OBS,
	SEQUENCE_FEATURE,
	PRIMER,
	PROVIRAL_REGION,
	METHYLATED_C,
	PROTEIN_CODING_PRIMARY_TRANSCRIPT,
	RIBOSOME_ENTRY_SITE,
	ATTENUATOR,
	TERMINATOR,
	ASSEMBLY_COMPONENT,
	EXON,
	SUPERCONTIG,
	CONTIG,
	READ,
	CLONE,
	DELETION,
	METHYLATED_A,
	SPLICE_SITE,
	FIVE_PRIME_SPLICE_SITE,
	THREE_PRIME_SPLICE_SITE,
	ENHANCER,
	PROMOTER,
	CROSS_GENOME_MATCH,
	OPERON,
	CLONE_INSERT_START,
	TRANSLATED_NUCLEOTIDE_MATCH,
	NON_TRANSCRIBED_REGION,
	PRIMARY_TRANSCRIPT,
	REPEAT_FAMILY,
	INTRON,
	RFLP_FRAGMENT,
	FIVE_PRIME_EXON_CODING_REGION,
	THREE_PRIME_EXON_CODING_REGION,
	NONCODING_EXON,
	UTR,
	FIVE_PRIME_UTR,
	THREE_PRIME_UTR,
	PROCESSED_TRANSCRIPT,
	MRNA,
	TF_BINDING_SITE,
	ORF,
	FLANKING_REGION,
	RRNA,
	TRNA,
	SNRNA,
	SNORNA,
	MIRNA,
	MICROSATELLITE,
	INVERTED_REPEAT,
	ORIGIN_OF_REPLICATION,
	CLIP,
	MODIFIED_BASE_SITE,
	METHYLATED_BASE_FEATURE,
	CPG_ISLAND,
	DIRECT_REPEAT,
	TRANSCRIPTION_START_SITE,
	CDS,
	START_CODON,
	STOP_CODON,
	TAG,
	SAGE_TAG,
	CONSERVED_REGION,
	STS,
	CODING_CONSERVED_REGION,
	EXON_JUNCTION,
	NC_CONSERVED_REGION,
	PSEUDOGENE,
	RNAI_REAGENT,
	CHROMOSOME,
	CHROMOSOME_BAND,
	MATCH,
	SPLICE_ENHANCER,
	EST,
	NUCLEOTIDE_MATCH,
	PROTEIN_MATCH,
	ASSEMBLY,
	CODON,
	INSERTION_SITE,
	TRANSPOSABLE_ELEMENT_INSERTION_SITE,
	SMALL_REGULATORY_NCRNA,
	ENZYMATIC_RNA,
	RIBOZYME,
	RRNA_5_8S,
	HAMMERHEAD_RIBOZYME,
	RNASE_MRP_RNA,
	RNASE_P_RNA,
	TELOMERASE_RNA,
	U1_SNRNA,
	U2_SNRNA,
	U4_SNRNA,
	U4ATAC_SNRNA,
	U5_SNRNA,
	U6_SNRNA,
	U6ATAC_SNRNA,
	U11_SNRNA,
	U12_SNRNA,
	U14_SNORNA,
	VAULT_RNA,
	Y_RNA,
	RRNA_18S,
	BINDING_SITE,
	RESTRICTION_FRAGMENT,
	SEQUENCE_DIFFERENCE,
	SIGNAL_PEPTIDE,
	MATURE_PROTEIN_REGION,
	ARS,
	RASIRNA,
	PSEUDOGENIC_REGION,
	DECAYED_EXON,
	GOLDEN_PATH_FRAGMENT,
	TILING_PATH,
	TILING_PATH_FRAGMENT,
	NC_PRIMARY_TRANSCRIPT,
	THREE_PRIME_CODING_EXON_NONCODING_REGION,
	FIVE_PRIME_CODING_EXON_NONCODING_REGION,
	VIRTUAL_SEQUENCE,
	TRANSCRIBED_REGION,
	POLYA_SIGNAL_SEQUENCE,
	POLYA_SITE,
	CENTROMERE,
	CAP,
	GROUP_I_INTRON,
	AUTOCATALYTICALLY_SPLICED_INTRON,
	SRP_RNA,
	GUIDE_RNA,
	GROUP_II_INTRON,
	INTERGENIC_REGION,
	POLYA_SEQUENCE,
	BRANCH_SITE,
	POLYPYRIMIDINE_TRACT,
	TRANSCRIPTION_END_SITE,
	TELOMERE,
	SILENCER,
	INSULATOR,
	CHROMOSOMAL_STRUCTURAL_ELEMENT,
	MINISATELLITE,
	ANTISENSE_RNA,
	ANTISENSE_PRIMARY_TRANSCRIPT,
	SIRNA,
	STRNA,
	SMALL_SUBUNIT_RRNA,
	LARGE_SUBUNIT_RRNA,
	RRNA_5S,
	RRNA_28S,
	NCRNA,
	REPEAT_REGION,
	DISPERSED_REPEAT,
	SPLICEOSOMAL_INTRON,
	INSERTION,
	EST_MATCH,
	TRANSCRIPT,
	NUCLEASE_SENSITIVE_SITE,
	DELETION_JUNCTION,
	GOLDEN_PATH,
	CDNA_MATCH,
	SNP,
	REAGENT,
	OLIGO,
	JUNCTION,
	REMARK,
	POSSIBLE_BASE_CALL_ERROR,
	POSSIBLE_ASSEMBLY_ERROR,
	EXPERIMENTAL_RESULT_REGION,
	GENE,
	TANDEM_REPEAT,
	TRANS_SPLICE_ACCEPTOR_SITE,
	NUCLEOTIDE_MOTIF,
	READING_FRAME,
	ULTRACONTIG,
	ORIT,
	TRANSIT_PEPTIDE,
	GAP,
	GENE_GROUP_REGULATORY_REGION,
	PSEUDOGENIC_RRNA,
	PSEUDOGENIC_TRNA,
	CHROMOSOME_PART,
	GENE_MEMBER_REGION,
	TRANSCRIPT_REGION,
	PROCESSED_TRANSCRIPT_REGION,
	PRIMARY_TRANSCRIPT_REGION,
	MRNA_REGION,
	UTR_REGION,
	POLYPEPTIDE_REGION,
	SPLICEOSOMAL_INTRON_REGION,
	GENE_COMPONENT_REGION,
	CDS_REGION,
	RRNA_16S,
	RRNA_23S,
	RRNA_25S,
	REGULATORY_REGION,
	GENE_GROUP,
	SUBSTITUTION,
	COMPLEX_SUBSTITUTION,
	POINT_MUTATION,
	INVERSION,
	REGULON,
	DATABANK_ENTRY,

	// Terms not in the sofa file (but specified in the GFF3 format)
	NUCLEOTIDE_TO_PROTEIN_MATCH,
	UNKNOWN
};
