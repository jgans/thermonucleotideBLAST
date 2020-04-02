#ifndef __NCBI_UTIL
#define __NCBI_UTIL

#include <stdlib.h>

/*
* ===========================================================================
*
*                            PUBLIC DOMAIN NOTICE                          
*               National Center for Biotechnology Information
*                                                                          
*  This software/database is a "United States Government Work" under the   
*  terms of the United States Copyright Act.  It was written as part of    
*  the author's official duties as a United States Government employee and 
*  thus cannot be copyrighted.  This software/database is freely available 
*  to the public for use. The National Library of Medicine and the U.S.    
*  Government have not placed any restriction on its use or reproduction.  
*                                                                          
*  Although all reasonable efforts have been taken to ensure the accuracy  
*  and reliability of the software and data, the NLM and the U.S.          
*  Government do not and cannot warrant the performance or results that    
*  may be obtained by using this software or data. The NLM and the U.S.    
*  Government disclaim all warranties, express or implied, including       
*  warranties of performance, merchantability or fitness for any particular
*  purpose.                                                                
*                                                                          
*  Please cite the author in any work or product based on this material.   
*
* ===========================================================================
*/

#ifndef PNTR
#define PNTR    *
#endif

#ifndef PNTR
#define PNTR	*
#endif
#ifndef HNDL
#define HNDL	*
#endif

#ifndef FnPtr
typedef int		(*Nlm_FnPtr)(void);
#define FnPtr		Nlm_FnPtr
#endif

#ifndef VoidPtr
typedef void PNTR	Nlm_VoidPtr;
#define VoidPtr		Nlm_VoidPtr
#endif
#ifndef Pointer
#define Pointer		Nlm_VoidPtr
#endif

#ifndef Handle
typedef void HNDL	Nlm_Handle;
#define Handle		Nlm_Handle
#endif

#ifndef Char
typedef char		Nlm_Char, PNTR Nlm_CharPtr;
#define Char		Nlm_Char
#define CharPtr		Nlm_CharPtr
#endif

#ifndef Uchar
typedef unsigned char	Nlm_Uchar, PNTR Nlm_UcharPtr;
#define Uchar		Nlm_Uchar
#define UcharPtr	Nlm_UcharPtr
#endif

#ifndef Boolean
typedef unsigned char	Nlm_Boolean, PNTR Nlm_BoolPtr;
#define Boolean		Nlm_Boolean
#define BoolPtr		Nlm_BoolPtr
#endif

#ifndef Byte
typedef unsigned char	Nlm_Byte, PNTR Nlm_BytePtr;
#define Byte		Nlm_Byte
#define BytePtr		Nlm_BytePtr
#define BYTE_MAX	UCHAR_MAX
#endif

#ifndef Int1
typedef signed char	Nlm_Int1, PNTR Nlm_Int1Ptr;
#define Int1		Nlm_Int1
#define Int1Ptr		Nlm_Int1Ptr
#define INT1_MIN	SCHAR_MIN
#define INT1_MAX	SCHAR_MAX
#endif

#ifndef Uint1
typedef unsigned char	Nlm_Uint1, PNTR Nlm_Uint1Ptr;
#define Uint1		Nlm_Uint1
#define Uint1Ptr	Nlm_Uint1Ptr
#define UINT1_MAX	UCHAR_MAX
#endif

#ifndef Int2
typedef short		Nlm_Int2, PNTR Nlm_Int2Ptr;
#define Int2		Nlm_Int2
#define Int2Ptr		Nlm_Int2Ptr
#define INT2_MIN	SHRT_MIN
#define INT2_MAX	SHRT_MAX
#endif

#ifndef Uint2
typedef unsigned short	Nlm_Uint2, PNTR Nlm_Uint2Ptr;
#define Uint2		Nlm_Uint2
#define Uint2Ptr	Nlm_Uint2Ptr
#define UINT2_MAX	USHRT_MAX
#endif

#ifndef Int4
typedef signed int  Nlm_Int4, PNTR Nlm_Int4Ptr;
#define Int4        Nlm_Int4
#define Int4Ptr     Nlm_Int4Ptr
#define INT4_MIN    (-2147483647-1)
#define INT4_MAX    2147483647
#endif

#ifndef Uint4
typedef unsigned int  Nlm_Uint4, PNTR Nlm_Uint4Ptr;
#define Uint4         Nlm_Uint4
#define Uint4Ptr      Nlm_Uint4Ptr
#define UINT4_MAX     4294967295U
#endif

#ifndef TIME_MAX
#define TIME_MAX	ULONG_MAX
#endif

#ifndef FloatLo
typedef float		Nlm_FloatLo, PNTR Nlm_FloatLoPtr;
#define FloatLo		Nlm_FloatLo
#define FloatLoPtr	Nlm_FloatLoPtr
#endif

#ifndef FloatHi
typedef double		Nlm_FloatHi, PNTR Nlm_FloatHiPtr;
#define FloatHi		Nlm_FloatHi
#define FloatHiPtr	Nlm_FloatHiPtr
#endif

#ifndef BigScalar
#define BigScalar long
#endif

#ifndef FAR
#define FAR
#endif

#define LIBCALL	FAR PASCAL

#define Int8  long
typedef Int8   Nlm_Int8;

typedef union dataval {
	Nlm_VoidPtr ptrvalue;
	Nlm_Int4 intvalue;
	Nlm_FloatHi realvalue;
	Nlm_Boolean boolvalue;
	Nlm_FnPtr	funcvalue;
	Nlm_Int8    bigintvalue;
}	DataVal, PNTR DataValPtr;

typedef struct valnode {
	Nlm_Uint1 choice;          /* to pick a choice */
	Nlm_Uint1 extended;        /* extra fields reserved to NCBI allocated in structure */
	DataVal data;              /* attached data */
	struct valnode PNTR next;  /* next in linked list */
} ValNode, PNTR ValNodePtr;

typedef ValNode SeqId, FAR *SeqIdPtr;

typedef struct textseqid {
	CharPtr name,
		accession,
		release;
	Int2 version;             /* INT2_MIN (ncbilcl.h) = not set */
} TextSeqId, PNTR TextSeqIdPtr;

typedef struct objid {
	Int4 id;
	CharPtr str;
} ObjectId, PNTR ObjectIdPtr;

typedef struct giim {
	Int4 id;
	CharPtr db,
		release;
} Giim, PNTR GiimPtr;

typedef struct idpat {
	CharPtr country,
		number,                          /** actually CHOICE of number or app_number */
		app_number;
	CharPtr doc_type;
} IdPat, PNTR IdPatPtr;

typedef struct patentseqid {
	Int4 seqid;
	IdPatPtr cit;
} PatentSeqId, PNTR PatentSeqIdPtr;

#define DatePtr 		NCBI_DatePtr
#define StringICmp		Nlm_StringICmp

typedef struct date {
	Uint1 data[8];      /* see box above */
	CharPtr str;            /* str or season or NULL */
} NCBI_Date, PNTR NCBI_DatePtr;

typedef struct pdbseqid {
	CharPtr mol;
	Uint1 chain;        /* 0 = no chain set.  default = 32 */
	DatePtr rel;
} PDBSeqId, PNTR PDBSeqIdPtr;

typedef struct dbtag {
	CharPtr db;
	ObjectIdPtr tag;
} Dbtag, PNTR DbtagPtr;

#define NUM_SEQID 20     /* total number of SeqId types */

#define SEQID_NOT_SET ( (Uint1)0)
#define SEQID_LOCAL ( (Uint1)1)
#define SEQID_GIBBSQ ( (Uint1)2)
#define SEQID_GIBBMT ( (Uint1)3)
#define SEQID_GIIM ( (Uint1)4)

#define SEQID_GENBANK ( (Uint1)5)
#define SEQID_EMBL ( (Uint1)6)
#define SEQID_PIR ( (Uint1)7)
#define SEQID_SWISSPROT ( (Uint1)8)


#define SEQID_PATENT ( (Uint1)9)
#define SEQID_OTHER ( (Uint1)10)
#define SEQID_GENERAL ( (Uint1)11)
#define SEQID_GI ( (Uint1)12)
#define SEQID_DDBJ ((Uint1)13)
#define SEQID_PRF ((Uint1)14)
#define SEQID_PDB ((Uint1)15)

#define SEQID_TPG ((Uint1)16)
#define SEQID_TPE ((Uint1)17)
#define SEQID_TPD ((Uint1)18)

#define SEQID_GPIPE  ((Uint1)19)

#define PRINTID_FASTA_SHORT ( (Uint1)1)
#define PRINTID_FASTA_LONG ( (Uint1)2)
#define PRINTID_TEXTID_LOCUS ( (Uint1)3)
#define PRINTID_TEXTID_ACCESSION ( (Uint1)4)
#define PRINTID_TEXTID_ACC_VER ( (Uint1)5)
#define PRINTID_TEXTID_ACC_ONLY ( (Uint1)6)
#define PRINTID_REPORT ( (Uint1)7)
#define PRINTID_FASTA_GENERAL ( (Uint1)8)
#define PRINTID_FASTA_ALL ( (Uint1)9)

#define SIC_DIFF 1
#define SIC_NO 0
#define SIC_YES 2

#define Seq_strand_unknown 0
#define Seq_strand_plus 1
#define Seq_strand_minus 2
#define Seq_strand_both 3
#define Seq_strand_both_rev 4
#define Seq_strand_other 255

enum _ErrSev { SEV_NONE=0, SEV_INFO, SEV_WARNING, SEV_ERROR, SEV_REJECT, SEV_FATAL, SEV_MAX };
typedef enum _ErrSev ErrSev;

SeqIdPtr SeqIdSetFree(SeqIdPtr sip);
SeqIdPtr SeqIdSetDup(SeqIdPtr seqid);

Uint1 SeqIdComp(SeqIdPtr a, SeqIdPtr b);

CharPtr SeqIdPrint(SeqIdPtr isip, CharPtr buf, Uint1 format);
CharPtr SeqIdWrite(SeqIdPtr isip, CharPtr buf, Uint1 format, Uint4 buflen);
SeqIdPtr SeqIdParse(CharPtr buf);
SeqIdPtr SeqIdDup(SeqIdPtr oldid);

ValNodePtr ValNodeNew(ValNodePtr vnp);
SeqIdPtr SeqIdFree(SeqIdPtr anp);

Boolean ObjectIdMatch(ObjectIdPtr a, ObjectIdPtr b);
Boolean IdPatMatch(IdPatPtr a, IdPatPtr b);
int Nlm_StringICmp(const char FAR *a, const char FAR *b);
int Nlm_StrICmp(const char FAR *a, const char FAR *b);
Boolean DbtagMatch(DbtagPtr a, DbtagPtr b);

ValNodePtr ValNodeCopyStr(ValNodePtr PNTR head, Nlm_Int2 choice, Nlm_CharPtr str);
ValNodePtr ValNodeAdd(ValNodePtr PNTR head);

#ifndef THIS_FILE
#define THIS_FILE  __FILE__
#endif
#ifndef THIS_MODULE
#define THIS_MODULE  ""
#endif

#define DBFLAG 0

#ifndef TRUE
/** bool replacment for C indicating true. */
#define TRUE 1
#endif

#ifndef FALSE
/** bool replacment for C indicating false. */
#define FALSE 0
#endif

#define StringHasNoText 	Nlm_StringHasNoText
#define StringLen       	Nlm_StringLen
#define StrLen   		Nlm_StrLen
#define Nlm_StrLen		strlen

Nlm_Boolean Nlm_StringHasNoText (const char FAR *str);
size_t Nlm_StringLen(const char *str);

#define MemNew(x)     Nlm_MemNew(x)

char* Nlm_MemNew(size_t size);

#define	StringCat		Nlm_StringCat
Nlm_CharPtr Nlm_StringCat(char FAR *to, const char FAR *from);

ValNodePtr ValNodeFreeData(ValNodePtr vnp);

#define	MemFree		Nlm_MemFree

void* Nlm_MemFree(void *ptr);

#define StringNCpy_0    Nlm_StringNCpy_0

Nlm_CharPtr Nlm_StringNCpy_0 (char FAR *to, const char FAR *from, size_t max);

Nlm_Uint4 Nlm_LabelCopy (Nlm_CharPtr to, const char* from, Nlm_Uint4 buflen);
void Nlm_LabelCopyNext(Nlm_CharPtr PNTR to, const char* from, Nlm_Uint4 PNTR buflen);

SeqIdPtr SeqIdSelect(SeqIdPtr sip, Uint1Ptr order, Int2 num);

#define SHOWVERSION 1

#define IS_DIGIT(c)     ('0'<=(c) && (c)<='9')
#define IS_UPPER(c)     ('A'<=(c) && (c)<='Z')
#define IS_LOWER(c)     ('a'<=(c) && (c)<='z')

#define TO_LOWER(c)     ((Nlm_Char)(IS_UPPER(c) ? (c)+' ' : (c)))
#define TO_UPPER(c)     ((Nlm_Char)(IS_LOWER(c) ? (c)-' ' : (c)))

#define StringCmp       Nlm_StringCmp
int Nlm_StringCmp(const char FAR *a, const char FAR *b);

ObjectIdPtr ObjectIdNew (void);

#define StringSave      Nlm_StringSave
Nlm_CharPtr Nlm_StringSave (const char FAR *from);
Nlm_CharPtr Nlm_StrSave (const char FAR *from);

GiimPtr GiimNew (void);
TextSeqIdPtr TextSeqIdNew (void);
PatentSeqIdPtr PatentSeqIdNew (void);
IdPatPtr IdPatNew (void);
DbtagPtr DbtagNew (void);
PDBSeqIdPtr PDBSeqIdNew (void);

#define StringNCpy      Nlm_StringNCpy

Nlm_CharPtr Nlm_StringNCpy (char FAR *to, const char FAR *from, size_t max);
Nlm_CharPtr Nlm_ClearDestString(Nlm_CharPtr to, size_t max);

ValNodePtr ValNodeFree (ValNodePtr vnp);

ObjectIdPtr ObjectIdDup(ObjectIdPtr oldid);
DbtagPtr DbtagDup (DbtagPtr oldtag);
NCBI_DatePtr DateDup(NCBI_DatePtr dp);
ObjectIdPtr ObjectIdFree (ObjectIdPtr oid);

#define	MemCopy	Nlm_MemCopy

void * Nlm_MemCopy(void *dst, const void *src, size_t bytes);

GiimPtr GiimFree(GiimPtr gip);
TextSeqIdPtr TextSeqIdFree(TextSeqIdPtr tsip);
PatentSeqIdPtr PatentSeqIdFree(PatentSeqIdPtr psip);

DbtagPtr DbtagFree(DbtagPtr dbt);
PDBSeqIdPtr PDBSeqIdFree(PDBSeqIdPtr pdbsip);

#define StrNCat  strncat

NCBI_DatePtr DateNew (void);
IdPatPtr IdPatFree(IdPatPtr idp);
NCBI_DatePtr DateFree (NCBI_DatePtr dp);

#endif // __NCBI_UTIL
