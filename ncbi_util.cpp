#include "ncbi_util.h"

#include <iostream>

#include <stdio.h>
#include <string.h>
#include <limits.h>

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

using namespace std;

#ifndef USE_NCBI

SeqIdPtr SeqIdSetFree (SeqIdPtr sip)
{
    SeqIdPtr next;

    while (sip != NULL)
    {
        next = sip->next;
        SeqIdFree(sip);
        sip = next;
    }
    return sip;
}

SeqIdPtr SeqIdSetDup(SeqIdPtr seqid)
{
   SeqIdPtr sid_head, sid, seqid_var;

   if (seqid == NULL)
      return seqid;
   else {
      seqid_var = seqid;
      sid_head = sid = SeqIdDup(seqid);
   }

   while ((seqid_var = seqid_var->next) != NULL) {
      sid->next = SeqIdDup(seqid_var);
      sid = sid->next;
   }

   return sid_head;
}

Uint1 SeqIdComp(SeqIdPtr a, SeqIdPtr b)
{
    Uint1 choice;
    TextSeqIdPtr at, bt;

    if ((a == NULL) || (b == NULL))
        return SIC_DIFF;

	choice = a->choice;
	if (choice != b->choice)
	{
		switch (choice)
		{
			case SEQID_GENBANK:          /* these could be confused */
			case SEQID_EMBL:
			case SEQID_DDBJ:
			case SEQID_TPG:
			case SEQID_TPE:
			case SEQID_TPD:
		    case SEQID_GPIPE:
				switch (b->choice)
				{
					case SEQID_GENBANK:   /* its ok */
					case SEQID_EMBL:
					case SEQID_DDBJ:
					case SEQID_TPG:
					case SEQID_TPE:
					case SEQID_TPD:
				    case SEQID_GPIPE:
						break;  
					default:
						return SIC_DIFF;
				}
				break;
			default:
				return SIC_DIFF;
		}
	}

    switch (choice)
    {
        case SEQID_NOT_SET:
            return SIC_DIFF;
        case SEQID_LOCAL:   
            if (ObjectIdMatch((ObjectIdPtr)a->data.ptrvalue, (ObjectIdPtr)b->data.ptrvalue))
				return SIC_YES;
			else
				return SIC_NO;
        case SEQID_GIBBSQ:   /* gibbsq */
        case SEQID_GIBBMT:   /* gibbmt */
        case SEQID_GI:  /* gi */
            if (a->data.intvalue == b->data.intvalue)
				return SIC_YES;
			else
				return SIC_NO;
        case SEQID_GIIM:   /* giim */
            if (((GiimPtr)a->data.ptrvalue)->id == ((GiimPtr)b->data.ptrvalue)->id)
				return SIC_YES;
			else
				return SIC_NO;
        case SEQID_PATENT:   /* patent seq */
            if (((PatentSeqIdPtr)a->data.ptrvalue)->seqid !=
                ((PatentSeqIdPtr)b->data.ptrvalue)->seqid)
                return SIC_NO;
            if (IdPatMatch(((PatentSeqIdPtr)a->data.ptrvalue)->cit,
                ((PatentSeqIdPtr)b->data.ptrvalue)->cit))
				return SIC_YES;
			else
				return SIC_NO;
		case SEQID_PDB:     /* pdb */
            if ( StringICmp(((PDBSeqIdPtr)a->data.ptrvalue)->mol,
                ((PDBSeqIdPtr)b->data.ptrvalue)->mol))
                return SIC_NO;
            /*
            if (TO_UPPER(((PDBSeqIdPtr)a->data.ptrvalue)->chain) !=
                TO_UPPER(((PDBSeqIdPtr)b->data.ptrvalue)->chain))
                return SIC_NO;
			*/
            if (((PDBSeqIdPtr)a->data.ptrvalue)->chain !=
                ((PDBSeqIdPtr)b->data.ptrvalue)->chain)
                return SIC_NO;
			return SIC_YES;
		case SEQID_GENERAL:  /* general */
			if (DbtagMatch((DbtagPtr)a->data.ptrvalue,
				(DbtagPtr)b->data.ptrvalue))
				return SIC_YES;
			else if (StringICmp(((DbtagPtr)a->data.ptrvalue)->db,
				((DbtagPtr)b->data.ptrvalue)->db))
				return SIC_DIFF; /* db strings do not match, okay */
			else
				return SIC_NO;

        case SEQID_GENBANK:
		case SEQID_EMBL:
		case SEQID_DDBJ:
		case SEQID_PIR:
		case SEQID_SWISSPROT:
		case SEQID_PRF:
		case SEQID_OTHER:
		case SEQID_TPG:
		case SEQID_TPE:
		case SEQID_TPD:
        case SEQID_GPIPE:

            at = (TextSeqIdPtr)a->data.ptrvalue;
            bt = (TextSeqIdPtr)b->data.ptrvalue;
            if ((at->accession != NULL) && (bt->accession != NULL))
            {
                if (! StringICmp(at->accession, bt->accession)) {
                    if (at->version > 0 &&
                        bt->version > 0 &&
                        at->version != bt->version) {
                        return SIC_NO;
                    }
                    return SIC_YES;
                } else {
                    return SIC_NO;
                }
            }
            else if ((at->name != NULL) && (bt->name != NULL))
            {
                if (! StringICmp(at->name, bt->name)) {
                    if (at->version > 0 &&
                        bt->version > 0 &&
                        at->version != bt->version) {
                        return SIC_NO;
                    }
                    return SIC_YES;
                } else {
                    return SIC_NO;
                }
            }
            else
                return SIC_DIFF;
		default:
			//ErrPostEx(SEV_ERROR, 0,0, "SeqIdComp: unsupported type [%d]",
			//	(int)choice);
			throw "Error in NCBI utility code";
			return SIC_DIFF;
     }
}

CharPtr SeqIdPrint(SeqIdPtr isip, CharPtr buf, Uint1 format)

{	
	return SeqIdWrite (isip, buf, format, 255); /* no knowledge of buffer size */
}

static const char * delim = "|";
static const char * txtid [NUM_SEQID] = {		  /* FASTA_LONG formats */
	"???" ,		/* not-set = ??? */
	"lcl",		/* local = lcl|integer or string */
	"bbs",      /* gibbsq = bbs|integer */
	"bbm",		/* gibbmt = bbm|integer */
	"gim",		/* giim = gim|integer */
	"gb",		/* genbank = gb|accession|locus */
	"emb",		/* embl = emb|accession|locus */
	"pir",		/* pir = pir|accession|name */
	"sp",		/* swissprot = sp|accession|name */
	"pat",		/* patent = pat|country|patent number (string)|seq number (integer) */
	"ref",		/* other = ref|accession|name|release - changed from oth to ref */
	"gnl",		/* general = gnl|database(string)|id (string or number) */
	"gi",		/* gi = gi|integer */
	"dbj",		/* ddbj = dbj|accession|locus */
	"prf",		/* prf = prf|accession|name */
	"pdb",		/* pdb = pdb|entry name (string)|chain id (char) */
	"tpg",      /* tpg = tpg|accession|name */
	"tpe",      /* tpe = tpe|accession|name */
	"tpd",      /* tpd = tpd|accession|name */
	"gpp"};     /* gpp = gpp|accession|name */

CharPtr SeqIdWrite(SeqIdPtr isip, CharPtr buf, Uint1 format, Uint4 buflen)

{
	SeqIdPtr sip;
	char localbuf[32];    /* for MS Windows */
	char *ldelim;
	char d [2];
	CharPtr tmp;
	
	static Uint1 fasta_order[NUM_SEQID] = {  /* order for other id FASTA_LONG */
		33, /* 0 = not set */
		20, /* 1 = local Object-id */
		15,  /* 2 = gibbsq */
		16,  /* 3 = gibbmt */
		30, /* 4 = giim Giimport-id */
		10, /* 5 = genbank */
		10, /* 6 = embl */
		10, /* 7 = pir */
		10, /* 8 = swissprot */
		15,  /* 9 = patent */
		12, /* 10 = other TextSeqId */
		13, /* 11 = general Dbtag */
		255,  /* 12 = gi */
		10, /* 13 = ddbj */
		10, /* 14 = prf */
		12,  /* 15 = pdb */
		10,  /* 16 = tpg */
		10,  /* 17 = tpe */
		10,  /* 18 = tpd */
		15   /* 19 = gpp */
    	};
	
	static Uint1 tmsmart_order[NUM_SEQID] = {  /* order for other id FASTA_LONG */
		33, /* 0 = not set */
		20, /* 1 = local Object-id */
		15,  /* 2 = gibbsq */
		16,  /* 3 = gibbmt */
		30, /* 4 = giim Giimport-id */
		10, /* 5 = genbank */
		10, /* 6 = embl */
		10, /* 7 = pir */
		10, /* 8 = swissprot */
		15,  /* 9 = patent */
		12, /* 10 = other TextSeqId */
		29, /* 11 = general Dbtag */
		255,  /* 12 = gi */
		10, /* 13 = ddbj */
		10, /* 14 = prf */
		12,  /* 15 = pdb */
		10,  /* 16 = tpg */
		10,  /* 17 = tpe */
		10,  /* 18 = tpd */
		15   /* 19 = gpp */
	};
	
	static Uint1 general_order[NUM_SEQID] = {  /* order for other id FASTA_LONG */
		33, /* 0 = not set */
		20, /* 1 = local Object-id */
		15,  /* 2 = gibbsq */
		16,  /* 3 = gibbmt */
		30, /* 4 = giim Giimport-id */
		10, /* 5 = genbank */
		10, /* 6 = embl */
		10, /* 7 = pir */
		10, /* 8 = swissprot */
		15,  /* 9 = patent */
		13, /* 10 = other TextSeqId */
		12, /* 11 = general Dbtag */
		255,  /* 12 = gi */
		10, /* 13 = ddbj */
		10, /* 14 = prf */
		12,  /* 15 = pdb */
		10,  /* 16 = tpg */
		10,  /* 17 = tpe */
		10,  /* 18 = tpd */
		15   /* 19 = gpp */
	};
	
	Boolean useGeneral = FALSE;
	TextSeqIdPtr tsip;
	PDBSeqIdPtr psip;
	ObjectIdPtr oip;
	PatentSeqIdPtr patsip;
	Boolean got_gi = FALSE;
	Boolean got_tmsmart = FALSE;
	DbtagPtr dbt;
	Char chainbuf[3];
	Char versionbuf[10];
	Int2 version = 0;

	buf[0] = '\0';
	buflen--;
	tmp = buf;
	if (isip == NULL)
		return tmp;

	d [0] = *delim;
	d [1] = '\0';
	ldelim = &(d [0]);
	
	if ((format >= ' ') && (format <= 127))  /* change delimiter */
	{
		if (format == 127)
			d [0] = '\t';
		else
			d [0] = (char) format;
		format = PRINTID_FASTA_SHORT;
	}

	if (format == PRINTID_FASTA_GENERAL) {
		useGeneral = TRUE;
		format = PRINTID_FASTA_LONG;
	}

	if (format == PRINTID_FASTA_ALL) {
	
		Char allbuf [41];
		ValNodePtr vnp, head = NULL;
		size_t len = 0;
		CharPtr str;
		Boolean notfirst;

		for (sip = isip; sip != NULL; sip = sip->next) {
			SeqIdWrite (sip, allbuf, PRINTID_FASTA_SHORT, sizeof (allbuf) - 1);
			ValNodeCopyStr (&head, 0, allbuf);
		}
		for (vnp = head; vnp != NULL; vnp = vnp->next) {
		  str = (CharPtr) vnp->data.ptrvalue;
		  if (! StringHasNoText (str)) {
		    len += StringLen (str) + 1;
		  }
		}
		if (len < 1) return buf;
		tmp = MemNew (len + 2);
		if (tmp == NULL) return buf;
		notfirst = FALSE;
		for (vnp = head; vnp != NULL; vnp = vnp->next) {
		  str = (CharPtr) vnp->data.ptrvalue;
		  if (! StringHasNoText (str)) {
		    if (notfirst) {
		      StringCat (tmp, "|");
		    }
		    StringCat (tmp, str);
		    notfirst = TRUE;
		  }
		}
		ValNodeFreeData (head);
		StringNCpy_0 (buf, tmp, buflen + 1);
		MemFree (tmp);
		return buf;
	}
	
	localbuf[0] = '\0';
							/* error on input, return ??? */
	if ( (! (isip -> choice)) || (format < PRINTID_FASTA_SHORT)
		|| (format > PRINTID_REPORT))
	{
		Nlm_LabelCopyNext(&tmp, txtid[0], &buflen);
		return tmp;
	}

	if (format == PRINTID_FASTA_LONG)   /* find the ids in the chain */
	{
		for (sip = isip; sip != NULL; sip = sip->next)  /* GI present? */
		{
			if (sip->choice == SEQID_GI)
			{
				sprintf(localbuf, "%s%s%ld", txtid[SEQID_GI], ldelim,
					(long)(sip->data.intvalue));
				Nlm_LabelCopyNext(&tmp, localbuf, &buflen);
				got_gi = TRUE;
			} else if (sip->choice == SEQID_GENERAL) {
				dbt = (DbtagPtr) sip->data.ptrvalue;
				if (dbt != NULL && StringICmp (dbt->db, "TMSMART") == 0) {
					got_tmsmart = TRUE;
				}
			}
		}
		if (useGeneral) {
			sip = SeqIdSelect(isip, general_order, NUM_SEQID);
		} else if (got_tmsmart) {
			sip = SeqIdSelect(isip, tmsmart_order, NUM_SEQID);
		} else {
			sip = SeqIdSelect(isip, fasta_order, NUM_SEQID);
		}
		if (sip == NULL)   /* only GI */
			return tmp;
		else if (got_gi)
		{
			if (sip->choice == SEQID_GIIM)   /* don't show GIIM with GI */
				return tmp;

			Nlm_LabelCopyNext(&tmp, ldelim, &buflen);
		}
		format = PRINTID_FASTA_SHORT; /* put on second (or only) SeqId in this format */
	}
	else
		sip = isip;          /* only one id processed */

							 /* deal with LOCUS and ACCESSION */
	if ((format == PRINTID_TEXTID_ACCESSION) || (format == PRINTID_TEXTID_LOCUS) ||
		(format == PRINTID_TEXTID_ACC_VER) || (format == PRINTID_TEXTID_ACC_ONLY))
	{
		if (format == PRINTID_TEXTID_ACCESSION) {
			format = PRINTID_TEXTID_ACC_ONLY;     /* current default */
		}
		switch (sip->choice)   /* get the real TextSeqId types */
		{
	        case SEQID_GENBANK:
    	    case SEQID_EMBL:
    	    case SEQID_DDBJ:
			case SEQID_PIR:
			case SEQID_SWISSPROT:
			case SEQID_PRF:
			case SEQID_OTHER:
			case SEQID_TPG:
			case SEQID_TPE:
			case SEQID_TPD:
		    case SEQID_GPIPE:
    	        tsip = (TextSeqIdPtr)sip->data.ptrvalue;
				if ((format == PRINTID_TEXTID_LOCUS) && (tsip->name != NULL)) {
					Nlm_LabelCopyNext(&tmp, tsip->name, &buflen);
					return tmp;
                } else if ((format == PRINTID_TEXTID_ACC_ONLY || format == PRINTID_TEXTID_LOCUS) 
					&& (tsip->accession != NULL)) {
					Nlm_LabelCopyNext(&tmp, tsip->accession, &buflen);
					return tmp;
                } else if ((format == PRINTID_TEXTID_ACC_VER) 
					&& (tsip->accession != NULL)) {
					if (tsip->version > 0 && tsip->release == NULL) {
						sprintf(localbuf, "%s.%d", tsip->accession,
							(int)(tsip->version));
					} else {
						sprintf(localbuf, "%s", tsip->accession);
					}
					Nlm_LabelCopyNext(&tmp, localbuf, &buflen);
					return tmp;
                }
				break;
			default:
				break;
		}
	}

	if (format == PRINTID_FASTA_SHORT)
	{
		Nlm_LabelCopyNext(&tmp, txtid[sip->choice], &buflen);
		Nlm_LabelCopyNext(&tmp, ldelim, &buflen);
	}

    switch (sip->choice) 
    {
        case SEQID_LOCAL:           /* object id */
            if ((((ObjectIdPtr)sip->data.ptrvalue)->str) == NULL)
			{
                sprintf(localbuf, "%ld", 
							(long)((ObjectIdPtr)sip->data.ptrvalue)->id);
				Nlm_LabelCopyNext(&tmp, localbuf, &buflen);
			}
            else
				Nlm_LabelCopyNext(&tmp, 
						((ObjectIdPtr)sip->data.ptrvalue)->str, &buflen);
            break;
        case SEQID_GIBBSQ:         
        case SEQID_GIBBMT:
		case SEQID_GI:
			sprintf(localbuf, "%ld", (long)sip->data.intvalue);
			Nlm_LabelCopyNext(&tmp, localbuf, &buflen);
            break;
        case SEQID_GIIM:
            sprintf(localbuf, "%ld", (long)((GiimPtr)sip->data.ptrvalue)->id);
			Nlm_LabelCopyNext(&tmp, localbuf, &buflen);
            break;
        case SEQID_GENBANK:
        case SEQID_EMBL:
        case SEQID_DDBJ:
        case SEQID_OTHER:
        case SEQID_TPG:
        case SEQID_TPE:
        case SEQID_TPD:
        case SEQID_GPIPE:
           tsip = (TextSeqIdPtr)(sip->data.ptrvalue);
	   if (((tsip->version > 0) && (tsip->release == NULL)) && SHOWVERSION)
		version = tsip->version;  /* show versions */
	   sprintf(versionbuf, ".%d", (int)version);
        case SEQID_PIR:
        case SEQID_SWISSPROT:
        case SEQID_PRF:
            tsip = (TextSeqIdPtr)sip->data.ptrvalue;
            if (tsip->accession != NULL)
			{
			   Nlm_LabelCopyNext(&tmp, tsip->accession, &buflen);
                           if (version)
				Nlm_LabelCopyNext(&tmp, versionbuf,&buflen);
			   if (format != PRINTID_FASTA_SHORT)
			 	break;
			}
			if (format == PRINTID_FASTA_SHORT)
				Nlm_LabelCopyNext(&tmp, ldelim, &buflen);
			if (tsip->name != NULL)
				Nlm_LabelCopyNext(&tmp, tsip->name, &buflen);
			/*
			if (sip->choice == SEQID_OTHER) {
				Nlm_LabelCopyNext(&tmp, ldelim, &buflen);
				if (tsip->release != NULL)
					Nlm_LabelCopyNext(&tmp, tsip->release, &buflen);
			}
			*/
            break;
        case SEQID_PATENT:
			patsip = (PatentSeqIdPtr)(sip->data.ptrvalue);
			Nlm_LabelCopyNext(&tmp, patsip->cit->country, &buflen);
			if (format == PRINTID_FASTA_SHORT)
				Nlm_LabelCopyNext(&tmp, ldelim, &buflen);
			Nlm_LabelCopyNext(&tmp, patsip->cit->number, &buflen);
			if (format == PRINTID_FASTA_SHORT)
				Nlm_LabelCopyNext(&tmp, ldelim, &buflen);
			else
				Nlm_LabelCopyNext(&tmp, "_", &buflen);
			sprintf(localbuf, "%d", (int)patsip->seqid);
			Nlm_LabelCopyNext(&tmp, localbuf, &buflen);
            break;
        case SEQID_GENERAL:
			oip = ((DbtagPtr)sip->data.ptrvalue)->tag;
			if((format == PRINTID_FASTA_SHORT) || (format == PRINTID_REPORT))
				Nlm_LabelCopyNext(&tmp, 
					((DbtagPtr)sip->data.ptrvalue)->db, &buflen);
			if (format == PRINTID_FASTA_SHORT)
				Nlm_LabelCopyNext(&tmp, ldelim, &buflen);
			else if (format == PRINTID_REPORT)
				Nlm_LabelCopyNext(&tmp, ":", &buflen);

			if (oip->str == NULL)
			{
				sprintf(localbuf, "%ld", (long) oip->id);
				Nlm_LabelCopyNext(&tmp, localbuf, &buflen);
			}
			else
				Nlm_LabelCopyNext(&tmp, oip->str, &buflen);
            break;
        case SEQID_PDB:
			psip = (PDBSeqIdPtr) sip->data.ptrvalue;
			chainbuf[0] = TO_UPPER (psip->chain);
			chainbuf[1] = '\0';
			chainbuf[2] = '\0';
			if (IS_LOWER (psip->chain)) {
			  chainbuf[1] = chainbuf [0];
			}
			Nlm_LabelCopyNext(&tmp, psip->mol, &buflen);
			if (format == PRINTID_FASTA_SHORT)
			{
				Nlm_LabelCopyNext(&tmp, ldelim, &buflen);
				if (chainbuf[0] == '|') /* special */
					Nlm_LabelCopyNext(&tmp, "VB",&buflen);
				else if (chainbuf[0] != '\0')
					Nlm_LabelCopyNext(&tmp,chainbuf, &buflen);
				else
					Nlm_LabelCopyNext(&tmp, " ", &buflen);
			}
			else if (psip->chain > ' ')
			{
				Nlm_LabelCopyNext(&tmp, "_", &buflen);
				Nlm_LabelCopyNext(&tmp,chainbuf, &buflen);
			}
            break;
		default:
			Nlm_LabelCopyNext(&tmp, txtid[0], &buflen);
			break;

    }
    return tmp;
}

#define SEQID_PARSE_BUF_SIZE 200
SeqIdPtr SeqIdParse(CharPtr buf)
{
	char localbuf[SEQID_PARSE_BUF_SIZE + 2];
	char * tmp, *strt, * tokens[6], *chain;
	char d;
	long num;
	CharPtr tp;
	Int2 numtoken, i, type = 0, j, ctr=0, numdigits; /* ctr is number of OK ids done */
	SeqIdPtr sip = NULL, head = NULL, last = NULL, tmpsip;
	ObjectIdPtr oip;
	DbtagPtr dp;
	TextSeqIdPtr tsip;
	PatentSeqIdPtr patsip;
	IdPatPtr ipp;
	PDBSeqIdPtr psip;
	GiimPtr gim;
	Boolean done = FALSE;
	static Uint1 expect_tokens[NUM_SEQID] = {  /* number of tokens to expect */
 	0, /* 0 = not set */
	1, /* 1 = local Object-id */
	1,  /* 2 = gibbsq */
	1,  /* 3 = gibbmt */
	1, /* 4 = giim Giimport-id */
	2, /* 5 = genbank */
	2, /* 6 = embl */
	2, /* 7 = pir */
	2, /* 8 = swissprot */
	3,  /* 9 = patent */
	3, /* 10 = other TextSeqId */
	2, /* 11 = general Dbtag */
	1,  /* 12 = gi */
	2, /* 13 = ddbj */
	2, /* 14 = prf */
	2,  /* 15 = pdb */
    2,  /* 16 = tpg */
    2,  /* 17 = tpe */
    2,  /* 18 = tpd */
    2   /* 19 = gpp */
	};

	if ((buf == NULL) || (*buf == '\0'))
		return NULL;

	d = *delim;   /* delimiter */
	while (! done)
	{
						/* set all tokens pointing to \0 */
		localbuf[SEQID_PARSE_BUF_SIZE + 1] = '\0';
		for (i = 0; i < 6; i++)
			tokens[i] = &localbuf[SEQID_PARSE_BUF_SIZE + 1];
		tp = buf;		/* save start of string */
						/* copy and tokenize - token\0token\0\n */
		for (tmp=localbuf, i=0; ((*buf != d) && (*buf != '\0') && (i < SEQID_PARSE_BUF_SIZE));
				i++,buf++,tmp++)
			*tmp = *buf;
		if (*buf != d) goto erret;  /* didn't get delimiter */
		*tmp = '\0';
		tmp++;
		buf++;
		for (j = 0, type = 0; j < NUM_SEQID; j++)
		{
			if (! StringCmp(localbuf, txtid[j]))
			{
				type = j;
				break;
			}
		}

		/* oth now ref, but still want to parse old style */
		if ((! type) && (! StringCmp(localbuf, "oth"))) {
			type = SEQID_OTHER;
		}

		if (! type) goto erret;

						/* copy and tokenize - token\0token\0\n */
		for (numtoken=0, strt=tmp;
			((i < SEQID_PARSE_BUF_SIZE) && (numtoken < (Int2)(expect_tokens[type])) && (! done));
			i++,buf++,tmp++)
		{
			if ((*buf == d) || (*buf == '\0'))
			{
				*tmp = '\0';
				tokens[numtoken] = strt;
				numtoken++;
				if (*buf == '\0')
				{
					if (type == SEQID_OTHER && (numtoken == 2 || numtoken == 1))
						done = TRUE;
					else if ((type == SEQID_GENBANK || type == SEQID_EMBL ||
							type == SEQID_DDBJ || type == SEQID_TPG ||
							type == SEQID_TPE || type == SEQID_TPD ||
							type == SEQID_GPIPE) && numtoken == 1)
						done = TRUE;
					else if (numtoken < (Int2)(expect_tokens[type]))
						goto erret;
					else
						done = TRUE;
				}
				strt = tmp+1;
			}
			else
				*tmp = *buf;
		}
		if (i == SEQID_PARSE_BUF_SIZE) goto erret;

		sip = ValNodeNew(head);
		if (head == NULL) head = sip;
		sip->choice = (Uint1) type;
		switch (type)
    	{
        	case SEQID_LOCAL:           /* object id */
				if (*tokens[0] == '\0') goto erret;
				oip = ObjectIdNew();
				sip->data.ptrvalue = oip;
				for (tmp = tokens[0], numdigits = 0; *tmp != '\0'; tmp++, numdigits++)
				{
					if (! IS_DIGIT(*tmp))   /* string type */
					{
						oip->str = StringSave(tokens[0]);
						break;
					}
				}
				if (oip->str == NULL)
				{
					sscanf(tokens[0], "%ld", &num);
					oip->id = (Int4)num;
					if (numdigits < 10 ||
						(numdigits == 10 && StringCmp (tokens [0], "2147483647") <= 0)) {
						sscanf(tokens[0], "%ld", &num);
						oip->id = (Int4)num;
					} else {
						oip->str = StringSave(tokens[0]);
					}
				}
				break;
	        case SEQID_GIBBSQ:         
    	    case SEQID_GIBBMT:
			case SEQID_GI:
				if (! IS_DIGIT(*tokens[0]))
					goto erret;
				sscanf(tokens[0], "%ld", &num);
				sip->data.intvalue = (Int4)num;
    	        break;
        	case SEQID_GIIM:
				if (! IS_DIGIT(*tokens[0])) goto erret;
				gim = GiimNew();
				sip->data.ptrvalue = gim;
				sscanf(tokens[0], "%ld", &num);
				gim->id = (Int4)num;
	            break;
    	    case SEQID_GENBANK:
        	case SEQID_EMBL:
	        case SEQID_PIR:
    	    case SEQID_SWISSPROT:
        	case SEQID_DDBJ:
			case SEQID_PRF:
    	    case SEQID_OTHER:
    	    case SEQID_TPG:
    	    case SEQID_TPE:
            case SEQID_TPD:
    	    case SEQID_GPIPE:
				if ((*tokens[0] == '\0') && (*tokens[1] == '\0'))
					goto erret;
	            tsip = TextSeqIdNew();
				sip->data.ptrvalue = tsip;
				if (*tokens[0] != '\0')
				{
                                        tmp = tokens[0]; /* check for version */
                                        while (*tmp != '\0')
					{
						if (*tmp == '.')
						{
                                                   if (IS_DIGIT(*(tmp+1)))
                                                   {
							*tmp = '\0';
                                                        sscanf((tmp+1),"%ld",&num);
                                                        tsip->version =(Int2)num;
						   }
						   else
							tmp++;
						}
						else
						  tmp++;
					}
					tsip->accession = StringSave(tokens[0]);
					*(tsip->accession) = TO_UPPER(*(tsip->accession));
				}
				if (*tokens[1] != '\0')
				{
					tsip->name = StringSave(tokens[1]);
					if (type != SEQID_OTHER) {
						tmp = tsip->name;
						while (*tmp != '\0')
						{
							*tmp = TO_UPPER(*tmp);
							tmp++;
						}
					}
				}
        	    break;
	        case SEQID_PATENT:
				if ((*tokens[0] == '\0') || (*tokens[1] == '\0')) goto erret;
				if (! IS_DIGIT(*tokens[2])) goto erret;
				patsip = PatentSeqIdNew();
				sip->data.ptrvalue = patsip;
				ipp = IdPatNew();
				patsip->cit = ipp;
				ipp->country = StringSave(tokens[0]);
				ipp->number = StringSave(tokens[1]);
				sscanf(tokens[2], "%ld", &num);
				patsip->seqid = (Int2)num;
            	break;
	        case SEQID_GENERAL:
				if ((*tokens[0] == '\0') || (*tokens[1] == '\0')) goto erret;
				dp = DbtagNew();
				sip->data.ptrvalue = dp;
				oip = ObjectIdNew();
				dp->tag = oip;
				dp->db = StringSave(tokens[0]);
				for (tmp = tokens[1], numdigits = 0; *tmp != '\0'; tmp++, numdigits++)
				{
					if (! IS_DIGIT(*tmp))   /* string type */
					{
						oip->str = StringSave(tokens[1]);
						break;
					}
				}
				if (oip->str == NULL)
				{
					if (numdigits < 10 ||
						(numdigits == 10 && StringCmp (tokens [1], "2147483647") <= 0)) {
						sscanf(tokens[1], "%ld", &num);
						oip->id = (Int4)num;
					} else {
						oip->str = StringSave(tokens[1]);
					}
				}
        	    break;
	        case SEQID_PDB:
				if (*tokens[0] == '\0') goto erret;
				psip = PDBSeqIdNew();
				sip->data.ptrvalue = psip;
				psip->mol = StringSave(tokens[0]);
				tmp = psip->mol;
				while (*tmp != '\0')
				{
					*tmp = TO_UPPER(*tmp);
					tmp++;
				}
				chain = tokens [1];
				if ((! StringICmp(tokens[1], "VB")) ||
                                    *(buf-1) == d)
					psip->chain = '|';
				else if (! StringHasNoText (tokens[1]))
					psip->chain = *tokens[1];
				/* double letter for chain indicates lower case */
				if (StringLen (chain) == 2 && TO_UPPER (chain [0]) == TO_UPPER (chain [1])) {
					psip->chain = TO_LOWER(psip->chain);
				} else {
					psip->chain = TO_UPPER(psip->chain);
				}
        	    break;
	    }
		last = sip;
		sip = NULL;
		ctr++;
	}
ret:
	return head;
erret:
	StringNCpy(localbuf, tp, SEQID_PARSE_BUF_SIZE);
	localbuf[SEQID_PARSE_BUF_SIZE] = '\0';
	//ErrPostEx(SEV_INFO, 0,0, "SeqIdParse Failure at %s", localbuf);
	throw "Error in NCBI utility code";
	
	if (sip == head)
		head = NULL;
	else
	{
		if (last != NULL)
			last->next = NULL;
		if (! ctr)     /* no good SeqIds */
			head = SeqIdSetFree(head);
		else	       /* at least one good SeqId.. keep it */
		{
			tmpsip = head;
			last = NULL;
			for (i = 0; i < ctr; i++)
			{
				last = tmpsip;
				tmpsip = tmpsip->next;
			}
			if (last != NULL)
				last->next = NULL;
			SeqIdSetFree(tmpsip);
		}
	}
	ValNodeFree(sip);
	goto ret;
}

SeqIdPtr SeqIdDup(SeqIdPtr oldid)
{
    TextSeqIdPtr at, bt;
	GiimPtr ga, gb;
	PatentSeqIdPtr psa, psb;
	PDBSeqIdPtr pdba, pdbb;
	SeqIdPtr newid = NULL;

    if (oldid == NULL)
        return oldid;

	newid = ValNodeNew(NULL);
	if (newid == NULL) return newid;
	MemCopy(newid, oldid, sizeof(ValNode));
	newid->next = NULL;    /* not in chain */
    switch (oldid->choice)
    {
        case SEQID_NOT_SET:
            break;
        case SEQID_LOCAL:
			newid->data.ptrvalue = ObjectIdDup((ObjectIdPtr)oldid->data.ptrvalue);
			break;
			                 /* integer types */
        case SEQID_GIBBSQ:   /* gibbsq */
        case SEQID_GIBBMT:   /* gibbmt */
        case SEQID_GI:  /* gi */
            break;

        case SEQID_GIIM:   /* giim */
			ga = (GiimPtr) oldid->data.ptrvalue;
			gb = GiimNew();
			if (gb == NULL) return NULL;
			gb->id = ga->id;
			gb->db = StringSave(ga->db);
			gb->release = StringSave(ga->release);
			newid->data.ptrvalue = gb;
			break;
        case SEQID_PATENT:   /* patent seq */
			psa = (PatentSeqIdPtr)oldid->data.ptrvalue;
			psb = PatentSeqIdNew();
			if (psb == NULL) return NULL;
			psb->seqid = psa->seqid;
			psb->cit = IdPatNew();
			psb->cit->country = StringSave(psa->cit->country);
			psb->cit->number = StringSave(psa->cit->number);
			psb->cit->app_number = StringSave(psa->cit->app_number);
			newid->data.ptrvalue = psb;
			break;
							/* TextSeqId Types */
		case SEQID_GENBANK:
		case SEQID_EMBL:
		case SEQID_PIR:
		case SEQID_SWISSPROT:
		case SEQID_OTHER:
		case SEQID_DDBJ:
		case SEQID_PRF:
		case SEQID_TPG:
		case SEQID_TPE:
		case SEQID_TPD:
        case SEQID_GPIPE:
			at = (TextSeqIdPtr)oldid->data.ptrvalue;
            bt = TextSeqIdNew();
			if (bt == NULL) return NULL;
			bt->name = StringSave(at->name);
			bt->accession = StringSave(at->accession);
			bt->release = StringSave(bt->release);
			bt->version = at->version;
			newid->data.ptrvalue = bt;
			break;
		case SEQID_GENERAL:
			newid->data.ptrvalue = DbtagDup((DbtagPtr)oldid->data.ptrvalue);
			break;
		case SEQID_PDB:
			pdba = (PDBSeqIdPtr)oldid->data.ptrvalue;
            pdbb = PDBSeqIdNew();
			if (pdbb == NULL) return NULL;
			newid->data.ptrvalue = pdbb;
			pdbb->mol = StringSave(pdba->mol);
			pdbb->chain = pdba->chain;
			pdbb->rel = DateDup(pdba->rel);
			break;
     }
	return newid;
}

ValNodePtr ValNodeNew(ValNodePtr vnp)
{
	ValNodePtr newnode;

	newnode = (ValNodePtr) Nlm_MemNew(sizeof(ValNode));
	if (vnp != NULL)
	{
		while (vnp->next != NULL)
			vnp = vnp->next;
			vnp->next = newnode;
		}
		return newnode;
}

SeqIdPtr SeqIdFree(SeqIdPtr anp)
{
    Pointer pnt;

    if (anp == NULL)
        return anp;
    
    pnt = anp->data.ptrvalue;
    switch (anp->choice)
    {
        case SEQID_LOCAL:      /* local */
            ObjectIdFree((ObjectIdPtr)pnt);
            break;
        case SEQID_GIBBSQ:      /* gibbseq */
        case SEQID_GIBBMT:      /* gibbmt */
            break;
        case SEQID_GIIM:      /* giimid */
            GiimFree((GiimPtr)pnt);
            break;
        case SEQID_GENBANK:      /* genbank */
        case SEQID_EMBL:      /* embl */
        case SEQID_PIR:      /* pir   */
        case SEQID_SWISSPROT:      /* swissprot */
        case SEQID_OTHER:     /* other */
        case SEQID_DDBJ:
		case SEQID_PRF:
    	case SEQID_TPG:
	    case SEQID_TPE:
	    case SEQID_TPD:
        case SEQID_GPIPE:
            TextSeqIdFree((TextSeqIdPtr)pnt);
            break;
        case SEQID_PATENT:      /* patent seq id */
            PatentSeqIdFree((PatentSeqIdPtr)pnt);
            break;
        case SEQID_GENERAL:     /* general */
            DbtagFree((DbtagPtr)pnt);
            break;
        case SEQID_GI:     /* gi */
            break;
		case SEQID_PDB:
			PDBSeqIdFree((PDBSeqIdPtr)pnt);
			break;
    }

	//ObjMgrDelete(OBJ_SEQID, (Pointer)anp);

	return (SeqIdPtr)MemFree(anp);
}

Boolean ObjectIdMatch(ObjectIdPtr a, ObjectIdPtr b)
{
	if (a == b)
		return TRUE;

    if ((a == NULL) || (b == NULL))   /* only one is null */
        return FALSE;

	if ((a->str != NULL) && (b->str != NULL))  /* same type */
	{
	    if (StringICmp(a->str, b->str))
    	    return FALSE;
		else
			return TRUE;
	}
    else if ((a->str == NULL) && (b->str == NULL))  /* must be same kind */
    {
        if (a->id == b->id)
            return TRUE;
        else
            return FALSE;
    }
    else                   /* different kinds */
        return FALSE;
}

Boolean IdPatMatch(IdPatPtr a, IdPatPtr b)
{
    if ((a == NULL) || (b == NULL))
        return FALSE;

    if (StringICmp(a->country, b->country))   /* countries must match */
        return FALSE;

    if ((a->number != NULL) && (b->number != NULL))
    {
        if (! StringICmp(a->number, b->number))
            return TRUE;
        else
            return FALSE;
    }
    else
    {
        if (! StringICmp(a->app_number, b->app_number))
            return TRUE;
        else
            return FALSE;
    }
}

int Nlm_StringICmp(const char FAR *a, const char FAR *b)
{
	return (a && b) ? Nlm_StrICmp(a, b) : (a ? 1 : (b ? -1 : 0));
}

int Nlm_StrICmp(const char FAR *a, const char FAR *b)
{
	int diff, done;

	if (a == b)   return 0;

	done = 0;
	while (! done)
	{
		diff = TO_UPPER(*a) - TO_UPPER(*b);
		if (diff)
			return (Nlm_Int2) diff;
		if (*a == '\0')
			done = 1;
		else
		{
			a++; b++;
		}
	}
	return 0;
}

Boolean DbtagMatch(DbtagPtr a, DbtagPtr b)
{
	if (a == b)
		return TRUE;

	if ((a == NULL) || (b == NULL))
		return FALSE;

	if (StringICmp(a->db, b->db))
		return FALSE;

	return ObjectIdMatch(a->tag, b->tag);
}

ValNodePtr ValNodeCopyStr(ValNodePtr PNTR head, Nlm_Int2 choice, Nlm_CharPtr str)
{
	ValNodePtr newnode;

	if (str == NULL) return NULL;

	newnode = ValNodeAdd(head);
	if (newnode != NULL)
	{
		newnode->choice = (Nlm_Uint1)choice;
		newnode->data.ptrvalue = StringSave(str);
	}

	return newnode;
}

ValNodePtr ValNodeAdd(ValNodePtr PNTR head)
{
	ValNodePtr newnode;

	if (head != NULL)
	{
		newnode = ValNodeNew(*head);
		if (*head == NULL)
			*head = newnode;
	}
	else
		newnode = ValNodeNew(NULL);

	return newnode;
}

Nlm_Boolean Nlm_StringHasNoText (const char FAR *str)
{
	if (str) {
		while (*str) {
			if ((unsigned char)(*str++) > ' ')
				return FALSE;
		}
	}
	return TRUE;
}

size_t Nlm_StringLen(const char *str)
{
	return str ? StrLen (str) : 0;
}

char* Nlm_MemNew(size_t size)
{
	if(size == 0){
		return NULL;
	}
	
	char* ptr = new char [size];
	
	if(!ptr){
		throw "Unable to allocate memory";
	}
	
	// Initialize this memory to zero
	memset(ptr, 0, size);
	
	return ptr;
}

Nlm_CharPtr Nlm_StringCat(char FAR *to, const char FAR *from)
{
	return (to && from) ? strcat(to, from) : to;
}

ValNodePtr ValNodeFreeData(ValNodePtr vnp)
{
	ValNodePtr next;

	while (vnp != NULL)
	{
		Nlm_MemFree(vnp->data.ptrvalue);
		next = vnp->next;
		Nlm_MemFree(vnp);
		vnp = next;
	}
	return NULL;
}

void* Nlm_MemFree(void *ptr)
{
	if(ptr == NULL){
		return NULL;
	}
	
	delete [] (char*)ptr;
	
	return NULL;
}

Nlm_CharPtr Nlm_StringNCpy_0 (char FAR *to, const char FAR *from, size_t max)
{
	if (to != NULL  &&  max > 0)
		to[0] = '\0';

	if (from != NULL)
		StrNCat(to, from, max - 1);

	return to;
}

Nlm_Uint4 Nlm_LabelCopy (Nlm_CharPtr to, const char* from, Nlm_Uint4 buflen)
{
	Nlm_Uint4 len;

	if( (to == NULL) || (from == NULL) ){
		return 0;
	}

	if (buflen == 0)         /* this is a sign of multiple writes */
	{
		*(to-1) = '>';
		return 0;
	}
	
	len = buflen;

	while ((*from != '\0') && (buflen))
	{
		*to = *from;
		from++; to++; buflen--;
	}

	if (*from != '\0')
	{
		*(to - 1) = '>';
	}

	*to = '\0';      /* buffer is bufferlen+1 */
	return (Nlm_Uint4)(len - buflen);
}

void Nlm_LabelCopyNext(Nlm_CharPtr PNTR to, const char* from, Nlm_Uint4 PNTR buflen)
{
	Nlm_Uint4 diff;

	diff = Nlm_LabelCopy(*to, from, *buflen);
	*buflen -= diff; *to += diff;
	
}

SeqIdPtr SeqIdSelect(SeqIdPtr sip, Uint1Ptr order, Int2 num)
{
    SeqIdPtr bestid;

	if ((sip == NULL) || (order == NULL))
		return NULL;

    for ( bestid = NULL; sip != NULL; sip = sip -> next) 
    {
		if ((Int2)sip->choice < num)
		{
			if (order[sip->choice] < 255)
			{
				if (bestid == NULL)
					bestid = sip;
				else if (order[sip->choice] < order[bestid->choice])
					bestid = sip;
			}
		} else {
			//ErrPostEx(SEV_ERROR, 0,0, "SeqIdSelect: choice [%d] out of range [%d]",
			//	(int)(sip->choice), (int)num);
			throw "SeqIdSelect: choice out of range";
			if(sip->choice > NUM_SEQID) /*** something is really wrong ***/
				return NULL;
		}
    }

    return bestid;
}

int Nlm_StringCmp(const char FAR *a, const char FAR *b)
{
	return (a && b) ? strcmp(a, b) : (a ? 1 : (b ? -1 : 0));
}

ObjectIdPtr ObjectIdNew (void)
{
	ObjectIdPtr oid;

	oid = (ObjectIdPtr)MemNew(sizeof(ObjectId));
	return oid;
}

Nlm_CharPtr Nlm_StringSave (const char FAR *from)
{
	return from ? Nlm_StrSave (from) : 0;
}

Nlm_CharPtr Nlm_StrSave (const char FAR *from)
{
	size_t len;
	Nlm_CharPtr to;

	len = Nlm_StringLen(from);
	if ((to = (Nlm_CharPtr) Nlm_MemNew(len + 1)) != NULL) {
		memcpy(to, from, len +1);
	}
	return to;
}

GiimPtr GiimNew (void)
{
	return (GiimPtr)MemNew(sizeof(Giim));
}

TextSeqIdPtr TextSeqIdNew (void)
{
	TextSeqIdPtr tsip;

	tsip = (TextSeqIdPtr)MemNew(sizeof(TextSeqId));
	if (tsip == NULL) return tsip;

	tsip->version = INT2_MIN;
	return tsip;
}

PatentSeqIdPtr PatentSeqIdNew (void)
{
	return (PatentSeqIdPtr)MemNew(sizeof(PatentSeqId));
}

IdPatPtr IdPatNew (void)
{
	IdPatPtr idp;

	idp = (IdPatPtr)MemNew(sizeof(IdPat));
	return idp;
}

DbtagPtr DbtagNew (void)
{
	DbtagPtr dbt;

	dbt = (DbtagPtr)MemNew(sizeof(Dbtag));
	return dbt;
}

PDBSeqIdPtr PDBSeqIdNew (void)
{
	PDBSeqIdPtr pdbsip;

	pdbsip = (PDBSeqIdPtr)MemNew(sizeof(PDBSeqId));
	if (pdbsip == NULL) return pdbsip;
	pdbsip->chain = (Uint1)32;
	return pdbsip;
}

Nlm_CharPtr Nlm_StringNCpy (char FAR *to, const char FAR *from, size_t max)
{
	return (to && from) ? strncpy(to, from, max) : Nlm_ClearDestString (to, max);
}

Nlm_CharPtr Nlm_ClearDestString(Nlm_CharPtr to, size_t max)
{
	if (to != NULL && max > 0) {
		memset (to, 0, max);
		*to = '\0';
	}
	return to;
}

ValNodePtr ValNodeFree (ValNodePtr vnp)
{
	ValNodePtr next;

	while (vnp != NULL)
	{
		next = vnp->next;
		Nlm_MemFree(vnp);
		vnp = next;
	}
	return NULL;
}

void * Nlm_MemCopy(void *dst, const void *src, size_t bytes)
{
	return (dst&&src) ? memcpy(dst, src, bytes) : NULL;
}

ObjectIdPtr ObjectIdDup(ObjectIdPtr oldid)
{
	ObjectIdPtr newid;

	if (oldid == NULL)
		return oldid;
   
	newid = ObjectIdNew();
	newid->id = oldid->id;
	newid->str = StringSave(oldid->str);
	return newid;
}

DbtagPtr DbtagDup (DbtagPtr oldtag)
{
	DbtagPtr newtag;

	if (oldtag == NULL)
		return oldtag;

	newtag = DbtagNew();
	if (newtag == NULL) return newtag;
	newtag->db = StringSave(oldtag->db);
	newtag->tag = ObjectIdDup(oldtag->tag);
	return newtag;
}

NCBI_DatePtr DateDup(NCBI_DatePtr dp)
{
	NCBI_DatePtr np;

	if (dp == NULL)
		return dp;

	np = DateNew();
	if (np == NULL) return np;
	MemCopy(np, dp, sizeof(NCBI_Date));
	if (dp->str != NULL)
		np->str = StringSave(dp->str);
	return np;
}

ObjectIdPtr ObjectIdFree (ObjectIdPtr oid)
{
	if (oid == NULL)
		return oid;

	MemFree(oid->str);
	return (ObjectIdPtr)MemFree(oid);
}

GiimPtr GiimFree(GiimPtr gip)
{
	if (gip == NULL)
		return gip;

	MemFree(gip->db);
	MemFree(gip->release);
	return (GiimPtr)MemFree(gip);
}

TextSeqIdPtr TextSeqIdFree(TextSeqIdPtr tsip)
{
	if (tsip == NULL)
		return tsip;

	MemFree(tsip->name);
	MemFree(tsip->accession);
	MemFree(tsip->release);
	return (TextSeqIdPtr)MemFree(tsip);
}

PatentSeqIdPtr PatentSeqIdFree(PatentSeqIdPtr psip)
{
	if (psip == NULL)
		return psip;

	IdPatFree(psip->cit);
	return (PatentSeqIdPtr)MemFree(psip);
}

DbtagPtr DbtagFree(DbtagPtr dbt)
{
	if (dbt == NULL)
		return dbt;

	dbt->tag = ObjectIdFree(dbt->tag);
	MemFree(dbt->db);
	return (DbtagPtr)MemFree(dbt);
}

PDBSeqIdPtr PDBSeqIdFree(PDBSeqIdPtr pdbsip)
{
	if (pdbsip == NULL)
		return pdbsip;

	MemFree(pdbsip->mol);
	DateFree(pdbsip->rel);
	return (PDBSeqIdPtr)MemFree(pdbsip);
}

NCBI_DatePtr DateNew (void)
{
	NCBI_DatePtr dp = NULL;
       
	dp = (NCBI_DatePtr)MemNew(sizeof(NCBI_Date));
	dp->data[4] = 255;  /* hour not set */
	dp->data[5] = 255;  /* minute not set */
	dp->data[6] = 255;  /* second not set */
	return dp;
}

IdPatPtr IdPatFree(IdPatPtr idp)
{
	if (idp == NULL)
		return idp;

	MemFree(idp->country);
	MemFree(idp->number);
	MemFree(idp->app_number);
	MemFree(idp->doc_type);
	return (IdPatPtr)MemFree(idp);
}

NCBI_DatePtr DateFree (NCBI_DatePtr dp)
{
	if (dp == NULL)
		return dp;

	MemFree(dp->str);
	return (NCBI_DatePtr)MemFree(dp);
}

#endif // USE_NCBI
