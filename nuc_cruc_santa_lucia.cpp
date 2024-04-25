#include "nuc_cruc.h"
#include <math.h>

// Linear interpolation for missing entropy values
#define	INTERPOLATE(S, BEGIN, END, X) ( S[BEGIN] + (S[END] - S[BEGIN])*float(X - BEGIN)/(END - BEGIN) )
	
void NucCruc::init_param_Santa_Lucia()
{
	tm_param = SANTA_LUCIA;
		
	// Initialize all parameters to 0 by default
	memset( param_H, 0, NUM_BASE_PAIR*NUM_BASE_PAIR*sizeof(float) );
	memset( param_S, 0, NUM_BASE_PAIR*NUM_BASE_PAIR*sizeof(float) );
	
	// Initialize forbidden interactions:
	// 
	#define	default_enthalpy	100.0f
	#define	default_entropy	0.0f
	for(char i = BASE::A;i <= BASE::I;i++){
		for(char j = BASE::A;j <= BASE::I;j++){
			
			// Gap bound to gap is not allowed
			int curr = BASE_PAIR(i, j);
			int prev = BASE_PAIR(BASE::GAP, BASE::GAP);
			
			param_H[SCORE_INDEX(curr, prev)] = param_H[SCORE_INDEX(prev, curr)] = default_enthalpy;
			param_S[SCORE_INDEX(curr, prev)] = param_S[SCORE_INDEX(prev, curr)] = default_entropy;
			
			// x-/-y is not allowed
			curr = BASE_PAIR(i, BASE::GAP);
			prev = BASE_PAIR(BASE::GAP, j);
			
			param_H[SCORE_INDEX(curr, prev)] = param_H[SCORE_INDEX(prev, curr)] = default_enthalpy;
			param_S[SCORE_INDEX(curr, prev)] = param_S[SCORE_INDEX(prev, curr)] = default_entropy;
			
			//-x/y- is not allowed
			curr = BASE_PAIR(BASE::GAP, i);
			prev = BASE_PAIR(j, BASE::GAP);
			
			param_H[SCORE_INDEX(curr, prev)] = param_H[SCORE_INDEX(prev, curr)] = default_enthalpy;
			param_S[SCORE_INDEX(curr, prev)] = param_S[SCORE_INDEX(prev, curr)] = default_entropy;
		}
	}

	// Exact match parameters are from:
	// "The thermodynamics of DNA structural motifs"
	// J. SantaLucia and D. Hicks
	// Annu. Rev. Biophys. Biomol. Struc. 2004. 33:415-440
	
	// Note that this article has an incorrect value for AT_AT
	// (it does not agree with all previous published values from the 
	// same group and does not match the Tm web server provided 
	// by the SantaLucia group).
	
	// Exact nearest-neighboor matches: Delta H 
	#define	AT_AT	-7.9f
	#define	AT_CG	-8.4f
	#define	AT_GC	-7.8f
	#define	AT_TA	-7.2f
	#define	CG_AT	-8.5f
	#define	CG_CG	-8.0f
	#define	CG_GC	-10.6f
	#define	GC_AT	-8.2f
	#define	GC_CG	-9.8f
	#define	TA_AT	-7.2f
	
	param_H[SCORE_INDEX(AT, AT)] = param_H[SCORE_INDEX(TA, TA)] = AT_AT; 	// AA/TT
	param_H[SCORE_INDEX(AT, CG)] = param_H[SCORE_INDEX(GC, TA)] = AT_CG; 	// AC/TG
	param_H[SCORE_INDEX(AT, GC)] = param_H[SCORE_INDEX(CG, TA)] = AT_GC; 	// AG/TC
	param_H[SCORE_INDEX(AT, TA)] = AT_TA; 				// AT/TA
	param_H[SCORE_INDEX(CG, AT)] = param_H[SCORE_INDEX(TA, GC)] = CG_AT; 	// CA/GT
	param_H[SCORE_INDEX(CG, CG)] = param_H[SCORE_INDEX(GC, GC)] = CG_CG; 	// CC/GG
	param_H[SCORE_INDEX(CG, GC)] = CG_GC;					// CG/GC
	param_H[SCORE_INDEX(GC, AT)] = param_H[SCORE_INDEX(TA, CG)] = GC_AT;	// GA/CT
	param_H[SCORE_INDEX(GC, CG)] = GC_CG;					// GC/CG
	param_H[SCORE_INDEX(TA, AT)] = TA_AT;					// TA/AT
	
	param_S[SCORE_INDEX(AT, AT)] = param_S[SCORE_INDEX(TA, TA)] = ENTROPY(-1.00f, AT_AT);	// AA/TT
	param_S[SCORE_INDEX(AT, CG)] = param_S[SCORE_INDEX(GC, TA)] = ENTROPY(-1.44f, AT_CG);	// AC/TG
	param_S[SCORE_INDEX(AT, GC)] = param_S[SCORE_INDEX(CG, TA)] = ENTROPY(-1.28f, AT_GC);	// AG/TC
	param_S[SCORE_INDEX(AT, TA)] = ENTROPY(-0.88f, AT_TA);				// AT/TA
	param_S[SCORE_INDEX(CG, AT)] = param_S[SCORE_INDEX(TA, GC)] = ENTROPY(-1.45f, CG_AT);	// CA/GT
	param_S[SCORE_INDEX(CG, CG)] = param_S[SCORE_INDEX(GC, GC)] = ENTROPY(-1.84f, CG_CG);	// CC/GG
	param_S[SCORE_INDEX(CG, GC)] = ENTROPY(-2.17f, CG_GC);				// CG/GC
	param_S[SCORE_INDEX(GC, AT)] = param_S[SCORE_INDEX(TA, CG)] = ENTROPY(-1.30f, GC_AT);	// GA/CT
	param_S[SCORE_INDEX(GC, CG)] = ENTROPY(-2.24f, GC_CG);				// GC/CG
	param_S[SCORE_INDEX(TA, AT)] = ENTROPY(-0.58f, TA_AT);				// TA/AT
	
	#define	AE_AT	0.2f
	#define	AE_CG	-6.3f
	#define	AE_GC	-3.7f
	#define	AE_TA	-2.9f
	#define	CE_AT	0.6f
	#define	CE_CG	-4.4f
	#define	CE_GC	-4.0f
	#define	CE_TA	-4.1f
	#define	GE_AT	-1.1f
	#define	GE_CG	-5.1f
	#define	GE_GC	-3.9f
	#define	GE_TA	-4.2f
	#define	TE_AT	-6.9f
	#define	TE_CG	-4.0f
	#define	TE_GC	-4.9f
	#define	TE_TA	-0.2f
	
	param_H[SCORE_INDEX(AE, AT)] = param_H[SCORE_INDEX(TA, EA)] = AE_AT;	// AA/ET
	param_H[SCORE_INDEX(AE, CG)] = param_H[SCORE_INDEX(GC, EA)] = AE_CG;	// AC/EG
	param_H[SCORE_INDEX(AE, GC)] = param_H[SCORE_INDEX(CG, EA)] = AE_GC;	// AG/EC
	param_H[SCORE_INDEX(AE, TA)] = param_H[SCORE_INDEX(AT, EA)] = AE_TA;	// AT/EA
	param_H[SCORE_INDEX(CE, AT)] = param_H[SCORE_INDEX(TA, EC)] = CE_AT;	// CA/ET
	param_H[SCORE_INDEX(CE, CG)] = param_H[SCORE_INDEX(GC, EC)] = CE_CG;	// CC/EG
	param_H[SCORE_INDEX(CE, GC)] = param_H[SCORE_INDEX(CG, EC)] = CE_GC;	// CG/EC
	param_H[SCORE_INDEX(CE, TA)] = param_H[SCORE_INDEX(AT, EC)] = CE_TA;	// CT/EA
	param_H[SCORE_INDEX(GE, AT)] = param_H[SCORE_INDEX(TA, EG)] = GE_AT;	// GA/ET
	param_H[SCORE_INDEX(GE, CG)] = param_H[SCORE_INDEX(GC, EG)] = GE_CG;	// GC/EG
	param_H[SCORE_INDEX(GE, GC)] = param_H[SCORE_INDEX(CG, EG)] = GE_GC;	// GG/EC
	param_H[SCORE_INDEX(GE, TA)] = param_H[SCORE_INDEX(AT, EG)] = GE_TA;	// GT/EA
	param_H[SCORE_INDEX(TE, AT)] = param_H[SCORE_INDEX(TA, ET)] = TE_AT;	// TA/ET
	param_H[SCORE_INDEX(TE, CG)] = param_H[SCORE_INDEX(GC, ET)] = TE_CG;	// TC/EG
	param_H[SCORE_INDEX(TE, GC)] = param_H[SCORE_INDEX(CG, ET)] = TE_GC;	// TG/EC
	param_H[SCORE_INDEX(TE, TA)] = param_H[SCORE_INDEX(AT, ET)] = TE_TA;	// TT/EA
	
	param_S[SCORE_INDEX(AE, AT)] = param_S[SCORE_INDEX(TA, EA)] = ENTROPY(-0.51f, AE_AT);	// AA/ET
	param_S[SCORE_INDEX(AE, CG)] = param_S[SCORE_INDEX(GC, EA)] = ENTROPY(-0.96f, AE_CG);	// AC/EG
	param_S[SCORE_INDEX(AE, GC)] = param_S[SCORE_INDEX(CG, EA)] = ENTROPY(-0.58f, AE_GC);	// AG/EC
	param_S[SCORE_INDEX(AE, TA)] = param_S[SCORE_INDEX(AT, EA)] = ENTROPY(-0.5f, AE_TA);	// AT/EA
	param_S[SCORE_INDEX(CE, AT)] = param_S[SCORE_INDEX(TA, EC)] = ENTROPY(-0.42f, CE_AT);	// CA/ET
	param_S[SCORE_INDEX(CE, CG)] = param_S[SCORE_INDEX(GC, EC)] = ENTROPY(-0.52f, CE_CG);	// CC/EG
	param_S[SCORE_INDEX(CE, GC)] = param_S[SCORE_INDEX(CG, EC)] = ENTROPY(-0.34f, CE_GC);	// CG/EC
	param_S[SCORE_INDEX(CE, TA)] = param_S[SCORE_INDEX(AT, EC)] = ENTROPY(-0.02f, CE_TA);	// CT/EA
	param_S[SCORE_INDEX(GE, AT)] = param_S[SCORE_INDEX(TA, EG)] = ENTROPY(-0.62f, GE_AT);	// GA/ET
	param_S[SCORE_INDEX(GE, CG)] = param_S[SCORE_INDEX(GC, EG)] = ENTROPY(-0.72f, GE_CG);	// GC/EG
	param_S[SCORE_INDEX(GE, GC)] = param_S[SCORE_INDEX(CG, EG)] = ENTROPY(-0.56f, GE_GC);	// GG/EC
	param_S[SCORE_INDEX(GE, TA)] = param_S[SCORE_INDEX(AT, EG)] = ENTROPY(0.48f, GE_TA);	// GT/EA
	param_S[SCORE_INDEX(TE, AT)] = param_S[SCORE_INDEX(TA, ET)] = ENTROPY(-0.71f, TE_AT);	// TA/ET
	param_S[SCORE_INDEX(TE, CG)] = param_S[SCORE_INDEX(GC, ET)] = ENTROPY(-0.58f, TE_CG);	// TC/EG
	param_S[SCORE_INDEX(TE, GC)] = param_S[SCORE_INDEX(CG, ET)] = ENTROPY(-0.61f, TE_GC);	// TG/EC
	param_S[SCORE_INDEX(TE, TA)] = param_S[SCORE_INDEX(AT, ET)] = ENTROPY(-0.10f, TE_TA);	// TT/EA
	
	#define	EA_AT	-0.7f
	#define	EA_CG	-2.1f
	#define	EA_GC	-5.9f
	#define	EA_TA	-0.5f
	#define	EC_AT	4.4f
	#define	EC_CG	-0.2f
	#define	EC_GC	-2.6f
	#define	EC_TA	4.7f
	#define	EG_AT	-1.6f
	#define	EG_CG	-3.9f
	#define	EG_GC	-3.2f
	#define	EG_TA	-4.1f
	#define	ET_AT	2.9f
	#define	ET_CG	-4.4f
	#define	ET_GC	-5.2f
	#define	ET_TA	-3.8f
	
	param_H[SCORE_INDEX(EA, AT)] = param_H[SCORE_INDEX(TA, AE)] = EA_AT;	// EA/AT
	param_H[SCORE_INDEX(EA, CG)] = param_H[SCORE_INDEX(GC, AE)] = EA_CG;	// EC/AG
	param_H[SCORE_INDEX(EA, GC)] = param_H[SCORE_INDEX(CG, AE)] = EA_GC;	// EG/AC
	param_H[SCORE_INDEX(EA, TA)] = param_H[SCORE_INDEX(AT, AE)] = EA_TA;	// ET/AA
	param_H[SCORE_INDEX(EC, AT)] = param_H[SCORE_INDEX(TA, CE)] = EC_AT;	// EA/CT
	param_H[SCORE_INDEX(EC, CG)] = param_H[SCORE_INDEX(GC, CE)] = EC_CG;	// EC/CG
	param_H[SCORE_INDEX(EC, GC)] = param_H[SCORE_INDEX(CG, CE)] = EC_GC;	// EG/CC
	param_H[SCORE_INDEX(EC, TA)] = param_H[SCORE_INDEX(AT, CE)] = EC_TA;	// ET/CA
	param_H[SCORE_INDEX(EG, AT)] = param_H[SCORE_INDEX(TA, GE)] = EG_AT;	// EA/GT
	param_H[SCORE_INDEX(EG, CG)] = param_H[SCORE_INDEX(GC, GE)] = EG_CG;	// EC/GG
	param_H[SCORE_INDEX(EG, GC)] = param_H[SCORE_INDEX(CG, GE)] = EG_GC;	// EG/GC
	param_H[SCORE_INDEX(EG, TA)] = param_H[SCORE_INDEX(AT, GE)] = EG_TA;	// ET/GA
	param_H[SCORE_INDEX(ET, AT)] = param_H[SCORE_INDEX(TA, TE)] = ET_AT;	// EA/TT
	param_H[SCORE_INDEX(ET, CG)] = param_H[SCORE_INDEX(GC, TE)] = ET_CG;	// EC/TG
	param_H[SCORE_INDEX(ET, GC)] = param_H[SCORE_INDEX(CG, TE)] = ET_GC;	// EG/TC
	param_H[SCORE_INDEX(ET, TA)] = param_H[SCORE_INDEX(AT, TE)] = ET_TA;	// ET/TA
	
	param_S[SCORE_INDEX(EA, AT)] = param_S[SCORE_INDEX(TA, AE)] = ENTROPY(-0.48f, EA_AT);	// EA/AT
	param_S[SCORE_INDEX(EA, CG)] = param_S[SCORE_INDEX(GC, AE)] = ENTROPY(-0.92f, EA_CG);	// EC/AG
	param_S[SCORE_INDEX(EA, GC)] = param_S[SCORE_INDEX(CG, AE)] = ENTROPY(-0.82f, EA_GC);	// EG/AC
	param_S[SCORE_INDEX(EA, TA)] = param_S[SCORE_INDEX(AT, AE)] = ENTROPY(-0.12f, EA_TA);	// ET/AA
	param_S[SCORE_INDEX(EC, AT)] = param_S[SCORE_INDEX(TA, CE)] = ENTROPY(-0.19f, EC_AT);	// EA/CT
	param_S[SCORE_INDEX(EC, CG)] = param_S[SCORE_INDEX(GC, CE)] = ENTROPY(-0.23f, EC_CG);	// EC/CG
	param_S[SCORE_INDEX(EC, GC)] = param_S[SCORE_INDEX(CG, CE)] = ENTROPY(-0.31f, EC_GC);	// EG/CC
	param_S[SCORE_INDEX(EC, TA)] = param_S[SCORE_INDEX(AT, CE)] = ENTROPY(0.28f, EC_TA);	// ET/CA
	param_S[SCORE_INDEX(EG, AT)] = param_S[SCORE_INDEX(TA, GE)] = ENTROPY(-0.50f, EG_AT);	// EA/GT
	param_S[SCORE_INDEX(EG, CG)] = param_S[SCORE_INDEX(GC, GE)] = ENTROPY(-0.44f, EG_CG);	// EC/GG
	param_S[SCORE_INDEX(EG, GC)] = param_S[SCORE_INDEX(CG, GE)] = ENTROPY(-0.01f, EG_GC);	// EG/GC
	param_S[SCORE_INDEX(EG, TA)] = param_S[SCORE_INDEX(AT, GE)] = ENTROPY(-0.01f, EG_TA);	// ET/GA
	param_S[SCORE_INDEX(ET, AT)] = param_S[SCORE_INDEX(TA, TE)] = ENTROPY(-0.29f, ET_AT);	// EA/TT
	param_S[SCORE_INDEX(ET, CG)] = param_S[SCORE_INDEX(GC, TE)] = ENTROPY(-0.35f, ET_CG);	// EC/TG
	param_S[SCORE_INDEX(ET, GC)] = param_S[SCORE_INDEX(CG, TE)] = ENTROPY(-0.52f, ET_GC);	// EG/TC
	param_S[SCORE_INDEX(ET, TA)] = param_S[SCORE_INDEX(AT, TE)] = ENTROPY(0.13f, ET_TA);	// ET/TA
	
	#define	AT_AG	-0.6f
	#define	AT_GA	-0.7f
	#define	CG_AG	-0.7f
	#define	CG_GA	-4.0f
	#define	GC_AG	-0.6f
	#define	GC_GA	0.5f
	#define	TA_AG	0.7f
	#define	TA_GA	3.0f
	
	// Single G-A mismatch
	// Allawi et. al. Biochemistry 1998, 37, 2170-2179
	param_H[SCORE_INDEX(AT, AG)] = param_H[SCORE_INDEX(GA, TA)] = AT_AG;	// Aa/Tg
	param_H[SCORE_INDEX(AT, GA)] = param_H[SCORE_INDEX(AG, TA)] = AT_GA;	// Ag/Ta
	param_H[SCORE_INDEX(CG, AG)] = param_H[SCORE_INDEX(GA, GC)] = CG_AG;	// Ca/Gg
	param_H[SCORE_INDEX(CG, GA)] = param_H[SCORE_INDEX(AG, GC)] = CG_GA;	// Cg/Ga
	param_H[SCORE_INDEX(GC, AG)] = param_H[SCORE_INDEX(GA, CG)] = GC_AG;	// Ga/Cg
	param_H[SCORE_INDEX(GC, GA)] = param_H[SCORE_INDEX(AG, CG)] = GC_GA;	// Gg/Ca
	param_H[SCORE_INDEX(TA, AG)] = param_H[SCORE_INDEX(GA, AT)] = TA_AG;	// Ta/Ag
	param_H[SCORE_INDEX(TA, GA)] = param_H[SCORE_INDEX(AG, AT)] = TA_GA;	// Tg/Aa
	
	param_S[SCORE_INDEX(AT, AG)] = param_S[SCORE_INDEX(GA, TA)] = ENTROPY(0.14f, AT_AG);	// Aa/Tg
	param_S[SCORE_INDEX(AT, GA)] = param_S[SCORE_INDEX(AG, TA)] = ENTROPY(0.02f, AT_GA);	// Ag/Ta
	param_S[SCORE_INDEX(CG, AG)] = param_S[SCORE_INDEX(GA, GC)] = ENTROPY(0.03f, CG_AG);	// Ca/Gg
	param_S[SCORE_INDEX(CG, GA)] = param_S[SCORE_INDEX(AG, GC)] = ENTROPY(0.11f, CG_GA);	// Cg/Ga
	param_S[SCORE_INDEX(GC, AG)] = param_S[SCORE_INDEX(GA, CG)] = ENTROPY(-0.25f, GC_AG);	// Ga/Cg
	param_S[SCORE_INDEX(GC, GA)] = param_S[SCORE_INDEX(AG, CG)] = ENTROPY(-0.52f, GC_GA);	// Gg/Ca
	param_S[SCORE_INDEX(TA, AG)] = param_S[SCORE_INDEX(GA, AT)] = ENTROPY(0.42f, TA_AG);	// Ta/Ag
	param_S[SCORE_INDEX(TA, GA)] = param_S[SCORE_INDEX(AG, AT)] = ENTROPY(0.74f, TA_GA);	// Tg/Aa
	
	#define	AT_CT	0.7f
	#define	AT_TC	-1.2f
	#define	CG_CT	-0.8f
	#define	CG_TC	-1.5f
	#define	GC_CT	2.3f
	#define	GC_TC	5.2f
	#define	TA_CT	1.2f
	#define	TA_TC	1.0f
	
	// Single C-T mismatch
	// Allawi et. al. Nucleic Acids Research, 1998, 26, 11, 2694-2701
	param_H[SCORE_INDEX(AT, CT)] = param_H[SCORE_INDEX(TC, TA)] = AT_CT;	// Ac/Tt
	param_H[SCORE_INDEX(AT, TC)] = param_H[SCORE_INDEX(CT, TA)] = AT_TC;	// At/Tc
	param_H[SCORE_INDEX(CG, CT)] = param_H[SCORE_INDEX(TC, GC)] = CG_CT;	// Cc/Gt
	param_H[SCORE_INDEX(CG, TC)] = param_H[SCORE_INDEX(CT, GC)] = CG_TC;	// Ct/Gc
	param_H[SCORE_INDEX(GC, CT)] = param_H[SCORE_INDEX(TC, CG)] = GC_CT;	// Gc/Ct
	param_H[SCORE_INDEX(GC, TC)] = param_H[SCORE_INDEX(CT, CG)] = GC_TC;	// Gt/Cc
	param_H[SCORE_INDEX(TA, CT)] = param_H[SCORE_INDEX(TC, AT)] = TA_CT;	// Tc/At
	param_H[SCORE_INDEX(TA, TC)] = param_H[SCORE_INDEX(CT, AT)] = TA_TC;	// Tt/Ac
	
	param_S[SCORE_INDEX(AT, CT)] = param_S[SCORE_INDEX(TC, TA)] = ENTROPY(0.64f, AT_CT);	// Ac/Tt
	param_S[SCORE_INDEX(AT, TC)] = param_S[SCORE_INDEX(CT, TA)] = ENTROPY(0.73f, AT_TC);	// At/Tc
	param_S[SCORE_INDEX(CG, CT)] = param_S[SCORE_INDEX(TC, GC)] = ENTROPY(0.62f, CG_CT);	// Cc/Gt
	param_S[SCORE_INDEX(CG, TC)] = param_S[SCORE_INDEX(CT, GC)] = ENTROPY(0.40f, CG_TC);	// Ct/Gc
	param_S[SCORE_INDEX(GC, CT)] = param_S[SCORE_INDEX(TC, CG)] = ENTROPY(0.62f, GC_CT);	// Gc/Ct
	param_S[SCORE_INDEX(GC, TC)] = param_S[SCORE_INDEX(CT, CG)] = ENTROPY(0.98f, GC_TC);	// Gt/Cc
	param_S[SCORE_INDEX(TA, CT)] = param_S[SCORE_INDEX(TC, AT)] = ENTROPY(0.97f, TA_CT);	// Tc/At
	param_S[SCORE_INDEX(TA, TC)] = param_S[SCORE_INDEX(CT, AT)] = ENTROPY(0.75f, TA_TC);	// Tt/Ac
	
	#define	AT_AC	2.3f
	#define	AT_CA	5.3f
	#define	CG_AC	1.9f
	#define	CG_CA	0.6f
	#define	GC_AC	5.2f
	#define	GC_CA	-0.7f
	#define	TA_AC	3.4f
	#define	TA_CA	7.6f
	
	// Single A-C mismatch
	// Allawi et. al. Biochemistry 1998, 37, 9435-9444
	param_H[SCORE_INDEX(AT, AC)] = param_H[SCORE_INDEX(CA, TA)] = AT_AC; 	// Aa/Tc
	param_H[SCORE_INDEX(AT, CA)] = param_H[SCORE_INDEX(AC, TA)] = AT_CA;	// Ac/Ta
	param_H[SCORE_INDEX(CG, AC)] = param_H[SCORE_INDEX(CA, GC)] = CG_AC;	// Ca/Gc
	param_H[SCORE_INDEX(CG, CA)] = param_H[SCORE_INDEX(AC, GC)] = CG_CA;	// Cc/Ga
	param_H[SCORE_INDEX(GC, AC)] = param_H[SCORE_INDEX(CA, CG)] = GC_AC;	// Ga/Cc
	param_H[SCORE_INDEX(GC, CA)] = param_H[SCORE_INDEX(AC, CG)] = GC_CA;	// Gc/Ca
	param_H[SCORE_INDEX(TA, AC)] = param_H[SCORE_INDEX(CA, AT)] = TA_AC;	// Ta/Ac
	param_H[SCORE_INDEX(TA, CA)] = param_H[SCORE_INDEX(AC, AT)] = TA_CA;	// Tc/Aa
	
	param_S[SCORE_INDEX(AT, AC)] = param_S[SCORE_INDEX(CA, TA)] = ENTROPY(0.88f, AT_AC); // Aa/Tc
	param_S[SCORE_INDEX(AT, CA)] = param_S[SCORE_INDEX(AC, TA)] = ENTROPY(0.77f, AT_CA);// Ac/Ta
	param_S[SCORE_INDEX(CG, AC)] = param_S[SCORE_INDEX(CA, GC)] = ENTROPY(0.75f, CG_AC);	// Ca/Gc
	param_S[SCORE_INDEX(CG, CA)] = param_S[SCORE_INDEX(AC, GC)] = ENTROPY(0.79f, CG_CA);// Cc/Ga
	param_S[SCORE_INDEX(GC, AC)] = param_S[SCORE_INDEX(CA, CG)] = ENTROPY(0.81f, GC_AC);// Ga/Cc
	param_S[SCORE_INDEX(GC, CA)] = param_S[SCORE_INDEX(AC, CG)] = ENTROPY(0.47f, GC_CA);// Gc/Ca
	param_S[SCORE_INDEX(TA, AC)] = param_S[SCORE_INDEX(CA, AT)] = ENTROPY(0.92f, TA_AC);	// Ta/Ac
	param_S[SCORE_INDEX(TA, CA)] = param_S[SCORE_INDEX(AC, AT)] = ENTROPY(1.33f, TA_CA);// Tc/Aa
	
	#define	AT_GT	1.0f
	#define	AT_TG	-2.5f
	#define	CG_GT	-4.1f
	#define	CG_TG	-2.8f
	#define	GC_GT	3.3f
	#define	GT_GT	5.8f
	#define	GC_TG	-4.4f
	#define	GT_TG	4.1f
	#define	TA_GT	-0.1f
	#define	TG_GT	-1.4f
	#define	TA_TG	-1.3f
	
	// Single G-T mismatch
	// Allawi et. al. Biochemistry 1997, 36, 10581-10594
	param_H[SCORE_INDEX(AT, GT)] = param_H[SCORE_INDEX(TG, TA)] = AT_GT;	// Ag/Tt
	param_H[SCORE_INDEX(AT, TG)] = param_H[SCORE_INDEX(GT, TA)] = AT_TG;	// At/Tg
	param_H[SCORE_INDEX(CG, GT)] = param_H[SCORE_INDEX(TG, GC)] = CG_GT;	// Cg/Gt
	param_H[SCORE_INDEX(CG, TG)] = param_H[SCORE_INDEX(GT, GC)] = CG_TG;	// Ct/Gg
	param_H[SCORE_INDEX(GC, GT)] = param_H[SCORE_INDEX(TG, CG)] = GC_GT;	// Gg/Ct
	param_H[SCORE_INDEX(GT, GT)] = param_H[SCORE_INDEX(TG, TG)] = GT_GT;	// gg/tt <-- double mismatch!
	param_H[SCORE_INDEX(GC, TG)] = param_H[SCORE_INDEX(GT, CG)] = GC_TG;	// Gt/Cg
	param_H[SCORE_INDEX(GT, TG)] 			    = GT_TG;	// gt/tg <-- double mismatch!
	param_H[SCORE_INDEX(TA, GT)] = param_H[SCORE_INDEX(TG, AT)] = TA_GT;	// Tg/At
	param_H[SCORE_INDEX(TG, GT)] 			    = TG_GT;	// tg/gt <-- double mismatch!
	param_H[SCORE_INDEX(TA, TG)] = param_H[SCORE_INDEX(GT, AT)] = TA_TG;	// Tt/Ag
	
	param_S[SCORE_INDEX(AT, GT)] = param_S[SCORE_INDEX(TG, TA)] = ENTROPY(0.71f, AT_GT);	// Ag/Tt
	param_S[SCORE_INDEX(AT, TG)] = param_S[SCORE_INDEX(GT, TA)] = ENTROPY(0.07f, AT_TG);	// At/Tg
	param_S[SCORE_INDEX(CG, GT)] = param_S[SCORE_INDEX(TG, GC)] = ENTROPY(-0.47f, CG_GT);	// Cg/Gt
	param_S[SCORE_INDEX(CG, TG)] = param_S[SCORE_INDEX(GT, GC)] = ENTROPY(-0.32f, CG_TG);	// Ct/Gg
	param_S[SCORE_INDEX(GC, GT)] = param_S[SCORE_INDEX(TG, CG)] = ENTROPY(0.08f, GC_GT);	// Gg/Ct
	param_S[SCORE_INDEX(GT, GT)] = param_S[SCORE_INDEX(TG, TG)] = ENTROPY(0.74f, GT_GT);	// gg/tt <-- double mismatch!
	param_S[SCORE_INDEX(GC, TG)] = param_S[SCORE_INDEX(GT, CG)] = ENTROPY(-0.59f, GC_TG);	// Gt/Cg
	param_S[SCORE_INDEX(GT, TG)] 			    = ENTROPY(1.15f, GT_TG);	// gt/tg <-- double mismatch!
	param_S[SCORE_INDEX(TA, GT)] = param_S[SCORE_INDEX(TG, AT)] = ENTROPY(0.43f, TA_GT);	// Tg/At
	param_S[SCORE_INDEX(TG, GT)] 			    = ENTROPY(0.52f, TG_GT);	// tg/gt <-- double mismatch!
	param_S[SCORE_INDEX(TA, TG)] = param_S[SCORE_INDEX(GT, AT)] = ENTROPY(0.34f, TA_TG);	// Tt/Ag
	
	#define	AT_AA	1.2f
	#define	CG_AA	-0.9f
	#define	GC_AA	-2.9f
	#define	TA_AA	4.7f
	
	// Single A-A, C-C, G-G or T-T mismatch
	// Peyret et. al. Biochemistry 1999, 38, 3468-3477
	param_H[SCORE_INDEX(AT, AA)] = param_H[SCORE_INDEX(AA, TA)] = AT_AA;	// Aa/Ta
	param_H[SCORE_INDEX(CG, AA)] = param_H[SCORE_INDEX(AA, GC)] = CG_AA;	// Ca/Ga
	param_H[SCORE_INDEX(GC, AA)] = param_H[SCORE_INDEX(AA, CG)] = GC_AA;	// Ga/Ca
	param_H[SCORE_INDEX(TA, AA)] = param_H[SCORE_INDEX(AA, AT)] = TA_AA;	// Ta/Aa
	
	param_S[SCORE_INDEX(AT, AA)] = param_S[SCORE_INDEX(AA, TA)] = ENTROPY(0.61f, AT_AA);	// Aa/Ta
	param_S[SCORE_INDEX(CG, AA)] = param_S[SCORE_INDEX(AA, GC)] = ENTROPY(0.43f, CG_AA);	// Ca/Ga
	param_S[SCORE_INDEX(GC, AA)] = param_S[SCORE_INDEX(AA, CG)] = ENTROPY(0.17f, GC_AA);	// Ga/Ca
	param_S[SCORE_INDEX(TA, AA)] = param_S[SCORE_INDEX(AA, AT)] = ENTROPY(0.69f, TA_AA);	// Ta/Aa
	
	#define	AT_CC	 0.0f
	#define	CG_CC	 -1.5f
	#define	GC_CC	 3.6f
	#define	TA_CC	 6.1f
	
	param_H[SCORE_INDEX(AT, CC)] = param_H[SCORE_INDEX(CC, TA)] = AT_CC;	// Ac/Tc
	param_H[SCORE_INDEX(CG, CC)] = param_H[SCORE_INDEX(CC, GC)] = CG_CC;	// Cc/Gc
	param_H[SCORE_INDEX(GC, CC)] = param_H[SCORE_INDEX(CC, CG)] = GC_CC;	// Gc/Cc
	param_H[SCORE_INDEX(TA, CC)] = param_H[SCORE_INDEX(CC, AT)] = TA_CC;	// Tc/Ac
	
	param_S[SCORE_INDEX(AT, CC)] = param_S[SCORE_INDEX(CC, TA)] = ENTROPY(1.33f, AT_CC);	// Ac/Tc
	param_S[SCORE_INDEX(CG, CC)] = param_S[SCORE_INDEX(CC, GC)] = ENTROPY(0.70f, CG_CC);	// Cc/Gc
	param_S[SCORE_INDEX(GC, CC)] = param_S[SCORE_INDEX(CC, CG)] = ENTROPY(0.79f, GC_CC);	// Gc/Cc
	param_S[SCORE_INDEX(TA, CC)] = param_S[SCORE_INDEX(CC, AT)] = ENTROPY(1.05f, TA_CC);	// Tc/Ac
	
	#define	AT_GG	 -3.1f
	#define	CG_GG	 -4.9f
	#define	GC_GG	 -6.0f
	#define	TA_GG	 1.6f
	
	param_H[SCORE_INDEX(AT, GG)] = param_H[SCORE_INDEX(GG, TA)] = AT_GG;	// Ag/Tg
	param_H[SCORE_INDEX(CG, GG)] = param_H[SCORE_INDEX(GG, GC)] = CG_GG;	// Cg/Gg
	param_H[SCORE_INDEX(GC, GG)] = param_H[SCORE_INDEX(GG, CG)] = GC_GG;	// Gg/Cg
	param_H[SCORE_INDEX(TA, GG)] = param_H[SCORE_INDEX(GG, AT)] = TA_GG;	// Tg/Ag
	
	param_S[SCORE_INDEX(AT, GG)] = param_S[SCORE_INDEX(GG, TA)] = ENTROPY(-0.13f, AT_GG);	// Ag/Tg
	param_S[SCORE_INDEX(CG, GG)] = param_S[SCORE_INDEX(GG, GC)] = ENTROPY(-0.11f, CG_GG);	// Cg/Gg
	param_S[SCORE_INDEX(GC, GG)] = param_S[SCORE_INDEX(GG, CG)] = ENTROPY(-1.11f, GC_GG);	// Gg/Cg
	param_S[SCORE_INDEX(TA, GG)] = param_S[SCORE_INDEX(GG, AT)] = ENTROPY(0.44f, TA_GG);	// Tg/Ag
	
	#define	AT_TT	 -2.7f
	#define	CG_TT	 -5.0f
	#define	GC_TT	 -2.2f
	#define	TA_TT	 0.2f
	
	param_H[SCORE_INDEX(AT, TT)] = param_H[SCORE_INDEX(TT, TA)] = AT_TT;	// At/Tt
	param_H[SCORE_INDEX(CG, TT)] = param_H[SCORE_INDEX(TT, GC)] = CG_TT;	// Ct/Gt
	param_H[SCORE_INDEX(GC, TT)] = param_H[SCORE_INDEX(TT, CG)] = GC_TT;	// Gt/Ct
	param_H[SCORE_INDEX(TA, TT)] = param_H[SCORE_INDEX(TT, AT)] = TA_TT;	// Tt/At
		
	param_S[SCORE_INDEX(AT, TT)] = param_S[SCORE_INDEX(TT, TA)] = ENTROPY(0.69f, AT_TT);	// At/Tt
	param_S[SCORE_INDEX(CG, TT)] = param_S[SCORE_INDEX(TT, GC)] = ENTROPY(-0.12f, CG_TT);	// Ct/Gt
	param_S[SCORE_INDEX(GC, TT)] = param_S[SCORE_INDEX(TT, CG)] = ENTROPY(0.45f, GC_TT);	// Gt/Ct
	param_S[SCORE_INDEX(TA, TT)] = param_S[SCORE_INDEX(TT, AT)] = ENTROPY(0.68f, TA_TT);	// Tt/At
	
	// Inosine base-pairs from
	// "Nearest-neighbor thermodynamics of deoxyinosine pairs in DNA duplexes"
	// Norman E. Watkins, Jr and John SantaLucia, Jr
	// Nucleic Acids Research, 2005, vol. 33, No. 19 pg. 6258-6267
	
	// Inosine IC pairs: Delta H
	#define	AT_IC	 -8.9f
	#define	TA_IC	 -5.9f
	#define	AT_CI	 -8.8f
	#define	TA_CI	 -4.9f
	#define	CG_IC	 -5.4f
	#define	GC_IC	 -6.8f
	#define	CG_CI	 -8.3f
	#define	GC_CI	 -5.0f
	
	param_H[SCORE_INDEX(AT, IC)] = param_H[SCORE_INDEX(CI, TA)] = AT_IC;	// AI/TC
	param_H[SCORE_INDEX(TA, IC)] = param_H[SCORE_INDEX(CI, AT)] = TA_IC;	// TI/AC
	param_H[SCORE_INDEX(AT, CI)] = param_H[SCORE_INDEX(IC, TA)] = AT_CI;	// AC/TI
	param_H[SCORE_INDEX(TA, CI)] = param_H[SCORE_INDEX(IC, AT)] = TA_CI;	// TC/AI
	param_H[SCORE_INDEX(CG, IC)] = param_H[SCORE_INDEX(CI, GC)] = CG_IC; 	// CI/GC
	param_H[SCORE_INDEX(GC, IC)] = param_H[SCORE_INDEX(CI, CG)] = GC_IC;	// GI/CC
	param_H[SCORE_INDEX(CG, CI)] = param_H[SCORE_INDEX(IC, GC)] = CG_CI;	// CC/GI
	param_H[SCORE_INDEX(GC, CI)] = param_H[SCORE_INDEX(IC, CG)] = GC_CI;	// GC/CI
	
	param_S[SCORE_INDEX(AT, IC)] = param_S[SCORE_INDEX(CI, TA)] = ENTROPY(-0.96f, AT_IC);	// AI/TC
	param_S[SCORE_INDEX(TA, IC)] = param_S[SCORE_INDEX(CI, AT)] = ENTROPY(-0.46f, TA_IC);	// AI/TC
	param_S[SCORE_INDEX(AT, CI)] = param_S[SCORE_INDEX(IC, TA)] = ENTROPY(-0.89f, AT_CI);	// AC/TI
	param_S[SCORE_INDEX(TA, CI)] = param_S[SCORE_INDEX(IC, AT)] = ENTROPY(-0.59f, TA_CI); // TC/AI
	param_S[SCORE_INDEX(CG, IC)] = param_S[SCORE_INDEX(CI, GC)] = ENTROPY(-1.14f, CG_IC);	// CI/GC
	param_S[SCORE_INDEX(GC, IC)] = param_S[SCORE_INDEX(CI, CG)] = ENTROPY(-0.86f, GC_IC);	// GI/CC
	param_S[SCORE_INDEX(CG, CI)] = param_S[SCORE_INDEX(IC, GC)] = ENTROPY(-0.88f, CG_CI);	// CC/GI
	param_S[SCORE_INDEX(GC, CI)] = param_S[SCORE_INDEX(IC, CG)] = ENTROPY(-1.07f, GC_CI);	// GC/CI
	
	// Inosine IA pairs: Delta H
	#define	AT_IA	 -8.3f
	#define	TA_IA	 -3.4f
	#define	AT_AI	 -0.7f
	#define	TA_AI	 -1.3f
	#define	CG_IA	  2.6f
	#define	GC_IA	 -7.8f
	#define	CG_AI	 -7.0f
	#define	GC_AI	 -7.6f
	
	param_H[SCORE_INDEX(AT, IA)] = param_H[SCORE_INDEX(AI, TA)] = AT_IA;	// AI/TA
	param_H[SCORE_INDEX(TA, IA)] = param_H[SCORE_INDEX(AI, AT)] = TA_IA;	// TI/AA
	param_H[SCORE_INDEX(AT, AI)] = param_H[SCORE_INDEX(IA, TA)] = AT_AI;	// AA/TI
	param_H[SCORE_INDEX(TA, AI)] = param_H[SCORE_INDEX(IA, AT)] = TA_AI;	// TA/AI
	param_H[SCORE_INDEX(CG, IA)] = param_H[SCORE_INDEX(AI, GC)] = CG_IA;	// CI/GA
	param_H[SCORE_INDEX(GC, IA)] = param_H[SCORE_INDEX(AI, CG)] = GC_IA;	// GI/CA
	param_H[SCORE_INDEX(CG, AI)] = param_H[SCORE_INDEX(IA, GC)] = CG_AI;   // CA/GI
	param_H[SCORE_INDEX(GC, AI)] = param_H[SCORE_INDEX(IA, CG)] = GC_AI;	// GA/CI
	
	param_S[SCORE_INDEX(AT, IA)] = param_S[SCORE_INDEX(AI, TA)] = ENTROPY(-0.51f, AT_IA);	// AI/TA
	param_S[SCORE_INDEX(TA, IA)] = param_S[SCORE_INDEX(AI, AT)] = ENTROPY(0.09f, TA_IA);	// TI/AA
	param_S[SCORE_INDEX(AT, AI)] = param_S[SCORE_INDEX(IA, TA)] = ENTROPY(0.12f, AT_AI);	// AA/TI
	param_S[SCORE_INDEX(TA, AI)] = param_S[SCORE_INDEX(IA, AT)] = ENTROPY(0.12f, TA_AI);	// TA/AI
	param_S[SCORE_INDEX(CG, IA)] = param_S[SCORE_INDEX(AI, GC)] = ENTROPY(-0.18f, CG_IA);	// CI/GA
	param_S[SCORE_INDEX(GC, IA)] = param_S[SCORE_INDEX(AI, CG)] = ENTROPY(-1.24f, GC_IA);	// GI/CA
	param_S[SCORE_INDEX(CG, AI)] = param_S[SCORE_INDEX(IA, GC)] = ENTROPY(-0.77f, CG_AI);	// CA/GI
	param_S[SCORE_INDEX(GC, AI)] = param_S[SCORE_INDEX(IA, CG)] = ENTROPY(-1.33f, GC_AI);	// GA/CI
	
	// Inosine IT pairs: Delta H
	#define	AT_IT	  0.49f
	#define	TA_IT	 -6.5f
	#define	AT_TI	 -5.6f
	#define	TA_TI	 -0.8f
	#define	CG_IT	 -1.0f
	#define	GC_IT	 -3.5f
	#define	CG_TI	  0.1f
	#define	GC_TI	 -4.3f
	
	param_H[SCORE_INDEX(AT, IT)] = param_H[SCORE_INDEX(TI, TA)] = AT_IT;	// AI/TT
	param_H[SCORE_INDEX(TA, IT)] = param_H[SCORE_INDEX(TI, AT)] = TA_IT;	// TI/AT
	param_H[SCORE_INDEX(AT, TI)] = param_H[SCORE_INDEX(IT, TA)] = AT_TI;	// AT/TI
	param_H[SCORE_INDEX(TA, TI)] = param_H[SCORE_INDEX(IT, AT)] = TA_TI;	// TT/AI
	param_H[SCORE_INDEX(CG, IT)] = param_H[SCORE_INDEX(TI, GC)] = CG_IT;	// CI/GT
	param_H[SCORE_INDEX(GC, IT)] = param_H[SCORE_INDEX(TI, CG)] = GC_IT;	// GI/CT
	param_H[SCORE_INDEX(CG, TI)] = param_H[SCORE_INDEX(IT, GC)] = CG_TI;	// CT/GI
	param_H[SCORE_INDEX(GC, TI)] = param_H[SCORE_INDEX(IT, CG)] = GC_TI;	// GT/CI
	
	param_S[SCORE_INDEX(AT, IT)] = param_S[SCORE_INDEX(TI, TA)] = ENTROPY(0.71f, AT_IT);	// AI/TT
	param_S[SCORE_INDEX(TA, IT)] = param_S[SCORE_INDEX(TI, AT)] = ENTROPY(0.36f, TA_IT);	// TI/AT
	param_S[SCORE_INDEX(AT, TI)] = param_S[SCORE_INDEX(IT, TA)] = ENTROPY(0.22f, AT_TI);	// AT/TI
	param_S[SCORE_INDEX(TA, TI)] = param_S[SCORE_INDEX(IT, AT)] = ENTROPY(0.54f, TA_TI);	// TT/AI
	param_S[SCORE_INDEX(CG, IT)] = param_S[SCORE_INDEX(TI, GC)] = ENTROPY(-0.26f, CG_IT);	// CI/GT
	param_S[SCORE_INDEX(GC, IT)] = param_S[SCORE_INDEX(TI, CG)] = ENTROPY(-0.19f, GC_IT);	// GI/CT
	param_S[SCORE_INDEX(CG, TI)] = param_S[SCORE_INDEX(IT, GC)] = ENTROPY(0.41f, CG_TI);	// CT/GI
	param_S[SCORE_INDEX(GC, TI)] = param_S[SCORE_INDEX(IT, CG)] = ENTROPY(-0.54f, GC_TI);	// GT/CI

	// Inosine IG pairs: Delta H
	#define	AT_IG	  -4.9f
	#define	TA_IG	  -1.9f
	#define	AT_GI	   0.1f
	#define	TA_GI	   1.0f
	#define	CG_IG	   7.1f
	#define	GC_IG	  -1.1f
	#define	CG_GI	   5.8f
	#define	GC_GI	  -7.6f
	
	param_H[SCORE_INDEX(AT, IG)] = param_H[SCORE_INDEX(GI, TA)] = AT_IG;	// AI/TG
	param_H[SCORE_INDEX(TA, IG)] = param_H[SCORE_INDEX(GI, AT)] = TA_IG;	// TI/AG
	param_H[SCORE_INDEX(AT, GI)] = param_H[SCORE_INDEX(IG, TA)] = AT_GI;	// AG/TI
	param_H[SCORE_INDEX(TA, GI)] = param_H[SCORE_INDEX(IG, AT)] = TA_GI;	// TG/AI
	param_H[SCORE_INDEX(CG, IG)] = param_H[SCORE_INDEX(GI, GC)] = CG_IG;	// CI/GG
	param_H[SCORE_INDEX(GC, IG)] = param_H[SCORE_INDEX(GI, CG)] = GC_IG;	// GI/CG
	param_H[SCORE_INDEX(CG, GI)] = param_H[SCORE_INDEX(IG, GC)] = CG_GI;	// CG/GI
	param_H[SCORE_INDEX(GC, GI)] = param_H[SCORE_INDEX(IG, CG)] = GC_GI;	// GG/CI
	
	param_S[SCORE_INDEX(AT, IG)] = param_S[SCORE_INDEX(GI, TA)] = ENTROPY(0.02f, AT_IG);	// AI/TG
	param_S[SCORE_INDEX(TA, IG)] = param_S[SCORE_INDEX(GI, AT)] = ENTROPY(0.76f, TA_IG);	// TI/AG
	param_S[SCORE_INDEX(AT, GI)] = param_S[SCORE_INDEX(IG, TA)] = ENTROPY(0.65f, AT_GI);	// AG/TI
	param_S[SCORE_INDEX(TA, GI)] = param_S[SCORE_INDEX(IG, AT)] = ENTROPY(0.70f, TA_GI);	// TG/AI
	param_S[SCORE_INDEX(CG, IG)] = param_S[SCORE_INDEX(GI, GC)] = ENTROPY(0.47f, CG_IG);	// CI/GG
	param_S[SCORE_INDEX(GC, IG)] = param_S[SCORE_INDEX(GI, CG)] = ENTROPY(-0.10f, GC_IG);	// GI/CG
	param_S[SCORE_INDEX(CG, GI)] = param_S[SCORE_INDEX(IG, GC)] = ENTROPY(0.54f, CG_GI);	// CG/GI
	param_S[SCORE_INDEX(GC, GI)] = param_S[SCORE_INDEX(IG, CG)] = ENTROPY(-0.74f, GC_GI);	// GG/CI
	
	// Inosine II pairs: Delta H
	#define	AT_II	  -3.3f
	#define	TA_II	   0.1f
	#define	CG_II	   1.3f
	#define	GC_II	  -0.5f
	
	param_H[SCORE_INDEX(AT, II)] = param_H[SCORE_INDEX(II, TA)] = AT_II;	// AI/TI
	param_H[SCORE_INDEX(TA, II)] = param_H[SCORE_INDEX(II, AT)] = TA_II;	// TI/AI
	param_H[SCORE_INDEX(CG, II)] = param_H[SCORE_INDEX(II, GC)] = CG_II;	// CI/GI
	param_H[SCORE_INDEX(GC, II)] = param_H[SCORE_INDEX(II, CG)] = GC_II;	// GI/CI
	
	param_S[SCORE_INDEX(AT, II)] = param_S[SCORE_INDEX(II, TA)] = ENTROPY(0.40f, AT_II);	// AI/TI
	param_S[SCORE_INDEX(TA, II)] = param_S[SCORE_INDEX(II, AT)] = ENTROPY(0.81f, TA_II);	// TI/AI
	param_S[SCORE_INDEX(CG, II)] = param_S[SCORE_INDEX(II, GC)] = ENTROPY(0.36f, CG_II);	// CI/GI
	param_S[SCORE_INDEX(GC, II)] = param_S[SCORE_INDEX(II, CG)] = ENTROPY(-0.09f, GC_II);	// GI/CI
	
	// Inosine tandem internal pairs: Delta H
	#define	IC_IC	-9.3f
	#define	IA_IC	-3.1f
	#define	IC_IA	-8.7f
	#define	IA_IA	-2.1f
	#define	IT_IA	 2.3f
	#define	IG_IA	 4.2f
	#define	IC_IT	-14.5f
	#define	IA_IT	-17.8f
	#define	IT_IT	-7.0f
	#define	IG_IT	-19.4f
	#define	IT_IG	 13.3f
	#define	IG_IG	 0.3f	// There is a typo in the paper (II/TG should read II/GG)
	#define	II_II	-10.65f	// Note that II/II is *not* a two state reaction!
							// Two values are published for II/II; -13.8 and -7.5. For
							// now, use the average -> -10.65
	
	param_H[SCORE_INDEX(IC, IC)] = param_H[SCORE_INDEX(CI, CI)] = IC_IC;	// II/CC
	param_H[SCORE_INDEX(IA, IC)] = param_H[SCORE_INDEX(CI, AI)] = IA_IC;	// II/AC
	param_H[SCORE_INDEX(IC, IA)] = param_H[SCORE_INDEX(AI, CI)] = IC_IA;	// II/CA
	param_H[SCORE_INDEX(IA, IA)] = param_H[SCORE_INDEX(AI, AI)] = IA_IA;	// II/AA
	param_H[SCORE_INDEX(IT, IA)] = param_H[SCORE_INDEX(AI, TI)] = IT_IA;	// II/TA
	param_H[SCORE_INDEX(IG, IA)] = param_H[SCORE_INDEX(AI, GI)] = IG_IA;	// II/GA
	param_H[SCORE_INDEX(IC, IT)] = param_H[SCORE_INDEX(TI, CI)] = IC_IT;	// II/CT
	param_H[SCORE_INDEX(IA, IT)] = param_H[SCORE_INDEX(TI, AI)] = IA_IT;	// II/AT
	param_H[SCORE_INDEX(IT, IT)] = param_H[SCORE_INDEX(TI, TI)] = IT_IT;	// II/TT
	param_H[SCORE_INDEX(IG, IT)] = param_H[SCORE_INDEX(TI, GI)] = IG_IT;	// II/GT
	param_H[SCORE_INDEX(IT, IG)] = param_H[SCORE_INDEX(GI, TI)] = IT_IG;	// II/TG
	param_H[SCORE_INDEX(IG, IG)] = param_H[SCORE_INDEX(GI, GI)] = IG_IG;	// II/GG
	param_H[SCORE_INDEX(II, II)]                   = II_II;	// II/II
	
	param_S[SCORE_INDEX(IC, IC)] = param_S[SCORE_INDEX(CI, CI)] = ENTROPY(-0.64f, IC_IC);	// II/CC
	param_S[SCORE_INDEX(IA, IC)] = param_S[SCORE_INDEX(CI, AI)] = ENTROPY(0.27f, IA_IC);	// II/AC
	param_S[SCORE_INDEX(IC, IA)] = param_S[SCORE_INDEX(AI, CI)] = ENTROPY(0.44f, IC_IA);	// II/CA
	param_S[SCORE_INDEX(IA, IA)] = param_S[SCORE_INDEX(AI, AI)] = ENTROPY(-0.27f, IA_IA);	// II/AA
	param_S[SCORE_INDEX(IT, IA)] = param_S[SCORE_INDEX(AI, TI)] = ENTROPY(0.83f, IT_IA);	// II/TA
	param_S[SCORE_INDEX(IG, IA)] = param_S[SCORE_INDEX(AI, GI)] = ENTROPY(0.30f, IG_IA);	// II/GA
	param_S[SCORE_INDEX(IC, IT)] = param_S[SCORE_INDEX(TI, CI)] = ENTROPY(0.33f, IC_IT);	// II/CT
	param_S[SCORE_INDEX(IA, IT)] = param_S[SCORE_INDEX(TI, AI)] = ENTROPY(0.19f, IA_IT);	// II/AT
	param_S[SCORE_INDEX(IT, IT)] = param_S[SCORE_INDEX(TI, TI)] = ENTROPY(1.69f, IT_IT);	// II/TT
	param_S[SCORE_INDEX(IG, IT)] = param_S[SCORE_INDEX(TI, GI)] = ENTROPY(0.13f, IG_IT);	// II/GT
	param_S[SCORE_INDEX(IT, IG)] = param_S[SCORE_INDEX(GI, TI)] = ENTROPY(0.03f, IT_IG);	// II/TG
	param_S[SCORE_INDEX(IG, IG)] = param_S[SCORE_INDEX(GI, GI)] = ENTROPY(-1.30f, IG_IG);	// II/GG
	param_S[SCORE_INDEX(II, II)]                   = ENTROPY(0.52f, II_II);	// II/II; Not two state!
									// Use the average DG:
									// (-0.7 and 1.74)
	
	// Inosine "other" tandem mismatch pairs: Delta H
	#define	IC_CI	-12.1f
	#define	CI_IC	-1.8f
	#define	IA_AI	-13.9f
	#define	AI_IA	-9.5f
	#define	IT_TI	-7.6f
	#define	TI_IT	-14.7f
	#define	IG_GI	 3.2f
	#define	GI_IG	-4.2f
	
	param_H[SCORE_INDEX(IC, CI)] = IC_CI;	// IC/CI
	param_H[SCORE_INDEX(CI, IC)] = CI_IC;	// CI/IC
	param_H[SCORE_INDEX(IA, AI)] = IA_AI;	// IA/AI
	param_H[SCORE_INDEX(AI, IA)] = AI_IA;	// AI/IA
	param_H[SCORE_INDEX(IT, TI)] = IT_TI;	// IT/TI
	param_H[SCORE_INDEX(TI, IT)] = TI_IT;	// TI/IT
	param_H[SCORE_INDEX(IG, GI)] = IG_GI;	// IG/GI
	param_H[SCORE_INDEX(GI, IG)] = GI_IG;	// GI/IG
	
	param_S[SCORE_INDEX(IC, CI)] =  ENTROPY(-0.85f, IC_CI);	// IC/CI
	param_S[SCORE_INDEX(CI, IC)] =  ENTROPY(0.06f, CI_IC);	// CI/IC
	param_S[SCORE_INDEX(IA, AI)] =  ENTROPY(-1.43f, IA_AI);	// IA/AI
	param_S[SCORE_INDEX(AI, IA)] =  ENTROPY(-0.56f, AI_IA);	// AI/IA
	param_S[SCORE_INDEX(IT, TI)] =  ENTROPY(2.03f, IT_TI);	// IT/TI
	param_S[SCORE_INDEX(TI, IT)] =  ENTROPY(0.61f, TI_IT);	// TI/IT
	param_S[SCORE_INDEX(IG, GI)] =  ENTROPY(1.18f, IG_GI);	// IG/GI
	param_S[SCORE_INDEX(GI, IG)] =  ENTROPY(1.12f, GI_IG);	// GI/IG
	
	///////////////////////////////////////////////////////////////////////
	// Internal loop and hairpin terminal mismatch parameters
	// These parameters have not been officially published. However, the parameters
	// are available from the unafold software package 
	// (http://frontend.bioinfo.rpi.edu/applications/hybrid/download.php).
	memcpy( param_loop_terminal_H, param_H, NUM_BASE_PAIR*NUM_BASE_PAIR*sizeof(float) );
	memcpy( param_loop_terminal_S, param_S, NUM_BASE_PAIR*NUM_BASE_PAIR*sizeof(float) );
	
	memcpy( param_hairpin_terminal_H, param_H, NUM_BASE_PAIR*NUM_BASE_PAIR*sizeof(float) );
	memcpy( param_hairpin_terminal_S, param_S, NUM_BASE_PAIR*NUM_BASE_PAIR*sizeof(float) );
	
	// ******************************************************************************
	// Include a subset of unpublished parameters from the unafold program.
	// Please note: Since unafold is not open source software, we can not 
	// distribute this code outside of LANL!
	// ******************************************************************************
	#include "nuc_cruc_santa_lucia_tstacki.cpp" // <-- Empty file
	#include "nuc_cruc_santa_lucia_tstackh.cpp" // <-- Empty file
	
	// The following parameters are from:
	// "The thermodynamics of DNA structural motifs"
	// J. SantaLucia and D. Hicks
	// Annu. Rev. Biophys. Biomol. Struc. 2004. 33:415-440
	param_init_H = 0.2f;
	param_init_S = ENTROPY(1.96f, 0.2f);
	
	param_AT_closing_H = 2.2f;
	param_AT_closing_S = ENTROPY(0.05f, 2.2f);
		
	param_symmetry_S = ENTROPY(0.43f, 0.0f);
	
	param_SALT = 0.368e-3f;
	
	param_asymmetric_loop_dS = ENTROPY(0.3f, 0.0f);
	
	param_bulge_AT_closing_S = ENTROPY(0.5f, 0.0f);
	
	///////////////////////////////////////////////////////
	// From Table 4 in the above reference:
	///////////////////////////////////////////////////////
	// Internal Loop
	param_loop_S[0] = 0.0f;
	param_loop_S[1] = 0.0f;
	param_loop_S[2] = 0.0f;
	param_loop_S[3] = ENTROPY(3.2f, 0.0f);
	param_loop_S[4] = ENTROPY(3.6f, 0.0f);
	param_loop_S[5] = ENTROPY(4.0f, 0.0f);
	param_loop_S[6] = ENTROPY(4.4f, 0.0f);
	param_loop_S[7] = ENTROPY(4.6f, 0.0f);
	param_loop_S[8] = ENTROPY(4.8f, 0.0f);
	param_loop_S[9] = ENTROPY(4.9f, 0.0f);
	param_loop_S[10] = ENTROPY(4.9f, 0.0f);
	
	param_loop_S[12] = ENTROPY(5.2f, 0.0f);
	
	param_loop_S[14] = ENTROPY(5.4f, 0.0f);
	
	param_loop_S[16] = ENTROPY(5.6f, 0.0f);
	
	param_loop_S[18] = ENTROPY(5.8f, 0.0f);
	
	param_loop_S[20] = ENTROPY(5.9f, 0.0f);
	
	param_loop_S[25] = ENTROPY(6.3f, 0.0f);
	
	param_loop_S[30] = ENTROPY(6.6f, 0.0f);
	
	param_loop_S[11] = INTERPOLATE(param_loop_S, 10, 12, 11);
	param_loop_S[13] = INTERPOLATE(param_loop_S, 12, 14, 13);
	param_loop_S[15] = INTERPOLATE(param_loop_S, 14, 16, 15);
	param_loop_S[17] = INTERPOLATE(param_loop_S, 16, 18, 17);
	param_loop_S[19] = INTERPOLATE(param_loop_S, 18, 20, 19);
	param_loop_S[21] = INTERPOLATE(param_loop_S, 20, 25, 21);
	param_loop_S[22] = INTERPOLATE(param_loop_S, 20, 25, 22);
	param_loop_S[23] = INTERPOLATE(param_loop_S, 20, 25, 23);
	param_loop_S[24] = INTERPOLATE(param_loop_S, 20, 25, 24);
	param_loop_S[26] = INTERPOLATE(param_loop_S, 25, 30, 26);
	param_loop_S[27] = INTERPOLATE(param_loop_S, 25, 30, 27);
	param_loop_S[28] = INTERPOLATE(param_loop_S, 25, 30, 28);
	param_loop_S[29] = INTERPOLATE(param_loop_S, 25, 30, 29);
	
	// Jacobson-Stockmayer entropy extrapolation
	// dS(n) = dS(x) - 2.44*R*ln(n/x)
	#define ANCHOR_LOOP_EXTRAP	30
	for(unsigned int i = 31;i < MAX_LOOP_LENGTH;i++){
		param_loop_S[i] = param_loop_S[ANCHOR_LOOP_EXTRAP] - 2.44f*NC_R*( log( double(i)/ANCHOR_LOOP_EXTRAP) );
	}
	
	///////////////////////////////////////////////////////
	// Bulge
	param_bulge_S[0] = 0.0f;
	param_bulge_S[1] = ENTROPY(4.0f, 0.0f);
	param_bulge_S[2] = ENTROPY(2.9f, 0.0f);
	param_bulge_S[3] = ENTROPY(3.1f, 0.0f);
	param_bulge_S[4] = ENTROPY(3.2f, 0.0f);
	param_bulge_S[5] = ENTROPY(3.3f, 0.0f);
	param_bulge_S[6] = ENTROPY(3.5f, 0.0f);
	param_bulge_S[7] = ENTROPY(3.7f, 0.0f);
	param_bulge_S[8] = ENTROPY(3.9f, 0.0f);
	param_bulge_S[9] = ENTROPY(4.1f, 0.0f);
	param_bulge_S[10] = ENTROPY(4.3f, 0.0f);
	
	param_bulge_S[12] = ENTROPY(4.5f, 0.0f);
	
	param_bulge_S[14] = ENTROPY(4.8f, 0.0f);
	
	param_bulge_S[16] = ENTROPY(5.0f, 0.0f);
	
	param_bulge_S[18] = ENTROPY(5.2f, 0.0f);
	
	param_bulge_S[20] = ENTROPY(5.3f, 0.0f);
	
	param_bulge_S[25] = ENTROPY(5.6f, 0.0f);
	
	param_bulge_S[30] = ENTROPY(5.9f, 0.0f);
	
	param_bulge_S[11] = INTERPOLATE(param_bulge_S, 10, 12, 11);
	param_bulge_S[13] = INTERPOLATE(param_bulge_S, 12, 14, 13);
	param_bulge_S[15] = INTERPOLATE(param_bulge_S, 14, 16, 15);
	param_bulge_S[17] = INTERPOLATE(param_bulge_S, 16, 18, 17);
	param_bulge_S[19] = INTERPOLATE(param_bulge_S, 18, 20, 19);
	param_bulge_S[21] = INTERPOLATE(param_bulge_S, 20, 25, 21);
	param_bulge_S[22] = INTERPOLATE(param_bulge_S, 20, 25, 22);
	param_bulge_S[23] = INTERPOLATE(param_bulge_S, 20, 25, 23);
	param_bulge_S[24] = INTERPOLATE(param_bulge_S, 20, 25, 24);
	param_bulge_S[26] = INTERPOLATE(param_bulge_S, 25, 30, 26);
	param_bulge_S[27] = INTERPOLATE(param_bulge_S, 25, 30, 27);
	param_bulge_S[28] = INTERPOLATE(param_bulge_S, 25, 30, 28);
	param_bulge_S[29] = INTERPOLATE(param_bulge_S, 25, 30, 29);
	
	// Jacobson-Stockmayer entropy extrapolation
	// dS(n) = dS(x) - 2.44*R*ln(n/x)
	#define ANCHOR_BULGE_EXTRAP	30
	
	for(unsigned int i = 31;i < MAX_BULGE_LENGTH;i++){
		param_bulge_S[i] = param_bulge_S[ANCHOR_BULGE_EXTRAP] - 2.44f*NC_R*( log( double(i)/ANCHOR_BULGE_EXTRAP) );
	}
	
	///////////////////////////////////////////////////////
	// Hairpin loop
	param_hairpin_S[0] = 0.0f;
	param_hairpin_S[1] = 0.0f;
	param_hairpin_S[2] = 0.0f;
	param_hairpin_S[3] = ENTROPY(3.5f, 0.0f);
	param_hairpin_S[4] = ENTROPY(3.5f, 0.0f);
	param_hairpin_S[5] = ENTROPY(3.3f, 0.0f);
	param_hairpin_S[6] = ENTROPY(4.0f, 0.0f);
	param_hairpin_S[7] = ENTROPY(4.2f, 0.0f);
	param_hairpin_S[8] = ENTROPY(4.3f, 0.0f);
	param_hairpin_S[9] = ENTROPY(4.5f, 0.0f);
	param_hairpin_S[10] = ENTROPY(4.6f, 0.0f);
	
	param_hairpin_S[12] = ENTROPY(5.0f, 0.0f);
	
	param_hairpin_S[14] = ENTROPY(5.1f, 0.0f);
	
	param_hairpin_S[16] = ENTROPY(5.3f, 0.0f);
	
	param_hairpin_S[18] = ENTROPY(5.5f, 0.0f);
	
	param_hairpin_S[20] = ENTROPY(5.7f, 0.0f);
	
	param_hairpin_S[25] = ENTROPY(6.1f, 0.0f);
	
	param_hairpin_S[30] = ENTROPY(6.3f, 0.0f);
	
	// Interpolate the missing values
	param_hairpin_S[11] = INTERPOLATE(param_hairpin_S, 10, 12, 11);
	param_hairpin_S[13] = INTERPOLATE(param_hairpin_S, 12, 14, 13);
	param_hairpin_S[15] = INTERPOLATE(param_hairpin_S, 14, 16, 15);
	param_hairpin_S[17] = INTERPOLATE(param_hairpin_S, 16, 18, 17);
	param_hairpin_S[19] = INTERPOLATE(param_hairpin_S, 18, 20, 19);
	param_hairpin_S[21] = INTERPOLATE(param_hairpin_S, 20, 25, 21);
	param_hairpin_S[22] = INTERPOLATE(param_hairpin_S, 20, 25, 22);
	param_hairpin_S[23] = INTERPOLATE(param_hairpin_S, 20, 25, 23);
	param_hairpin_S[24] = INTERPOLATE(param_hairpin_S, 20, 25, 24);
	param_hairpin_S[26] = INTERPOLATE(param_hairpin_S, 25, 30, 26);
	param_hairpin_S[27] = INTERPOLATE(param_hairpin_S, 25, 30, 27);
	param_hairpin_S[28] = INTERPOLATE(param_hairpin_S, 25, 30, 28);
	param_hairpin_S[29] = INTERPOLATE(param_hairpin_S, 25, 30, 29);
	
	// Jacobson-Stockmayer entropy extrapolation
	// dS(n) = dS(x) - 2.44*R*ln(n/x)
	#define ANCHOR_HAIRPIN_EXTRAP	30
	
	for(unsigned int i = 31;i < MAX_HAIRPIN_LENGTH;i++){
		param_hairpin_S[i] = param_hairpin_S[ANCHOR_HAIRPIN_EXTRAP] - 2.44f*NC_R*( log( double(i)/ANCHOR_HAIRPIN_EXTRAP) );
	}
	
	// From:
	// Supplemental Material: Annu.Rev.Biophs.Biomol.Struct.33:415-40
	// doi: 10.1146/annurev.biophys.32.110601.141800
	// The Termodynamicso f DNA Structural Motifs
	// SantaLucia and Hicks, 2004
	
	// The table has been saved as a text file in: nuc_cruc_santa_lucia_data.txt
	// The code below was generated with the command:
	// cat nuc_cruc_santa_lucia_data.txt | awk '{print "param_hairpin_special_H["$1"] = "$3";\nparam_hairpin_special_S["$1"] = " ($3-$2)/310.15";\n"}'
	//
	// This command does *not* append "f" to each number to avoid compiler errors when the number
	// is printed as an integer.
	
	param_hairpin_special_H[AAAAAT] = 0.5;
	param_hairpin_special_S[AAAAAT] = -0.000644849;

	param_hairpin_special_H[AAAACT] = 0.7;
	param_hairpin_special_S[AAAACT] = 0.00161212;

	param_hairpin_special_H[AAACAT] = 1;
	param_hairpin_special_S[AAACAT] = 0.00161212;

	param_hairpin_special_H[ACTTGT] = 0;
	param_hairpin_special_S[ACTTGT] = 0.00419152;

	param_hairpin_special_H[AGAAAT] = -1.1;
	param_hairpin_special_S[AGAAAT] = 0.00161212;

	param_hairpin_special_H[AGAAT] = -1.5;
	param_hairpin_special_S[AGAAT] = 0;

	param_hairpin_special_H[AGAGAT] = -1.1;
	param_hairpin_special_S[AGAGAT] = 0.00161212;

	param_hairpin_special_H[AGATAT] = -1.5;
	param_hairpin_special_S[AGATAT] = 0.00161212;

	param_hairpin_special_H[AGCAAT] = -1.6;
	param_hairpin_special_S[AGCAAT] = 0.00161212;

	param_hairpin_special_H[AGCAT] = -1.5;
	param_hairpin_special_S[AGCAT] = 0;

	param_hairpin_special_H[AGCGAT] = -1.1;
	param_hairpin_special_S[AGCGAT] = 0.00161212;

	param_hairpin_special_H[AGCTTT] = 0.2;
	param_hairpin_special_S[AGCTTT] = 0.00161212;

	param_hairpin_special_H[AGGAAT] = -1.1;
	param_hairpin_special_S[AGGAAT] = 0.00161212;

	param_hairpin_special_H[AGGAT] = -1.5;
	param_hairpin_special_S[AGGAT] = 0;

	param_hairpin_special_H[AGGGAT] = -1.1;
	param_hairpin_special_S[AGGGAT] = 0.00161212;

	param_hairpin_special_H[AGGGGT] = 0.5;
	param_hairpin_special_S[AGGGGT] = 0.000644849;

	param_hairpin_special_H[AGTAAT] = -1.6;
	param_hairpin_special_S[AGTAAT] = 0.00161212;

	param_hairpin_special_H[AGTAT] = -1.5;
	param_hairpin_special_S[AGTAT] = 0;

	param_hairpin_special_H[AGTGAT] = -1.1;
	param_hairpin_special_S[AGTGAT] = 0.00161212;

	param_hairpin_special_H[AGTTCT] = 0.8;
	param_hairpin_special_S[AGTTCT] = 0.00161212;

	param_hairpin_special_H[ATTCGT] = -0.2;
	param_hairpin_special_S[ATTCGT] = 0.00161212;

	param_hairpin_special_H[ATTTGT] = 0;
	param_hairpin_special_S[ATTTGT] = 0.00161212;

	param_hairpin_special_H[ATTTTT] = -0.5;
	param_hairpin_special_S[ATTTTT] = 0.00161212;

	param_hairpin_special_H[CAAAAG] = 0.5;
	param_hairpin_special_S[CAAAAG] = -0.0012897;

	param_hairpin_special_H[CAAACG] = 0.7;
	param_hairpin_special_S[CAAACG] = 0;

	param_hairpin_special_H[CAACAG] = 1;
	param_hairpin_special_S[CAACAG] = 0;

	param_hairpin_special_H[CAACCG] = 0;
	param_hairpin_special_S[CAACCG] = 0;

	param_hairpin_special_H[CCTTGG] = 0;
	param_hairpin_special_S[CCTTGG] = 0.0025794;

	param_hairpin_special_H[CGAAAG] = -1.1;
	param_hairpin_special_S[CGAAAG] = 0;

	param_hairpin_special_H[CGAAG] = -2.0;
	param_hairpin_special_S[CGAAG] = 0;

	param_hairpin_special_H[CGAGAG] = -1.1;
	param_hairpin_special_S[CGAGAG] = 0;

	param_hairpin_special_H[CGATAG] = -1.5;
	param_hairpin_special_S[CGATAG] = 0;

	param_hairpin_special_H[CGCAAG] = -1.6;
	param_hairpin_special_S[CGCAAG] = 0;

	param_hairpin_special_H[CGCAG] = -2.0;
	param_hairpin_special_S[CGCAG] = 0;

	param_hairpin_special_H[CGCGAG] = -1.1;
	param_hairpin_special_S[CGCGAG] = 0;

	param_hairpin_special_H[CGCTTG] = 0.2;
	param_hairpin_special_S[CGCTTG] = 0;

	param_hairpin_special_H[CGGAAG] = -1.1;
	param_hairpin_special_S[CGGAAG] = 0;

	param_hairpin_special_H[CGGAG] = -2.0;
	param_hairpin_special_S[CGGAG] = 0;

	param_hairpin_special_H[CGGGAG] = -1;
	param_hairpin_special_S[CGGGAG] = 0;

	param_hairpin_special_H[CGGGGG] = 0.5;
	param_hairpin_special_S[CGGGGG] = -0.000967274;

	param_hairpin_special_H[CGTAAG] = -1.6;
	param_hairpin_special_S[CGTAAG] = 0;

	param_hairpin_special_H[CGTAG] = -2.0;
	param_hairpin_special_S[CGTAG] = 0;

	param_hairpin_special_H[CGTGAG] = -1.1;
	param_hairpin_special_S[CGTGAG] = 0;

	param_hairpin_special_H[CGTTCG] = 0.8;
	param_hairpin_special_S[CGTTCG] = 0;

	param_hairpin_special_H[CTTCGG] = -0.2;
	param_hairpin_special_S[CTTCGG] = 0;

	param_hairpin_special_H[CTTTGG] = 0;
	param_hairpin_special_S[CTTTGG] = 0;

	param_hairpin_special_H[CTTTTG] = -0.5;
	param_hairpin_special_S[CTTTTG] = 0;

	param_hairpin_special_H[GAAAAC] = 0.5;
	param_hairpin_special_S[GAAAAC] = -0.00322425;

	param_hairpin_special_H[GAAAAT] = 0.5;
	param_hairpin_special_S[GAAAAT] = -0.00322425;

	param_hairpin_special_H[GAAACC] = 0.7;
	param_hairpin_special_S[GAAACC] = 0;

	param_hairpin_special_H[GAAACT] = 1;
	param_hairpin_special_S[GAAACT] = 0;

	param_hairpin_special_H[GAACAC] = 1;
	param_hairpin_special_S[GAACAC] = 0;

	param_hairpin_special_H[GAACAT] = 1;
	param_hairpin_special_S[GAACAT] = 0;

	param_hairpin_special_H[GCTTGC] = 0;
	param_hairpin_special_S[GCTTGC] = 0.0025794;

	param_hairpin_special_H[GCTTGT] = 0;
	param_hairpin_special_S[GCTTGT] = 0.00161212;

	param_hairpin_special_H[GGAAAC] = -1.1;
	param_hairpin_special_S[GGAAAC] = 0;

	param_hairpin_special_H[GGAAAT] = -1.1;
	param_hairpin_special_S[GGAAAT] = 0;

	param_hairpin_special_H[GGAAC] = -2.0;
	param_hairpin_special_S[GGAAC] = 0;

	param_hairpin_special_H[GGAGAC] = -1.1;
	param_hairpin_special_S[GGAGAC] = 0;

	param_hairpin_special_H[GGAGAT] = -1.1;
	param_hairpin_special_S[GGAGAT] = 0;

	param_hairpin_special_H[GGATAC] = -1.6;
	param_hairpin_special_S[GGATAC] = 0;

	param_hairpin_special_H[GGATAT] = -1.6;
	param_hairpin_special_S[GGATAT] = 0;

	param_hairpin_special_H[GGCAAC] = -1.6;
	param_hairpin_special_S[GGCAAC] = 0;

	param_hairpin_special_H[GGCAAT] = -1.6;
	param_hairpin_special_S[GGCAAT] = 0;

	param_hairpin_special_H[GGCAC] = -2.0;
	param_hairpin_special_S[GGCAC] = 0;

	param_hairpin_special_H[GGCGAC] = -1.1;
	param_hairpin_special_S[GGCGAC] = 0;

	param_hairpin_special_H[GGCGAT] = -1.1;
	param_hairpin_special_S[GGCGAT] = 0;

	param_hairpin_special_H[GGCTTC] = 0.2;
	param_hairpin_special_S[GGCTTC] = 0;

	param_hairpin_special_H[GGCTTT] = -0.1;
	param_hairpin_special_S[GGCTTT] = 0;

	param_hairpin_special_H[GGGAAC] = -1.1;
	param_hairpin_special_S[GGGAAC] = 0;

	param_hairpin_special_H[GGGAAT] = -1.1;
	param_hairpin_special_S[GGGAAT] = 0;

	param_hairpin_special_H[GGGAC] = -2.0;
	param_hairpin_special_S[GGGAC] = 0;

	param_hairpin_special_H[GGGGAC] = -1.1;
	param_hairpin_special_S[GGGGAC] = 0;

	param_hairpin_special_H[GGGGAT] = -1.1;
	param_hairpin_special_S[GGGGAT] = 0;

	param_hairpin_special_H[GGGGGC] = 0.5;
	param_hairpin_special_S[GGGGGC] = -0.000967274;

	param_hairpin_special_H[GGGGGT] = 0.5;
	param_hairpin_special_S[GGGGGT] = -0.000967274;

	param_hairpin_special_H[GGTAAC] = -1.6;
	param_hairpin_special_S[GGTAAC] = 0;

	param_hairpin_special_H[GGTAAT] = -1.6;
	param_hairpin_special_S[GGTAAT] = 0;

	param_hairpin_special_H[GGTAC] = -2.0;
	param_hairpin_special_S[GGTAC] = 0;

	param_hairpin_special_H[GGTGAC] = -1.1;
	param_hairpin_special_S[GGTGAC] = 0;

	param_hairpin_special_H[GGTGAT] = -1.1;
	param_hairpin_special_S[GGTGAT] = 0;

	param_hairpin_special_H[GGTTCC] = 0.8;
	param_hairpin_special_S[GGTTCC] = 0;

	param_hairpin_special_H[GTATAT] = -0.5;
	param_hairpin_special_S[GTATAT] = 0;

	param_hairpin_special_H[GTTCGC] = -0.2;
	param_hairpin_special_S[GTTCGC] = 0;

	param_hairpin_special_H[GTTCGT] = -0.4;
	param_hairpin_special_S[GTTCGT] = 0;

	param_hairpin_special_H[GTTTGC] = 0;
	param_hairpin_special_S[GTTTGC] = 0;

	param_hairpin_special_H[GTTTGT] = -0.4;
	param_hairpin_special_S[GTTTGT] = 0;

	param_hairpin_special_H[GTTTTC] = -0.5;
	param_hairpin_special_S[GTTTTC] = 0;

	param_hairpin_special_H[GTTTTT] = -0.5;
	param_hairpin_special_S[GTTTTT] = 0;

	param_hairpin_special_H[TAAAAA] = 0.5;
	param_hairpin_special_S[TAAAAA] = 0.000322425;

	param_hairpin_special_H[TAAAAG] = 0.5;
	param_hairpin_special_S[TAAAAG] = -0.00161212;

	param_hairpin_special_H[TAAACA] = 0.7;
	param_hairpin_special_S[TAAACA] = 0.00161212;

	param_hairpin_special_H[TAAACG] = 1;
	param_hairpin_special_S[TAAACG] = 0.00161212;

	param_hairpin_special_H[TAACAA] = 1;
	param_hairpin_special_S[TAACAA] = 0.00161212;

	param_hairpin_special_H[TAACAG] = 1;
	param_hairpin_special_S[TAACAG] = 0.00161212;

	param_hairpin_special_H[TCTTGA] = 0;
	param_hairpin_special_S[TCTTGA] = 0.00419152;

	param_hairpin_special_H[TCTTGG] = 0;
	param_hairpin_special_S[TCTTGG] = 0.00322425;

	param_hairpin_special_H[TGAAA] = -1.5;
	param_hairpin_special_S[TGAAA] = 0;

	param_hairpin_special_H[TGAAAA] = -1.1;
	param_hairpin_special_S[TGAAAA] = 0.00161212;

	param_hairpin_special_H[TGAAAG] = -1;
	param_hairpin_special_S[TGAAAG] = 0.00161212;

	param_hairpin_special_H[TGAGAA] = -1.1;
	param_hairpin_special_S[TGAGAA] = 0.00161212;

	param_hairpin_special_H[TGAGAG] = -1;
	param_hairpin_special_S[TGAGAG] = 0.00161212;

	param_hairpin_special_H[TGATAA] = -1.6;
	param_hairpin_special_S[TGATAA] = 0.00161212;

	param_hairpin_special_H[TGATAG] = -1.5;
	param_hairpin_special_S[TGATAG] = 0.00161212;

	param_hairpin_special_H[TGCAA] = -1.5;
	param_hairpin_special_S[TGCAA] = 0;

	param_hairpin_special_H[TGCAAA] = -1.6;
	param_hairpin_special_S[TGCAAA] = 0.00161212;

	param_hairpin_special_H[TGCAAG] = -1.5;
	param_hairpin_special_S[TGCAAG] = 0.00161212;

	param_hairpin_special_H[TGCGAA] = -1.1;
	param_hairpin_special_S[TGCGAA] = 0.00161212;

	param_hairpin_special_H[TGCGAG] = -1;
	param_hairpin_special_S[TGCGAG] = 0.00161212;

	param_hairpin_special_H[TGCTTA] = 0.2;
	param_hairpin_special_S[TGCTTA] = 0.00161212;

	param_hairpin_special_H[TGCTTG] = -0.1;
	param_hairpin_special_S[TGCTTG] = 0.00161212;

	param_hairpin_special_H[TGGAA] = -1.5;
	param_hairpin_special_S[TGGAA] = 0;

	param_hairpin_special_H[TGGAAA] = -1.1;
	param_hairpin_special_S[TGGAAA] = 0.00161212;

	param_hairpin_special_H[TGGAAG] = -1;
	param_hairpin_special_S[TGGAAG] = 0.00161212;

	param_hairpin_special_H[TGGGAA] = -1.1;
	param_hairpin_special_S[TGGGAA] = 0.00161212;

	param_hairpin_special_H[TGGGAG] = -1;
	param_hairpin_special_S[TGGGAG] = 0.00161212;

	param_hairpin_special_H[TGGGGA] = 0.5;
	param_hairpin_special_S[TGGGGA] = 0.000644849;

	param_hairpin_special_H[TGGGGG] = 0.5;
	param_hairpin_special_S[TGGGGG] = 0.000644849;

	param_hairpin_special_H[TGTAA] = -1.5;
	param_hairpin_special_S[TGTAA] = 0;

	param_hairpin_special_H[TGTAAA] = -1.6;
	param_hairpin_special_S[TGTAAA] = 0.00161212;

	param_hairpin_special_H[TGTAAG] = -1.5;
	param_hairpin_special_S[TGTAAG] = 0.00161212;

	param_hairpin_special_H[TGTGAA] = -1.1;
	param_hairpin_special_S[TGTGAA] = 0.00161212;

	param_hairpin_special_H[TGTGAG] = -1;
	param_hairpin_special_S[TGTGAG] = 0.00161212;

	param_hairpin_special_H[TGTTCA] = 0.8;
	param_hairpin_special_S[TGTTCA] = 0.00161212;

	param_hairpin_special_H[TTTCGA] = -0.2;
	param_hairpin_special_S[TTTCGA] = 0.00161212;

	param_hairpin_special_H[TTTCGG] = -0.4;
	param_hairpin_special_S[TTTCGG] = 0.00161212;

	param_hairpin_special_H[TTTTAG] = -1;
	param_hairpin_special_S[TTTTAG] = 0.00161212;

	param_hairpin_special_H[TTTTGA] = 0;
	param_hairpin_special_S[TTTTGA] = 0.00161212;

	param_hairpin_special_H[TTTTGG] = -0.4;
	param_hairpin_special_S[TTTTGG] = 0.00161212;

	param_hairpin_special_H[TTTTTA] = -0.5;
	param_hairpin_special_S[TTTTTA] = 0.00161212;

	param_hairpin_special_H[TTTTTG] = -0.5;
	param_hairpin_special_S[TTTTTG] = 0.00161212;
}

