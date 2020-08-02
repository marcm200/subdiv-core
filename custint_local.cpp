/*

	simple interval struct with end-points double
	based on constant upward rounding, no
	rounding mode switch necessary during operations
	
	Marc Meidlinger
	May-July 2020
	
	NO WARRANTY OF CORRECT RESULTS (although I did my best to
	ensure correctness of the implementation)
	
	CAUTION: Every arithmetics routine cust^NNN_ZAB()
			 ASSUMES rounding mode is set to UPWARDS
			 
			 It is at the discretion of the user to
			 ensure this by defining a struct of
			 type CustRoundingUpwards
	
	based on the article:

	Interval arithmetic with fixed rounding modw
	by S. Rump et al., 2016

*/

#ifndef _CUSTOMIZED_INTERVAL
#define _CUSTOMIZED_INTERVAL

#include "custint.h"
#include "stdio.h"
#include "stdlib.h"
#include "stdint.h"
#include "math.h"
#include "fenv.h"
#include "float.h"

// get the current rounding mode, store it and
// set to upward infinity
struct CustRoundingUpwards {
	int32_t mode;
	
	CustRoundingUpwards();
	virtual ~CustRoundingUpwards();
};

struct CustRoundingDownwards {
	int32_t mode;
	
	CustRoundingDownwards();
	virtual ~CustRoundingDownwards();
};

// interval arithmetics with double and constant
// rounding upwards (downwards simulated by sign change)

struct CustInterval {
	double left,right;
	
	CustInterval();
	CustInterval(const double);
	CustInterval(const double,const double);
};

struct CustComplex {
	CustInterval re,im;
	
	CustComplex(const CustInterval&,const CustInterval&);
	CustComplex(const double,const double);
	CustComplex();
};

struct CustIntervalMatrix2x2 {
	CustInterval x1y1,x2y1,x1y2,x2y2;
};


// forward

// general functions
void terminate_program(void);
inline double minimumDouble(const double,const double,const double,const double);
inline double maximumDouble(const double,const double,const double,const double);

// CustInterval operations
int32_t custAdd_ZAB(CustInterval&,CustInterval&,CustInterval&);
int32_t custSub_ZAB(CustInterval&,CustInterval&,CustInterval&);
int32_t custMul_ZAB(CustInterval&,CustInterval&,CustInterval&);
int32_t custDiv_ZAB(CustInterval&,CustInterval&,CustInterval&);
int32_t custSqrt_ZA_switch(CustInterval&,CustInterval&);
int32_t custSqrt_ZA(CustInterval&,CustInterval&);
int32_t custPow2_ZA(CustInterval&,CustInterval&);
void custPrint(FILE*,CustInterval&);
int32_t custMidpoint_ZA(CustInterval&,CustInterval&);

// CustComplex operations
int32_t custCplxAdd_ZAB(CustComplex&,CustComplex&,CustComplex&);
int32_t custCplxSub_ZAB(CustComplex&,CustComplex&,CustComplex&);
int32_t custCplxMul_ZAB(CustComplex&,CustComplex&,CustComplex&);
int32_t custCplxMul_ZArealint64(CustComplex&,CustComplex&,const int64_t);
int32_t custCplxDiv_ZAB(CustComplex&,CustComplex&,CustComplex&);
int32_t custCplxNormQ_ZA(CustInterval&,CustComplex&);
int32_t custCplxSqrt_ZA(CustComplex&,CustComplex&);
int32_t custCplxMidpoint_ZA(CustComplex&,CustComplex&);

// custMatrix routines
int32_t custMatrixAdd_ZAB(CustIntervalMatrix2x2&,CustIntervalMatrix2x2&,CustIntervalMatrix2x2&);
int32_t custMatrixSub_ZAB(CustIntervalMatrix2x2&,CustIntervalMatrix2x2&,CustIntervalMatrix2x2&);
int32_t custMatrixMul_ZAB(CustIntervalMatrix2x2&,CustIntervalMatrix2x2&,CustIntervalMatrix2x2&);
int32_t custMatrixCplxMul_ZMA(CustComplex&,CustIntervalMatrix2x2&,CustComplex&);
int32_t custMatrixInverse_ZA(CustIntervalMatrix2x2&,CustIntervalMatrix2x2&);


// general function

void terminate_program(void) {
	fprintf(stderr,"Program terminates.\n");
	exit(99);
}

double minimumDouble(
	const double a,const double b,
	const double c,const double d
) {
	double m=a;
	if (b < m) m=b;
	if (c < m) m=c;
	if (d < m) m=d;
	
	return m;
}
	
double maximumDouble(
	const double a,const double b,
	const double c,const double d
) {
	double m=a;
	if (b > m) m=b;
	if (c > m) m=c;
	if (d > m) m=d;
	
	return m;
}


// struct CustRoundingUpwards 
// the constructor retrieves the current rounding mode and
// saves it, then sets to upward infinity
// the destructor resets to the old rounding mode

CustRoundingUpwards::CustRoundingUpwards() {
	mode=fegetround();
	if (mode<0) {
		fprintf(stderr,"CustRoundingUpwards::construct. Cannot get current rounding mode.\n");
		terminate_program();
	}
	if (fesetround(FE_UPWARD) != 0) {
		fprintf(stderr,"CustRoundingUpwards::constructor. Cannot set rounding mode to upwards.\n");
		terminate_program();
	}
	
	//printf("\nCustRoundingUpwards::set %i\n",fegetround());
}

CustRoundingUpwards::~CustRoundingUpwards() {
	if (fesetround(mode) != 0) {
		fprintf(stderr,"CustRoundingUpwards::destructor. Cannot reset rounding mode to %i.\n",mode);
		terminate_program();
	}
}

// CustRoundingDownards
CustRoundingDownwards::CustRoundingDownwards() {
	mode=fegetround();
	if (mode<0) {
		fprintf(stderr,"CustRoundingDownwards::construct. Cannot get current rounding mode.\n");
		terminate_program();
	}
	if (fesetround(FE_DOWNWARD) != 0) {
		fprintf(stderr,"CustRoundingDownwards::constructor. Cannot set rounding mode to upwards.\n");
		terminate_program();
	}
}

CustRoundingDownwards::~CustRoundingDownwards() {
	if (fesetround(mode) != 0) {
		fprintf(stderr,"CustRoundingDownwards::destructor. Cannot reset rounding mode to %i.\n",mode);
		terminate_program();
	}
}

// CustInterval

CustInterval::CustInterval() {
	// do nothing for speed reasons
}

CustInterval::CustInterval(const double a) {
	left=right=a;
	int32_t fpa=fpclassify(a);
	if (
		(fpa != FP_NORMAL) &&
		(fpa != FP_ZERO)
	) {
		fprintf(stderr,"CustInterval::constructor. Error, only normal floating-points (and 0) are allowed.\n");
		terminate_program();
	}
}

CustInterval::CustInterval(const double a,const double b) {
	left=a;
	right=b;
	int32_t fpa=fpclassify(a);
	int32_t fpb=fpclassify(b);
	if (
		(
			(fpa != FP_NORMAL) &&
			(fpa != FP_ZERO)
		) ||
		(
			(fpb != FP_NORMAL) &&
			(fpb != FP_ZERO)
		)
	) {
		fprintf(stderr,"CustInterval::constructor. Error, only normal floating-points (and 0) are allowed.\n");
		fprintf(stderr,": %.20lg..%.20lg\n",a,b);
		terminate_program();
	}
	
	if (a > b) {
		fprintf(stderr,"CustInterval::constructor. Not an interval. End points in wrong order.\n");
		terminate_program();
	}
}


// arithmetics

// all routines assume that the given parameters
// are normal or zero (i.e. no Inf, NaN or subnormals)
// (except division, that checks for the dividend
// not containing zero)
// resulting intervals are always normals or zero at
// end points (returning value 0).
// If not possible, the return value is -1
// and should be handled by the calling function appropriately

// a subnormal resulting interval end point is
// moved towards a circumferencing next normal (or zero)
#define IAVALIDATE_LEFT(VAR) \
{\
	switch (fpclassify(VAR)) {\
		case FP_INFINITE:\
		case FP_NAN: return -1; /* error */\
		case FP_SUBNORMAL: {\
			if (VAR > 0.0) {\
				/* left end, positive subnormal */\
				/* move to 0 */\
				VAR=0.0;\
			} else if (VAR < 0.0) {\
				/* left end, negative subnormal */\
				/* move to -FMIN */\
				VAR=-DBL_MIN;\
			}\
		}\
		default: {\
			break;\
		}\
	}\
}

#define IAVALIDATE_RIGHT(VAR) \
{\
	switch (fpclassify(VAR)) {\
		case FP_INFINITE:\
		case FP_NAN: return -1; /* error */\
		case FP_SUBNORMAL: {\
			if (VAR > 0.0) {\
				/* right end, positive subnormal */\
				/* move to smallest double */\
				VAR=DBL_MIN;\
			} else if (VAR < 0.0) {\
				/* right end, negative subnormal */\
				/* move to 0 */\
				VAR=0.0;\
			}\
		}\
		default: {\
			break;\
		}\
	}\
}

int32_t custAdd_ZAB(
	CustInterval& res,
	CustInterval& A,
	CustInterval& B
) {
	// /////////////////////////////////////////
	// rounding mode is assumed UPWARD
	// /////////////////////////////////////////

	// it is assumed that res is not the same object as A or B

	// left end: simulates downwards rounding by
	// consecutive negation (see articel by S.Rump)
	double a1=A.left; a1=-a1; // no rounding occurs
	double b1=B.left; b1=-b1; // no rounding occurs
	
	// ROUNDING upwards
	res.left=a1+b1; 
	
	res.left=-res.left; // no rounding
	// now the result is DOWNWARD-rounded A.left+B.left
	
	// right end-point: uses upward rounding
	res.right=A.right+B.right;
	
	// if result is subnormal => enlarge interval or
	// return error value -1 if infinity or NaN
	IAVALIDATE_LEFT(res.left);
	IAVALIDATE_RIGHT(res.right);
	
	// everything worked, valid interval in variable res
	return 0;
}
	
int32_t custSub_ZAB(
	CustInterval& res,
	CustInterval& A,
	CustInterval& B
) {
	// /////////////////////////////////////////
	// rounding mode is assumed UPWARD
	// /////////////////////////////////////////

	// it is assumed that res is not the same object as A or B

	// [aL..aR] - [bL..bR] = [aL-bR..aR-bL]
	
	// left end: simulates downwards rounding by
	// consecutive negation (see articel by S.Rump)
	// aL -DOWN bR = -( (-aL) UP- (-bR) )
	double a1=A.left;  a1=-a1; // no rounding occurs
	double b1=B.right; b1=-b1; // no rounding occurs
	
	// ROUNDING upwards
	res.left=a1-b1; 
	
	res.left=-res.left; // no rounding
	// now the result is DOWNWARD-rounded A.left-B.right
	
	// right end-point: uses upward rounding
	res.right=A.right-B.left;
	
	// if result is subnormal => enlarge interval or
	// return error value -1 if infinity or NaN
	IAVALIDATE_LEFT(res.left);
	IAVALIDATE_RIGHT(res.right);
	
	// everything worked, valid interval in variable res
	return 0;

}

int32_t custMul_ZAB(
	CustInterval& res,
	CustInterval& A,
	CustInterval& B
) {
	// /////////////////////////////////////////
	// rounding mode is assumed UPWARD
	// /////////////////////////////////////////

	// it is assumed that res is not the same object as A or B

	// could be implemented using sign checks
	
	// downward-mul: -( upward(a*(-b)) )
	
	// left end

	// forced execution order to 
	// obtain the desired ROUNDING effect
	double w1=B.left; w1=-w1; // no rounding
	
	// ROUNDING upward
	w1 *= A.left;
	
	// no rounding
	w1=-w1;
	// result is DOWNWARD-rounded A.left*B.left;
	
	// and the others
	double w2=B.right; w2=-w2; w2 *= A.left;  w2=-w2;
	double w3=B.left;  w3=-w3; w3 *= A.right; w3=-w3;
	double w4=B.right; w4=-w4; w4 *= A.right; w4=-w4;
	res.left=minimumDouble(w1,w2,w3,w4);
	
	// right end: only binary operation, so execution
	// order is clear and temporary results can be passed directly
	// to maximumDouble
	res.right=maximumDouble(
		A.left*B.left,
		A.left*B.right,
		A.right*B.left,
		A.right*B.right
	);
	
	// if resulting end point(s) are subnormal => adjust
	// if inf,nan => return error value -1
	IAVALIDATE_LEFT(res.left);
	IAVALIDATE_RIGHT(res.right);
	
	// everything worked, valid interval in variable res
	return 0;
}

int32_t custDiv_ZAB(
	CustInterval& res,
	CustInterval& A,
	CustInterval& B
) {
	// /////////////////////////////////////////
	// rounding mode is assumed UPWARD
	// /////////////////////////////////////////

	// it is assumed that res is not the same object as A or B

	// [a0..a1] / [b0..b1] = [a0..b0] * [1/b1..1/b0]
	// b must not contain zero
	
	// end-point test:
	// here also a check for inf/nan is possible without
	// speed-loss as the check for zero at the end-points
	// has to be done anyways
	if (fpclassify(B.left) != FP_NORMAL) return -1; // also checks for zero
	if (fpclassify(B.right)!= FP_NORMAL) return -1; // also checks for zero
	
	// zero within?
	if ( (B.left < 0.0) && (B.right > 0.0) ) return -1;
	
	CustInterval binv;
	// left end: downward 1/b.right = -( UPWARD( (-1)/B.right ) )
	binv.left=-1.0; // no rounding
	
	// ROUNDING upward
	binv.left /= B.right;
	
	binv.left = -binv.left; // no rounding
	// now binv.left is downward rounded 1/B.right
	
	// right end point: upward rounding
	binv.right=1.0/B.left;
	
	// subnormals can be moved to zero or DBL_min
	// i.e. binv can be a superset of the true 1/B
	// as basic interval operations are inclusion monotone,
	// multiplication below is a valid result
	IAVALIDATE_LEFT(binv.left)
	IAVALIDATE_RIGHT(binv.right)
	
	return custMul_ZAB(res,A,binv);
}

int32_t custSqrt_ZA_switch(
	CustInterval& res,
	CustInterval& A
) {
	// currently uses switching rounding mode
	//printf("CustIA::sqrt: experimental\n");
	
	{
		CustRoundingUpwards rdup;
		if (A.right < 0.0) return -1; // error
		res.right=sqrt(A.right);
		
		// implicit destruction
	}
	
	{
		CustRoundingDownwards rddown;
		if (A.left < 0.0) return -1; // error
		res.left=sqrt(A.left);
		
		// implicit destruction
	}
	
	// now rounding mode is set back to up
	if (fegetround() != FE_UPWARD) {
		fprintf(stderr,"custSqrt: error, not upward rounding\n");
		return -1;
	}

	IAVALIDATE_LEFT(res.left);
	IAVALIDATE_RIGHT(res.right);
	
	return 0;
}

int32_t custSqrt_ZA(
	CustInterval& res,
	CustInterval& A
) {
	fprintf(stderr,"CustIA::sqrt: noswitch experimental\n");
	return -1;
	
	// uses upward computed sqrt-function
	// and substracts (downard) 1 ulp
	// checks for underflow 
	
	
	return 0;
}

int32_t custPow2_ZA(CustInterval& erg,CustInterval& A) {
	// if straddling zero => result: 0..max(A.left^2,A.right^2)
	// else use multiplication routine
	
	int32_t error=0;
	
	if ( (A.left < 0.0) && (A.right > 0.0 ) ) {
		erg.left=0.0;
		
		// UPWARD rounding
		double l2=A.left*A.left;
		double r2=A.right*A.right;
		
		if (l2 > r2) erg.right=l2; 
		else erg.right=r2;
	} else {
		error += custMul_ZAB(erg,A,A);
	}
	
	return error;
}
	
// CustComplex

CustComplex::CustComplex(const CustInterval& v1,const CustInterval& v2) {
	re.left=v1.left;
	re.right=v1.right;
	im.left=v2.left;
	im.right=v2.right;
}

CustComplex::CustComplex(const double a,const double b) {
	// complex point
	re.left=re.right=a;
	im.left=im.right=b;
}

CustComplex::CustComplex() {
	// no initialization for speed up
}

int32_t custCplxAdd_ZAB(CustComplex& erg,CustComplex& A,CustComplex& B) {
	int32_t error=0;
	// real and imag interval are added separately
	error += custAdd_ZAB(erg.re,A.re,B.re);
	error += custAdd_ZAB(erg.im,A.im,B.im);
	
	return error;
}

int32_t custCplxSub_ZAB(CustComplex& erg,CustComplex& A,CustComplex& B) {
	int32_t error=0;
	// real and imag interval are subtracted separately
	error += custSub_ZAB(erg.re,A.re,B.re);
	error += custSub_ZAB(erg.im,A.im,B.im);
	
	return error;
}

int32_t custCplxMul_ZAB(CustComplex& erg,CustComplex& A,CustComplex& B) {
	int32_t error=0;
	// (Are+Aim*i) * (Bre+Bim*i) = 
	// (Are*Bre-Aim*Bim + (Are*Bim+Aim*Bre)*i)
	CustInterval arebre,arebim,aimbre,aimbim;
	
	error += custMul_ZAB(arebre,A.re,B.re);
	error += custMul_ZAB(arebim,A.re,B.im);
	error += custMul_ZAB(aimbre,A.im,B.re);
	error += custMul_ZAB(aimbim,A.im,B.im);
	
	error += custSub_ZAB(erg.re,arebre,aimbim);
	error += custAdd_ZAB(erg.im,arebim,aimbre);

	return error;
}

int32_t custCplxDiv_ZAB(CustComplex& erg,CustComplex& A,CustComplex& B) {
	int32_t error=0;
	// (a+b*i)/(c+d*i) =: (e+f*i)
	// by maxima:
	// realpart: (b*d+a*c)/(d^2+c^2)
	// imagpart: (b*c-a*d)/(d^2+c^2)
	CustInterval d2c2,d2,c2;
	error += custPow2_ZA(c2,B.re); // interval does not straddle 0, but may end-point contain it
	error += custPow2_ZA(d2,B.im);
	error += custAdd_ZAB(d2c2,d2,c2);

	// realpart
	CustInterval bd,ac,bdac;
	error += custMul_ZAB(bd,A.im,B.im);
	error += custMul_ZAB(ac,A.re,B.re);
	error += custAdd_ZAB(bdac,bd,ac);
	error += custDiv_ZAB(erg.re,bdac,d2c2);

	// imagpart: (b*c-a*d)/d2c2
	CustInterval bc,ad,bcad;
	error += custMul_ZAB(bc,A.im,B.re);
	error += custMul_ZAB(ad,A.re,B.im);
	error += custSub_ZAB(bcad,bc,ad);
	error += custDiv_ZAB(erg.im,bcad,d2c2);

	return error;
}

int32_t custCplxNormQ_ZA(CustInterval& normq,CustComplex& A) {
	// coputes the SQUARE of the norm [a..b]^2+[c..d]^2
	// from the input complex numer [a..b]x[c..d]
	int32_t error=0;
	
	CustInterval realq,imagq;
	// ATTN: re,im CAN contain zero, so
	// re*re is NOT applicable as it might result
	// in an norm-interval with negative endpoint
	// one has to use a power-definition
	error += custPow2_ZA(realq,A.re);
	error += custPow2_ZA(imagq,A.im);
	error += custAdd_ZAB(normq,realq,imagq);
	
	return error;
}

#define OUTIA(TEXT,WW) \
{\
	printf("%s; %.20lg..%.20lg\n",TEXT,WW.left,WW.right);\
}

int32_t custCplxSqrt_ZA(
	CustComplex& erg,
	CustComplex& val
) {
	// jump back with error as soon as possible
	
	if (
		(fpclassify(val.im.left) == FP_ZERO) &&
		(fpclassify(val.im.right) == FP_ZERO)
	) {
		// pure real root
		erg.im.left=erg.im.right=0.0;

		return custSqrt_ZA_switch(erg.re,val.re);
	}
	// using principal root formula from
	// wikipedia:
	
	// sqrt(a+b*i) =: e+f*i
	// e=sqrt( 0.5*( sqrt(a^2+b^2) + a) )
	// f=sign(b)*sqrt( 0.5*(sqrt(a^2+b^2) - a) )
	
	// so if val.im contains 0 at an end or within
	// return with error. Those have to be computed
	// using a different IA library
	// currently not applicable
	#define CUSTCONTAINSZERO(II) \
	(\
		(II.left <= 0.0) && (II.right >= 0.0)\
	)
	
	if (CUSTCONTAINSZERO(val.im) > 0) {
		return -1;
	}
	
	CustInterval a2,b2,a2b2,sqrta2b2;
	
	// direct return to reduce number of
	// operations, as sqrt is then compüted by
	// a different IA library, e.g. kv
	if (custMul_ZAB(a2,val.re,val.re) != 0) return -1;
	if (custMul_ZAB(b2,val.im,val.im) != 0) return -1;
	if (custAdd_ZAB(a2b2,a2,b2) != 0) return -1;
	if (custSqrt_ZA_switch(sqrta2b2,a2b2) != 0) return -1;
	CustInterval t1,t2;
	CustInterval cust05(0.5,0.5);
	CustInterval custminus1(-1.0,-1.0);
	if (custSub_ZAB(t1,sqrta2b2,val.re) != 0) return -1;
	if (custMul_ZAB(t2,t1,cust05) != 0) return -1;
	
	if (val.im.right < 0.0) {
		CustInterval t3;
		if (custSqrt_ZA_switch(t3,t2) != 0) return -1;
		if (custMul_ZAB(erg.im,t3,custminus1) != 0) return -1;;
	} else {
		if (custSqrt_ZA_switch(erg.im,t2) != 0) return -1;
	}
	
	// real part
	// e=sqrt( t5 )
	CustInterval t4,t5;
	if (custAdd_ZAB(t4,sqrta2b2,val.re) != 0) return -1;
	if (custMul_ZAB(t5,t4,cust05) != 0) return -1;
	if (custSqrt_ZA_switch(erg.re,t5) != 0) return -1;

	return 0;
}

void custPrint(FILE* f,CustInterval& a) {
	fprintf(f,"[%.20lg..%.20lg]",a.left,a.right);
}

/* custIntervalMatrix2x2 */

int32_t custMatrixAdd_ZAB(
	CustIntervalMatrix2x2& erg,
	CustIntervalMatrix2x2& A,
	CustIntervalMatrix2x2& B
) {
	int32_t error=0;
	
	error += custAdd_ZAB(erg.x1y1,A.x1y1,B.x1y1);
	error += custAdd_ZAB(erg.x2y1,A.x2y1,B.x2y1);
	error += custAdd_ZAB(erg.x1y2,A.x1y2,B.x1y2);
	error += custAdd_ZAB(erg.x2y2,A.x2y2,B.x2y2);
	
	return error;
}

int32_t custMatrixSub_ZAB(
	CustIntervalMatrix2x2& erg,
	CustIntervalMatrix2x2& A,
	CustIntervalMatrix2x2& B
) {
	int32_t error=0;
	
	error += custSub_ZAB(erg.x1y1,A.x1y1,B.x1y1);
	error += custSub_ZAB(erg.x2y1,A.x2y1,B.x2y1);
	error += custSub_ZAB(erg.x1y2,A.x1y2,B.x1y2);
	error += custSub_ZAB(erg.x2y2,A.x2y2,B.x2y2);
	
	return error;
}

int32_t custMatrixMul_ZAB(
	CustIntervalMatrix2x2& erg,
	CustIntervalMatrix2x2& A,
	CustIntervalMatrix2x2& B
) {
	int32_t error=0;
	
	/*
		(
			A11[a1..a2]	A21[b1..b2]
			A12[c1..c2]	A22[d1..d2]
		)
		
		*
		
		(
			B11[e1..e2]	B21[f1..f2]
			B12[g1..g2]	B22[h1..h2]
		)
		
		=
		
		(
			ae+bg		af+bh
			ce+dg		cf+dh
		)
		
		=
		
		(
			A11B11+A21B12		A11B21+A21B22
			A12B11+A22B12		A12B21+A22B22
		)
		
	*/

	CustInterval t1,t2;
	error += custMul_ZAB(t1,A.x1y1,B.x1y1);
	error += custMul_ZAB(t2,A.x2y1,B.x1y2);
	error += custAdd_ZAB(erg.x1y1,t1,t2);

	error += custMul_ZAB(t1,A.x1y1,B.x2y1);
	error += custMul_ZAB(t2,A.x2y1,B.x2y2);
	error += custAdd_ZAB(erg.x2y1,t1,t2);

	error += custMul_ZAB(t1,A.x1y2,B.x1y1);
	error += custMul_ZAB(t2,A.x2y2,B.x1y2);
	error += custAdd_ZAB(erg.x1y2,t1,t2);

	error += custMul_ZAB(t1,A.x1y2,B.x2y1);
	error += custMul_ZAB(t2,A.x2y2,B.x2y2);
	error += custAdd_ZAB(erg.x2y2,t1,t2);

	return error;
}

int32_t custMatrixCplxMul_ZMA(
	CustComplex& erg,
	CustIntervalMatrix2x2& M,
	CustComplex& A
) {
	// erg=M*A; matrix of CustInterval to the left
	
	/*

			M11	M21		Ar		M11*Ar+M21*Ai
					*		=
			M12	M22		Ai		M12*Ar+M22*Ai

	*/
	
	int32_t error=0;
	
	CustInterval t1,t2;
	error += custMul_ZAB(t1,M.x1y1,A.re);
	error += custMul_ZAB(t2,M.x2y1,A.im);
	error += custAdd_ZAB(erg.re,t1,t2);

	CustInterval t3,t4;
	error += custMul_ZAB(t3,M.x1y2,A.re);
	error += custMul_ZAB(t4,M.x2y2,A.im);
	error += custAdd_ZAB(erg.im,t3,t4);

			
	return error;
	
}

int32_t custMidpoint_ZA(CustInterval& erg,CustInterval& A) {
	// compue via IAmul
	CustInterval tx0,tx1,halb,t1;
	halb.left=halb.right=0.5;
	
	tx0.left=tx0.right=A.left;
	tx1.left=tx1.right=A.right;
	
	int32_t error=0;
	
	error += custAdd_ZAB(t1,tx0,tx1);
	error += custMul_ZAB(erg,halb,t1);
	
	// erg must be a subset (or equal) to A
	// so if rounding results in the "midpoint"
	// been shifted outside => error
	
	if (
		(erg.left < A.left) ||
		(erg.right > A.right)
	) error+=-1;
	
	return error;
}

int32_t custCplxMidpoint_ZA(
	CustComplex& erg,
	CustComplex& A
) {
	int32_t error=0;
	
	error += custMidpoint_ZA(erg.re,A.re);
	error += custMidpoint_ZA(erg.im,A.im);
	
	return error;
}

int32_t custMatrixInverse_ZA(
	CustIntervalMatrix2x2& erg,
	CustIntervalMatrix2x2& A
) {
	/*
					)-1
		x11		x21	)		1	x22		-x21
					)	=  ---
		x12		x22	)	   det	-x12	x11
					)
					
		det=x11*x22 - x21*x12
					
	*/
	
	int32_t error=0;
	
	CustInterval det,t1,t2;
	error += custMul_ZAB(t1,A.x1y1,A.x2y2);
	error += custMul_ZAB(t2,A.x2y1,A.x1y2);
	error += custSub_ZAB(det,t1,t2);
	
	// 0 in detrminant ?
	if (
		(fpclassify(det.left) == FP_ZERO) ||
		(fpclassify(det.right) == FP_ZERO) ||
		(
			(det.left <= 0.0) &&
			(det.right >= 0.0)
		)
	) {
		return -1; // not invertible
	}
	
	CustInterval t3,t4,minus1;
	minus1.left=minus1.right=-1.0;

	error += custDiv_ZAB(erg.x1y1,A.x2y2,det);
	error += custDiv_ZAB(t3,A.x2y1,det);
	error += custMul_ZAB(erg.x2y1,minus1,t3);
	error += custDiv_ZAB(t4,A.x1y2,det);
	error += custMul_ZAB(erg.x1y2,minus1,t4);
	error += custDiv_ZAB(erg.x2y2,A.x1y1,det);

	return error;
}

#endif

