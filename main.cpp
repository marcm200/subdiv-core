/*

	subdivision for root finding using
	the interval Newton operator and the custInt
	library
	
	based on:
	
	N. Kamath
	Subdivision Algorithms for Complex Root 
	Isolation: Empirical Comparisons

*/

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "stdint.h"
#include "string.h"

// custInt and function for polynomial of degree-10
// are included below


// consts

// 2^53 to output double precision numbers as rationals
const int32_t TWO53BITS=53;
const int64_t TWO53=(int64_t)1 << TWO53BITS;

// maximum degree (ANZCOEFF-1) of polynomial allocated
const int32_t ANZCOEFF=11;

// maximum number of boxes memory is allocated
// for in one list
const int32_t MAXBOXES=(int32_t)1 << 16;

// box types and colors
const uint8_t ROOTTYP_UNCLEAR=0;
const uint8_t ROOTTYP_NOROOT=1;
const uint8_t ROOTTYP_DEFINITEROOT=2;
const uint8_t COLORRED=3;
const uint8_t COLORBLUE=4;

// how deep the subdivision goes for a found root
const uint8_t ROOTSTOP_WHENFOUND=0;
const uint8_t ROOTSTOP_ATWIDTH=1;
const uint8_t ROOTSTOP_DEEPEST=2;

// custInt library
#include "custint_local.cpp"


// structs

struct CoeffParam {
	CustComplex Ai[ANZCOEFF];
};

struct RootBox {
	uint8_t roottyp;
	CustComplex rect;
};

typedef unsigned char BYTE;

// palette entry
struct RGB3 {
	BYTE R,G,B;
};

// bitmap-like image allocated in one huge memory block
struct Charmap {
	int32_t xlen,ylen; 
	BYTE* dataYX;
	RGB3 palette[256];
		
	Charmap();
	virtual ~Charmap();
		
	void setlenxy(const int32_t,const int32_t);
	void save(const char*);
	void lineHorizVert(const int32_t,const int32_t,const int32_t,const int32_t,const BYTE);
	void fillrect(const int32_t,const int32_t,const int32_t,const int32_t,const BYTE);
	void setPaletteRGB(const int32_t,const BYTE,const BYTE,const BYTE);
	void setPunkt(const int32_t,const int32_t,const BYTE);
	BYTE getPunkt(const int32_t,const int32_t);
};

// parameters for the subdivision routine
struct SubDivParam {
	// nbr of boxes maximally storable in
	// variables tmpboxesN and returnedboxes
	int32_t allocatedboxes;
	
	int32_t ctrreturnedboxes;
	// memory has to be set externally and deallocated
	RootBox *returnedboxes;
	
	// complex coefficients of the polynomial
	CoeffParam coeff;

	// temporary, has to be set and deallocated
	// externally
	RootBox *tmpboxes1;
	RootBox *tmpboxes2;

	// how deep to go
	int32_t maxrefinementdepth;
	
	// start box - may encompass all roots
	// or just a region of interest
	RootBox startbox;
	
	// function pointers
	// compute image of an input box
	int32_t (*fA_ZA)(CoeffParam&,CustComplex&,CustComplex&);
	// partial derivatives
	int32_t (*ux_ZA)(CoeffParam&,CustInterval&,CustComplex&);
	int32_t (*uy_ZA)(CoeffParam&,CustInterval&,CustComplex&);
	// Cauchyidentity: ux=vy, uy=-vx
	
	// image pointer, memory for the image has to already
	// be allocated
	Charmap *img;
	
	virtual ~SubDivParam();

};


// defines as small functions

#define SETCOEFFCPLX2(VV,NR,RR0,RR1,II0,II1) \
{\
	VV.Ai[NR].re.left=RR0;\
	VV.Ai[NR].re.right=RR1;\
	VV.Ai[NR].im.left=II0;\
	VV.Ai[NR].im.right=II1;\
}

#define CLEARCOEFF(VV) \
{\
	for(int32_t i=0;i<ANZCOEFF;i++) {\
		VV.Ai[i].re.left=\
		VV.Ai[i].re.right=\
		VV.Ai[i].im.left=\
		VV.Ai[i].im.right=0.0;\
	}\
}

#define LOGMSG(TT) \
{\
	fprintf(flog,TT); fflush(flog);\
	printf(TT);\
}

#define LOGMSG2(TT,AA) \
{\
	fprintf(flog,TT,AA); fflush(flog);\
	printf(TT,AA);\
}

#define LOGMSG3(TT,AA,BB) \
{\
	fprintf(flog,TT,AA,BB); fflush(flog);\
	printf(TT,AA,BB);\
}

#define LOGMSG4(TT,AA,BB,CC) \
{\
	fprintf(flog,TT,AA,BB,CC); fflush(flog);\
	printf(TT,AA,BB,CC);\
}

#define LOGMSG5(TT,AA,BB,CC,DD) \
{\
	fprintf(flog,TT,AA,BB,CC,DD); fflush(flog);\
	printf(TT,AA,BB,CC,DD);\
}

#define ADDBOX(PTR,CTR,TT,BB) \
{\
	if (CTR > (subject.allocatedboxes-8)) {\
		/* no more memory */\
		globalerror=-1;\
	} else {\
		PTR[CTR].roottyp=TT;\
		PTR[CTR].rect.re.left=BB.re.left;\
		PTR[CTR].rect.re.right=BB.re.right;\
		PTR[CTR].rect.im.left=BB.im.left;\
		PTR[CTR].rect.im.right=BB.im.right;\
		CTR++;\
	}\
}
	
#define ADDBOX2(PTR,CTR,TT,XX0,XX1,YY0,YY1) \
{\
	if (CTR > (subject.allocatedboxes-8)) {\
		/* no more memory */\
		globalerror=-1;\
	} else {\
		PTR[CTR].roottyp=TT;\
		PTR[CTR].rect.re.left=XX0;\
		PTR[CTR].rect.re.right=XX1;\
		PTR[CTR].rect.im.left=YY0;\
		PTR[CTR].rect.im.right=YY1;\
		CTR++;\
	}\
}

#define AUSGABECPLX(CC) \
{\
	LOGMSG5("approx [%.20lg..%.20lg] x [%.20lg..%.20lg]",\
		CC.re.left,\
		CC.re.right,\
		CC.im.left,\
		CC.im.right\
	);\
}

#define GET64FROMDOUBLE(ERG,VAL) \
{\
	double sig;\
	int32_t exponent;\
	sig=frexp(VAL,&exponent);\
	double sig64d=TWO53*sig;\
	if (floor(sig64d) == sig64d) {\
		sprintf(ERG,"%.20lg*2^%i",sig64d,exponent-TWO53BITS);\
	} else {\
		sprintf(ERG,"(?)");\
	}\
}

#define AUSGABECPLX64(CC) \
{\
	LOGMSG(" exact [");\
	char tt[1024];\
	GET64FROMDOUBLE(tt,CC.re.left);\
	LOGMSG2("%s..",tt);\
	GET64FROMDOUBLE(tt,CC.re.right);\
	LOGMSG2("%s]+i*[",tt);\
	GET64FROMDOUBLE(tt,CC.im.left);\
	LOGMSG2("%s..",tt);\
	GET64FROMDOUBLE(tt,CC.im.right);\
	LOGMSG2("%s])",tt);\
}

// globals

FILE *flog=NULL;
// overview subdivision image. Must be divisible by 4
// to save a viewable image
int32_t SCREENWIDTH=(int32_t)1 << 10; 

// root box parameter
double stopwidth=1.0;
uint8_t STOPATROOT=ROOTSTOP_WHENFOUND;

// polynomial of degree 10
#include "_func_poly10_local.cpp"


// routines

char* chomp(char* s) {
	if (!s) return 0;
	for(int i=strlen(s);i>=0;i--) if (s[i]<32) s[i]=0; else break;
	return s;
}

char* upper(char* s) {
	if (!s) return NULL;
	
	for(int32_t i=0;i<(int32_t)strlen(s);i++) {
		if ((s[i]>='a')&&(s[i]<='z')) s[i]=s[i]-'a'+'A';
	}

	return s;
}

void setPaletteTo(Charmap& md) {
	// image colors
	for(int32_t i=0;i<256;i++) md.setPaletteRGB(i,255,0,0);
	md.setPaletteRGB(ROOTTYP_UNCLEAR,127,127,127);
	md.setPaletteRGB(ROOTTYP_DEFINITEROOT,0,0,0);
	md.setPaletteRGB(ROOTTYP_NOROOT,255,255,255);
	md.setPaletteRGB(COLORRED,255,0,0);
	md.setPaletteRGB(COLORBLUE,0,0,255);
}

// struct Charmap

void write2(FILE *f,const BYTE a,const BYTE b) {
	fwrite(&a,1,sizeof(a),f);
	fwrite(&b,1,sizeof(b),f);
}

void write4(FILE *f,const BYTE a,const BYTE b,const BYTE c,const BYTE d) {
	fwrite(&a,1,sizeof(a),f);
	fwrite(&b,1,sizeof(b),f);
	fwrite(&c,1,sizeof(c),f);
	fwrite(&d,1,sizeof(d),f);
}

void Charmap::lineHorizVert(const int32_t ax,const int32_t ay,const int32_t bx,const int32_t by,const BYTE awert) {
	if (!dataYX) return;
	
	if (ax == bx) {
		// vertical line
		int32_t y0,y1;
		if (ay < by) { y0=ay; y1=by; } else { y0=by; y1=ay; }
		for(int32_t y=y0;y<=y1;y++) {
			dataYX[y*xlen+ax]=awert;
		} // y
	} else if (ay == by) {
		// horizontal line
		int32_t x0,x1;
		if (ax < bx) { x0=ax; x1=bx; } else { x0=bx; x1=ax; }
		for(int32_t x=x0;x<=x1;x++) {
			dataYX[ay*xlen+x]=awert;
		} // y
	} else {
		// geht nicht
		printf("Charmap: Diagonal line not implemented.\n");
		return;
	}
}

void Charmap::save(const char* afn) {
	// data is stored in the format of an 8 bit Bitmap
	
	FILE *fbmp=fopen(afn,"wb");
	write2(fbmp,66,77); // BM
	
	uint32_t off
		=		14 // FILEHeader
			+	40 // Bitmapheader
			+	256*4 // ColorPalette
		;
	
	// filelen will overflow if image width is too large
	// but external viewers can display the image nonetheless
	uint32_t filelen
			=	off
			+	(ylen*xlen);
		;
			
	fwrite(&filelen,1,sizeof(filelen),fbmp);
	write4(fbmp,0,0,0,0);
	fwrite(&off,1,sizeof(off),fbmp); // offset, ab da beginnen die PIXEL
	write4(fbmp,40,0,0,0);
	
	uint32_t w = xlen;
	fwrite(&w,sizeof(w),1,fbmp);
	w = ylen;
	fwrite(&w,sizeof(w),1,fbmp);
	write2(fbmp,1,0); 
	write2(fbmp,8,0); // 8 bits per pixel
	write4(fbmp,0,0,0,0);
	write4(fbmp,0,0,0,0);
	write4(fbmp,19,10,0,0);
	write4(fbmp,19,10,0,0);
	write4(fbmp,0,1,0,0); // number of palette entries
	write4(fbmp,0,0,0,0);
	BYTE puffer[4];
	for(int32_t i=0;i<256;i++) {
		puffer[0]=palette[i].B;
		puffer[1]=palette[i].G;
		puffer[2]=palette[i].R;
		puffer[3]=0;
		fwrite(puffer,4,sizeof(BYTE),fbmp);
	}
	
	fwrite(&dataYX[0],xlen*ylen,sizeof(BYTE),fbmp);
	
	fclose(fbmp);	
}

void Charmap::fillrect(
	const int32_t ax,const int32_t ay,
	const int32_t bx,const int32_t by,
	const BYTE aw
) {
	int32_t lx=ax,rx=bx;
	if (ax > bx) { lx=bx; rx=ax; }
	int32_t ly=ay,ry=by;
	if (ay > by) { ly=by; ry=ay; }
	
	for(int32_t y=ly;y<=ry;y++) {
		for(int32_t x=lx;x<=rx;x++) {
			dataYX[y*xlen+x]=aw;
		}
	}
}

void Charmap::setlenxy(const int32_t ax,const int32_t ay) {
	// already allocated: so take that memory
	if ( 
		(ax == xlen) && 
		(ay == ylen) &&
		(dataYX)
	) {
		// same size => use already allocated memory
		return;
	}
	
	if ((xlen>0)&&(dataYX)) {
		delete[] dataYX;
	}
	
	xlen=ax;
	ylen=ay;
	dataYX=new BYTE[ay*ax];
	
}

Charmap::Charmap() {
	xlen=ylen=0;
	dataYX=NULL;
}

Charmap::~Charmap() {
	if ((xlen>0)&&(dataYX)) {
		delete[] dataYX;
	}
}

void Charmap::setPaletteRGB(const int32_t pos,const BYTE ar,const BYTE ag,const BYTE ab) {
	if ((pos<0)||(pos>255)) return;
	palette[pos].R=ar;
	palette[pos].G=ag;
	palette[pos].B=ab;
}

void Charmap::setPunkt(const int32_t ax,const int32_t ay,const BYTE awert) {
	if (!dataYX) return;
	if (
		(ax<0) || (ax >= xlen) || (ay<0) || (ay>=ylen)
	) return;
	
	dataYX[ay*xlen+ax]=awert;
}

BYTE Charmap::getPunkt(const int32_t ax,const int32_t ay) {
	if (!dataYX) {
		LOGMSG("Error. Image not set\n");
		exit(99);
	}
	if (
		(ax<0) || (ax >= xlen) || (ay<0) || (ay>=ylen)
	) {
		LOGMSG3("Error. Out of screen %i,%i\n",ax,ay);
		exit(99);
	}

	return dataYX[ay*xlen+ax];
}

// polynom10

int32_t compute_polynom10_fA_ZA(
	CoeffParam& coeff,
	CustComplex& fA,
	CustComplex& A
) {
	int32_t error=0;
	
	// real-part component function
	error += get_poly10_u(coeff,A,fA.re);
	// imaginary-part component function
	error += get_poly10_v(coeff,A,fA.im);
	
	return error;
}

int32_t compute_partial_polynom10_ux_ZA(
	CoeffParam& coeff,
	CustInterval& fA,
	CustComplex& A
) {
	return get_poly10_u_dx(coeff,A,fA);
}

int32_t compute_partial_polynom10_uy_ZA(
	CoeffParam& coeff,
	CustInterval& fA,
	CustComplex& A
) {
	return get_poly10_u_dy(coeff,A,fA);
}

int32_t iacontainszero(CustInterval& a) {
	if (
		(a.left <= 0.0) &&
		(a.right >= 0.0)
	) return 1;
	
	return 0;
}

int32_t IAoperator_Newton_ZAS(
	CustComplex& NA,
	CustComplex& A,
	SubDivParam& subject
) {

	// Newton operator: 
	// N(A) = M - iajac(A)^-1 * f(M)
	// pp49 of article in preface
	
	int32_t error=0;
	
	CustComplex midp,fmidp;
	error += custCplxMidpoint_ZA(midp,A);
	error += subject.fA_ZA(subject.coeff,fmidp,midp);
	
	CustInterval eux,euy,evx,evy;
	error += subject.ux_ZA(subject.coeff,eux,A);
	error += subject.uy_ZA(subject.coeff,euy,A);
	// ux=vy
	evy.left=eux.left;
	evy.right=eux.right;
	// -uy=vx
	CustInterval minus1;
	minus1.left=minus1.right=-1.0;
	error += custMul_ZAB(evx,euy,minus1);
	
	// Jacobian:	ux	uy
	//				vx	vy
	CustIntervalMatrix2x2 jac;
	jac.x1y1.left=eux.left;
	jac.x1y1.right=eux.right;
	jac.x2y1.left=euy.left;
	jac.x2y1.right=euy.right;
	jac.x1y2.left=evx.left;
	jac.x1y2.right=evx.right;
	jac.x2y2.left=evy.left;
	jac.x2y2.right=evy.right;
	
	CustIntervalMatrix2x2 invjac;
	
	error += custMatrixInverse_ZA(invjac,jac);
	// not invertible or general error here or before
	if (error != 0) return -1; 

	// N(A) = M - iajac(A)^-1 * f(M)
	CustComplex e1;
	error += custMatrixCplxMul_ZMA(e1,invjac,fmidp);
	// N(A) = M - e1
	error += custCplxSub_ZAB(NA,midp,e1);
	
	return error;
}

int32_t custIA_subdivision(SubDivParam& subject) {
	/* 

		return-value: -1 < memory error occured, DEFINITE boxes still valid 
		0: no error, but also no definite root regions
		1: no error and definite root regions

	*/
	
	int32_t globalerror=0;
	
	double zx0=1,zx1=1,zy0=1,zy1=1,skx=1,sky=1;
	if (subject.img) {
		subject.img->fillrect(0,0,subject.img->xlen-1,subject.img->ylen-1,ROOTTYP_UNCLEAR);
		zx0=subject.startbox.rect.re.left;
		zx1=subject.startbox.rect.re.right;
		zy0=subject.startbox.rect.im.left;
		zy1=subject.startbox.rect.im.right;
		skx = (zx1-zx0) / subject.img->xlen;
		sky = (zy1-zy0) / subject.img->ylen;
	}
	
	// if N(A) ntersect A = empty => no root in A, discard
	// if N(A) subset or equal to A => root 
	
	subject.ctrreturnedboxes=0;
	int32_t ctrin=0,ctrout=0;
	RootBox *pin=subject.tmpboxes1,*pout=subject.tmpboxes2;
	
	ADDBOX(pin,ctrin,ROOTTYP_UNCLEAR,subject.startbox.rect);
	
	#define TRIMX(XX) \
	{\
		if (XX < 0) XX=0; else if (XX >= subject.img->xlen) XX=subject.img->xlen-1;\
	}
	
	#define TRIMY(YY) \
	{\
		if (YY < 0) YY=0; else if (YY >= subject.img->ylen) YY=subject.img->ylen-1;\
	}

	#define DRAWRECTANGLE(CC,FF,BORDER) \
	{\
		if (\
			(CC.re.right < subject.startbox.rect.re.left) ||\
			(CC.re.left > subject.startbox.rect.re.right) ||\
			(CC.im.right < subject.startbox.rect.im.left) ||\
			(CC.im.left > subject.startbox.rect.im.right) \
		) {\
			LOGMSG("\nImplementation error. A fully outside\n");\
			exit(99);\
		}\
		/* just an indicator, so rounding not relevant currently */\
		int32_t xx0=(int32_t)floor( (CC.re.left - zx0) / skx);\
		int32_t xx1=(int32_t)floor( (CC.re.right - zx0) / skx);\
		int32_t yy0=(int32_t)floor( (CC.im.left - zy0) / sky);\
		int32_t yy1=(int32_t)floor( (CC.im.right - zy0) / sky);\
		/* rectangle drawable */\
		TRIMX(xx0)\
		TRIMX(xx1)\
		TRIMY(yy0)\
		TRIMY(yy1)\
		subject.img->fillrect(xx0,yy0,xx1,yy1,FF);\
		if (\
			( (xx1-xx0) >= 3 ) &&\
			( (yy1-yy0) >= 3 ) &&\
			(BORDER > 0)\
		) {\
			subject.img->lineHorizVert(xx0,yy0,xx0,yy1,COLORRED);\
			subject.img->lineHorizVert(xx0,yy1,xx1,yy1,COLORRED);\
			subject.img->lineHorizVert(xx1,yy1,xx1,yy0,COLORRED);\
			subject.img->lineHorizVert(xx1,yy0,xx0,yy0,COLORRED);\
		}\
	}
	
	int32_t def=0;

	for(int32_t level=1;level <= subject.maxrefinementdepth;level++) {
		
		if (
			(ctrin >= (MAXBOXES-16)) ||
			(globalerror != 0)
		) {
			// too many => nothing to copy
			// (arbitrary choice to speed up)
			
			// returnedboxes may already contain some
			// definite roots
			
			ctrout=0;
			ctrin=0;
			break;
		}
		ctrout=0;
		for(int32_t inc=0;inc<ctrin;inc++) {
			RootBox *A=&pin[inc];
			
			CustComplex Asub;
			Asub.re.left=A->rect.re.left;
			Asub.re.right=A->rect.re.right;
			Asub.im.left=A->rect.im.left;
			Asub.im.right=A->rect.im.right;
			
			// some baxes can be discarded before computing 
			// the operator 
			
			CustComplex f;
			if (subject.fA_ZA(subject.coeff,f,A->rect) == 0) {
				// intersection with origin ?

				if (
					(iacontainszero(f.re) <= 0) ||
					(iacontainszero(f.im) <= 0)
				) {
					// no => discard A
					if (subject.img) {
						DRAWRECTANGLE(A->rect,ROOTTYP_NOROOT,1)
					} // draw white rectangle
					
					continue;
				}
			} // f

			// analyze boxes
			CustComplex KA; // operator result
			
			if (IAoperator_Newton_ZAS(KA,A->rect,subject) == 0) {
				if (
					(KA.re.left > A->rect.re.right) ||
					(KA.re.right < A->rect.re.left) ||
					(KA.im.left > A->rect.im.right) ||
					(KA.im.right < A->rect.im.left)
				) {
					// no intersection
					// A can be discarded
					if (subject.img) {
						DRAWRECTANGLE(A->rect,ROOTTYP_NOROOT,1)
					} // draw white rectangle
				
					continue;
				} // no intersection A,KA
			
				// does A contain a root
				if (
					(A->rect.re.left <= KA.re.left) &&
					(KA.re.right <= A->rect.re.right) &&
					(A->rect.im.left <= KA.im.left) &&
					(KA.im.right <= A->rect.im.right)
				) {
					int8_t subdivide=1;
					
					if (
						(level >= 2)
					) {
						if (STOPATROOT == ROOTSTOP_WHENFOUND) subdivide=0;
						else 
						if (STOPATROOT == ROOTSTOP_ATWIDTH) {
							// the computation of width here is
							// NOT rounding-controlled as it
							// does not affect the validity of
							// the result, just the level at which
							// it is reported
							
							if (
								((A->rect.re.right - A->rect.re.left) <= stopwidth) ||
								((A->rect.im.right - A->rect.im.left) <= stopwidth)
							) subdivide=0;
						}
						
					} // level
					
					if (level == subject.maxrefinementdepth) subdivide=0;

					if (subdivide <= 0) {
						ADDBOX(
							subject.returnedboxes,
							subject.ctrreturnedboxes,
							ROOTTYP_DEFINITEROOT,
							A->rect
						);
						if (
							(subject.img)
						) {
							DRAWRECTANGLE(A->rect,ROOTTYP_DEFINITEROOT,0)
						}
						continue;
					} // no subdivision
						
				} // containing a root 
				
			} // IAop returned valid value
			
			// subdivide Asub into 2x2 squares
			// taking the LEFT of the midpoint interval for every
			// position, so the subdivided parts
			// share an edge even in case of rounded
			// numbers
				
			CustComplex midpsub;
			if (custCplxMidpoint_ZA(midpsub,Asub) != 0) {
				// error cannot subdivide
				LOGMSG("\nError. Midpoint not consistent.\n");
				exit(99);
			}
			
			// midp can be a small interval
			
			double mx=midpsub.re.left;
			double my=midpsub.im.left;
				
			// lower left
			ADDBOX2(pout,ctrout,ROOTTYP_UNCLEAR,
				Asub.re.left,mx,
				Asub.im.left,my
			);
			
			// lower right
			ADDBOX2(pout,ctrout,ROOTTYP_UNCLEAR,
				mx,Asub.re.right,
				Asub.im.left,my
			);
			
			// upper left
			ADDBOX2(pout,ctrout,ROOTTYP_UNCLEAR,
				Asub.re.left,mx,
				my,Asub.im.right
			);
			
			// upper right
			ADDBOX2(pout,ctrout,ROOTTYP_UNCLEAR,
				mx,Asub.re.right,
				my,Asub.im.right
			);
			
		} // while
		
		// swap lists
		RootBox *p=pin;
		pin=pout;
		pout=p;
		ctrin=ctrout;
		ctrout=0;
		
		// are there any left to analyze ?
		if (ctrin <= 0) {
			break; // nothing to do
		}
		
	} // for level
	
	// now copy, what's left in pin into returnedboxes

	// two sweeps: one for copying and drawing
	// gray UNCLEAR regions

	// there can already be definite root boxes
	// in list returnedboxes
	
	for(int32_t i=0;i<ctrin;i++) {
		if (pin[i].roottyp == ROOTTYP_DEFINITEROOT) def=1;
		
		if (
			(subject.img) &&
			(
				(pin[i].roottyp == ROOTTYP_UNCLEAR)
			)
		) {
			DRAWRECTANGLE(pin[i].rect,ROOTTYP_UNCLEAR,0);
		} 
		
		ADDBOX(
			subject.returnedboxes,
			subject.ctrreturnedboxes,
			pin[i].roottyp,
			pin[i].rect
		)
	} // i
	
	// 2nd pass: draw BLACK pixels so they take
	// precedence
	if (subject.img) {
		for(int32_t i=0;i<subject.ctrreturnedboxes;i++) {
			if (
				(subject.returnedboxes[i].roottyp == ROOTTYP_DEFINITEROOT)
			) {
				DRAWRECTANGLE(
					subject.returnedboxes[i].rect,
					subject.returnedboxes[i].roottyp,0);
			}
		} //
	} // if

	if (globalerror != 0) return globalerror;
	
	return def;
}

void ausgabePolynom(CoeffParam& coeff) {
	for(int32_t i=(ANZCOEFF-1);i>=0;i--) {
		if (
			(coeff.Ai[i].re.left != 0.0) ||
			(coeff.Ai[i].re.right != 0.0) ||
			(coeff.Ai[i].im.left != 0.0) ||
			(coeff.Ai[i].im.right != 0.0)
		) {
			LOGMSG5("([%.20lg..%.20lg]+i*[%.20lg..%.20lg])",
				coeff.Ai[i].re.left,
				coeff.Ai[i].re.right,
				coeff.Ai[i].im.left,
				coeff.Ai[i].im.right
			);
		
			if (i>0) {
				LOGMSG("*z");
				if (i > 1) {
					LOGMSG2("^%i",i);
				}
				LOGMSG(" + ");
			}
		} // coefficient non-zero
	}
}

void ausgabeRoots(SubDivParam& subject) {
	LOGMSG("\nundecided regions\n");
	int8_t da=0;
	for(int32_t i=0;i<subject.ctrreturnedboxes;i++) {
		if (subject.returnedboxes[i].roottyp != ROOTTYP_UNCLEAR) continue;
		
		da=1;
		LOGMSG2("\n#%i ",i);
		LOGMSG("undecided\n  ");
		AUSGABECPLX(subject.returnedboxes[i].rect);
		LOGMSG("\n");
	}
	if (da <= 0) {
		LOGMSG("  none\n");
	}
	
	LOGMSG("\n\nDEFINITE root containing regions");
	da=0;
	for(int32_t i=0;i<subject.ctrreturnedboxes;i++) {
		if (subject.returnedboxes[i].roottyp != ROOTTYP_DEFINITEROOT) continue;
		
		da=1;
		LOGMSG2("\n\n#%i DEFINITE:\n  ",i);
		AUSGABECPLX64(subject.returnedboxes[i].rect);
		LOGMSG("\n  ");
		AUSGABECPLX(subject.returnedboxes[i].rect);
		CustInterval no;
		if (custCplxNormQ_ZA(no,subject.returnedboxes[i].rect) != 0) {
			LOGMSG(" Error norm.\n");
			exit(99);
		}

		LOGMSG("\n");
		if (no.right < 1.0) {
			LOGMSG("  *** fully in unit circle *** ");
		}

		LOGMSG("  ||squared=");
		char tt[2048];
		GET64FROMDOUBLE(tt,no.left)
		LOGMSG2("[%s..",tt);
		GET64FROMDOUBLE(tt,no.right)
		LOGMSG4("..%s] (approx [%.5lg..%.5lg]",tt,no.left,no.right);

	} // for
	if (da <= 0) {
		LOGMSG("\n  none\n");
	}
}

void find_roots(SubDivParam& subject) {
	// coefficients have already been set

	// polynom 10
	subject.ux_ZA=compute_partial_polynom10_ux_ZA;
	subject.uy_ZA=compute_partial_polynom10_uy_ZA;
	subject.fA_ZA=compute_polynom10_fA_ZA;
	
	subject.allocatedboxes=MAXBOXES;
	subject.tmpboxes1=new RootBox[subject.allocatedboxes];
	subject.tmpboxes2=new RootBox[subject.allocatedboxes];
	subject.returnedboxes=new RootBox[subject.allocatedboxes];
	
	// if image is not wanted, set point to NULL instead
	subject.img=new Charmap;
	subject.img->setlenxy(SCREENWIDTH,SCREENWIDTH);
	setPaletteTo(*subject.img);
	
	LOGMSG("\nATTN: coefficient output is approximate\n");
	LOGMSG("\n  f(z)= ");
	ausgabePolynom(subject.coeff);
	LOGMSG("\n");
	
	int32_t ret = custIA_subdivision(subject);

	if (ret < 0) {
		ausgabeRoots(subject);
		LOGMSG("\n\nMethod returned box list memory error. Only DEFINITEROOTS are valid.\n");
		// returned boxes are still valid if they
		// have flag DEFINITEROOT
		// however the union of all returned boxes does not
		// necessarily cover all roots
		exit(99);
	} else {
		LOGMSG2("\n\n%i boxes returned in total\n",subject.ctrreturnedboxes);
		ausgabeRoots(subject);
	}

	if (subject.img) {
		// draw unit-circle
		double skx=(subject.startbox.rect.re.right-subject.startbox.rect.re.left)/subject.img->xlen;
		double sky=(subject.startbox.rect.im.right-subject.startbox.rect.im.left)/subject.img->ylen;
		for(double d=0.0;d<=(2*M_PI);d+=0.001) {
			int32_t xx=(int32_t)floor( ( cos(d) - subject.startbox.rect.re.left) / skx );
			int32_t yy=(int32_t)floor( ( sin(d) - subject.startbox.rect.im.left) / sky );
			if (
				(xx>=0) && (xx<subject.img->xlen) &&
				(yy>=0) && (yy<subject.img->ylen)
			) {
				subject.img->setPunkt(xx,yy,COLORBLUE);
			}
		} // d
		
		subject.img->save("_roots.bmp");
	}
}

// struct SubDivParam

SubDivParam::~SubDivParam() {
	if (img) delete img;
	if (tmpboxes1) delete[] tmpboxes1;
	if (tmpboxes2) delete[] tmpboxes2;
}

int32_t loadCoefficients(
	const char* afn,
	CoeffParam& coeff
) {
	FILE *f=fopen(afn,"rt");
	if (!f) return 0;
	
	char tmp[1024];
	while (!feof(f)) {
		fgets(tmp,1000,f);
		chomp(tmp); upper(tmp);
		
		if (tmp[0] == '/') continue; // comment
		if (tmp[0] == '.') break; // end of file
		
		if (tmp[0] == 'R') {
			int32_t i;
			double a,b,c,d;
			if (sscanf(&tmp[1],"%i,%lf,%lf,%lf,%lf",&i,&a,&b,&c,&d) == 5) {
				if (
					(i < 0) ||
					(i >= ANZCOEFF)
				) return -1; // not supported
				
				SETCOEFFCPLX2(coeff,i,a,b,c,d);
			} else
			if (sscanf(&tmp[1],"%i,%lf,%lf",&i,&a,&b) == 3) {
				if (
					(i < 0) ||
					(i >= ANZCOEFF)
				) return -1; // not supported
				
				SETCOEFFCPLX2(coeff,i,a,a,b,b);
			} else {
				LOGMSG2("\nError. Coefficient file: Unknown line |%s|\n",
					tmp);
				exit(99);
			}
		} else
		if (tmp[0] == 'D') {
			// line: Dindex,denominator,a,b,c,d
			// values are then a/denominator etc...
			// it is at the discretion of the user to
			// see that those values are accurately 
			// representable in C++ double (if accuracy is needed)
			int32_t i;
			int64_t denom,a,b,c,d;

			if (sscanf(&tmp[1],"%i,%I64d,%I64d,%I64d,%I64d,%I64d",&i,&denom,&a,&b,&c,&d) == 6) {
				if (
					(i < 0) ||
					(i >= ANZCOEFF)
				) return -1; // not supported
				
				double da=a; da /= denom;
				double db=b; db /= denom;
				double dc=c; dc /= denom;
				double dd=d; dd /= denom;
				SETCOEFFCPLX2(coeff,i,da,db,dc,dd);
			} else
			if (sscanf(&tmp[1],"%i,%I64d,%I64d,%I64d",&i,&denom,&a,&b) == 4) {
				if (
					(i < 0) ||
					(i >= ANZCOEFF)
				) return -1; // not supported
				
				double da=a; da /= denom;
				double db=b; db /= denom;
				SETCOEFFCPLX2(coeff,i,da,da,db,db);
			} else {
				LOGMSG2("\nError. Coefficient file: Unknown line |%s|\n",
					tmp);
				exit(99);
			}
			
		} else {
			LOGMSG2("\nError. Coefficient file: Unknown line |%s|\n",
				tmp);
			exit(99);
		}
	} // while
	
	fclose(f);
	
	return 1;
}


// main

int32_t main(int32_t argc,char** argv) {
	// rounding mode necessary to be UPWARD
	// for use of custInt library
	CustRoundingUpwards *rd=new CustRoundingUpwards;
	
	flog=fopen("subdiv-core.log","at");
	fprintf(flog,"\n-----------------\n");
	printf("subdiv-core\n");
	
	char fn[1024];
	sprintf(fn,"_infkt.txt");
	
	// some standard parameters for the root finding method
	SubDivParam subject;
	CLEARCOEFF(subject.coeff)
	subject.startbox.rect.re.left =-2.0;
	subject.startbox.rect.re.right= 2.0;
	subject.startbox.rect.im.left =-2.0;
	subject.startbox.rect.im.right= 2.0;
	subject.maxrefinementdepth=16;
	STOPATROOT=ROOTSTOP_WHENFOUND;
	
	// parsing command-line parameters
	for(int32_t i=1;i<argc;i++) {
		upper(argv[i]);
		
		if (strstr(argv[i],"FILE=")==argv[i]) {
			sprintf(fn,&argv[i][5]);
		} else
		if (strstr(argv[i],"STOP=DEEPEST")==argv[i]) {
			STOPATROOT=ROOTSTOP_DEEPEST;
		} else
		if (strstr(argv[i],"STOP=WIDTH,")==argv[i]) {
			double a;
			if (sscanf(&argv[i][11],"%lf",&a) == 1) {
				STOPATROOT=ROOTSTOP_ATWIDTH;
				stopwidth=a;
			}
		} else
		if (strstr(argv[i],"RANGE=")==argv[i]) {
			double a,b,c,d;
			if (sscanf(&argv[i][6],"%lf,%lf,%lf,%lf",&a,&b,&c,&d) ==4) {
				subject.startbox.rect.re.left=a;
				subject.startbox.rect.re.right=b;
				subject.startbox.rect.im.left=c;
				subject.startbox.rect.im.right=d;
			}
		} else
		if (strstr(argv[i],"DEPTH=")==argv[i]) {
			int32_t a;
			if (sscanf(&argv[i][6],"%i",&a) == 1) {
				subject.maxrefinementdepth=a;
			}
		}
	} // i
	
	// print configuration
	LOGMSG("\nINFO");
	LOGMSG2("\n  max subdivision depth %i\n",subject.maxrefinementdepth);
	LOGMSG5("  range x%.20lg..%.20lg,y%.20lg..%.20lg\n",
		subject.startbox.rect.re.left,
		subject.startbox.rect.re.right,
		subject.startbox.rect.im.left,
		subject.startbox.rect.im.right
	);
	LOGMSG2("  input file %s\n",fn);
	
	if (STOPATROOT == ROOTSTOP_WHENFOUND) {
		LOGMSG("  stop box subdivision as soon as root is found\n");
	} else
	if (STOPATROOT == ROOTSTOP_DEEPEST) {
		LOGMSG("  subdivide till deepest level\n");
	} else
	if (STOPATROOT == ROOTSTOP_ATWIDTH) {
		LOGMSG2("  stop box subdivision at root and width smaller than %.5lg\n",stopwidth);
	}
	
	// load coefficients from file
	if (loadCoefficients(fn,subject.coeff) <= 0) {
		LOGMSG("\nError. Cooefficient file not found.\n");
		exit(99);
	}
	

	// /////////////////////////
	// main routine
	// /////////////////////////

	find_roots(subject);
	
	// /////////////////////////


	fclose(flog);
	delete rd;
	
    return 0;
}

