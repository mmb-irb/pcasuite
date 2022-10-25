/* matfit.f -- translated by f2c (version 20061008).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include <math.h>

/* Table of constant values */

static int c__3 = 3;
static int c__0 = 0;

/* Subroutine */ int matfit(int *n, float *xa, float *xb, float *r__, float *
	v, float *rmse, int *entry__)
{
    /* System generated locals */
    int i__1;

    /* Builtin functions */
    double sqrt(double);

    /* Local variables */
    static int i__, j, k;
    static double t, xn, cma[3], cmb[3], xni, umat[9]	/* was [3][3]
	    */, xasq, xbsq;
    extern /* Subroutine */ int qkfit_(double *, double *, float *,
	    int *);
    static double rtsum;


/*     SUBROUTINE TO FIT THE COORD SET XA(3,N) TO THE SET XB(3,N) */
/*     IN THE SENSE OF: */
/*            XA= R*XB +V */
/*     R IS A UNITARY 3.3 RIGHT HANDED ROTATION MATRIX */
/*     AND V IS THE OFFSET VECTOR. THIS IS AN EXACT SOLUTION */

/*     IF ENTRY IS LOGICALLY FALSE ONLY THE RMS COORDINATE ERROR */
/*     WILL BE RETURNED (BUT QUICKLY) */

/*    THIS SUBROUTINE IS A COMBINATION OF MCLACHLAN'S AND KABSCH'S */
/*    TECHNIQUES. SEE */
/*     KABSCH, W. ACTA CRYST A34, 827,1978 */
/*     MCLACHAN, A.D., J. MOL. BIO. NNN, NNNN 1978 */
/*     WRITTEN BY S.J. REMINGTON 11/78. */

/*     THIS SUBROUTINE USES THE IBM SSP EIGENVALUE ROUTINE 'EIGEN' */

    /* Parameter adjustments */
    xb -= 4;
    xa -= 4;
    r__ -= 4;
    --v;

    /* Function Body */
    xn = (double) (*n);
    xasq = 0.f;
    xbsq = 0.f;
    xni = 1.f / xn;

/*     ACCUMULATE UNCORRECTED (FOR C.M.) SUMS AND SQUARES */

    for (i__ = 1; i__ <= 3; ++i__) {
	cma[i__ - 1] = 0.f;
	cmb[i__ - 1] = 0.f;
	for (j = 1; j <= 3; ++j) {
/* L10: */
	    umat[i__ + j * 3 - 4] = 0.f;
	}

	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    for (k = 1; k <= 3; ++k) {
/* L20: */
		umat[i__ + k * 3 - 4] += xa[i__ + j * 3] * xb[k + j * 3];
	    }

	    t = xa[i__ + j * 3];
	    cma[i__ - 1] += t;
	    xasq += t * t;
	    t = xb[i__ + j * 3];
	    cmb[i__ - 1] += t;
	    xbsq += t * t;
/* L30: */
	}
/* L40: */
    }

/*     SUBTRACT CM OFFSETS */

    for (i__ = 1; i__ <= 3; ++i__) {
	xasq -= cma[i__ - 1] * cma[i__ - 1] * xni;
	xbsq -= cmb[i__ - 1] * cmb[i__ - 1] * xni;
	for (j = 1; j <= 3; ++j) {
	    umat[i__ + j * 3 - 4] = (umat[i__ + j * 3 - 4] - cma[i__ - 1] *
		    cmb[j - 1] * xni) * xni;
/* L50: */
	}
    }

/*     FIT IT */

    qkfit_(umat, &rtsum, &r__[4], entry__);
    *rmse = (xasq + xbsq) * xni - rtsum * 2.f;
    if (*rmse < 0.f) {
	*rmse = 0.f;
    }
    *rmse = sqrt(*rmse);

/*      CALCULATE OFFSET IF ENTRY=.TRUE. */

    if (! (*entry__)) {
	return 0;
    }

    for (i__ = 1; i__ <= 3; ++i__) {
	t = 0.f;
	for (j = 1; j <= 3; ++j) {
/* L60: */
	    t += r__[i__ + j * 3] * cmb[j - 1];
	}
	v[i__] = (cma[i__ - 1] - t) * xni;
/* L70: */
    }
    return 0;
} /* matfit_ */

double d_sign(double *a, double *b) {
  double ret;

  if (*a == 0.0) {
    ret = 0.0;
  } else {
    if (*b>=0.0) {
      ret = fabs (*a);
    } else {
      ret = -fabs (*a);
    }
  }

  return ret;
}

/* Subroutine */ int qkfit_(double *umat, double *rtsum, float *r__,
	int *entry__)
{
    /* Initialized data */

    static double eps = 1e-12;
    static double pi = 3.14159265358979;
    static double one = 1.;
    static double two = 2.;
    static double three = 3.;
    static double half = .5;
    static double third = .333333333;
    static double forthr = 1.333333333;
    static struct {
	double fill_1[1];
	double e_2[2];
	double fill_3[2];
	double e_4;
	double fill_5[3];
	} equiv_6 = { {0.}, {0., 0.}, {0., 0.}, 0., {0., 0., 0.} };


    /* System generated locals */
    int i__1;
    static double equiv_7[9];

    /* Builtin functions */
    double sqrt(double), d_sign(double *, double *), atan(
	    double), cos(double);

    /* Local variables */
#define a ((double *)&equiv_6)
#define b (equiv_7)
    static int i__, j, k;
    static double s, t, b1, b2, cc, b13, dd, b23;
    static int ia;
    static double b33, qq, rt;
#define aam ((double *)&equiv_6)
#define bam ((double *)&equiv_6 + 4)
#define cam ((double *)&equiv_6 + 8)
#define fam ((double *)&equiv_6 + 7)
#define gam ((double *)&equiv_6 + 6)
#define ham ((double *)&equiv_6 + 3)
    static double du11, du21, du31;
#define utr (equiv_7)
    static double diff;
    static int isig;
    static double detu, root[3];
    extern /* Subroutine */ int eigen_(double *, double *, int *,
	    int *);
    static double digav, theta, argsq;
    extern /* Subroutine */ int esort_(double *, double *, int *,
	    int *);
    static double cos3th;
#define usqmat ((double *)&equiv_6)


/*     THE 'EIGENVALUE ONLY' ENTRY WAS */
/*     ADAPTED FROM PROGRAM BY A.D. MCLACHAN 7/78 */



    /* Parameter adjustments */
    r__ -= 4;
    umat -= 4;

    /* Function Body */
    isig = 1;

/*      IF ENTRY IS .TRUE. GET OUT THE ROTATION MATRIX */

    if (*entry__) {
	goto L200;
    }

/*     CALC DET OF UMAT */

    du11 = umat[8] * umat[12] - umat[11] * umat[9];
    du21 = umat[11] * umat[6] - umat[5] * umat[12];
    du31 = umat[5] * umat[9] - umat[8] * umat[6];
    detu = umat[4] * du11 + umat[7] * du21 + umat[10] * du31;
/* %    WRITE(6,999) DETU */
/* %999 FORMAT(/(3F12.5)) */
    if (detu < 0.f) {
	isig = -1;
    }

/*     FORM USQMAT AS POSITIVE SEMI DEFINITE MATRIX */

    for (j = 1; j <= 3; ++j) {
	i__1 = j;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    usqmat[i__ + j * 3 - 4] = umat[i__ * 3 + 1] * umat[j * 3 + 1] +
		    umat[i__ * 3 + 2] * umat[j * 3 + 2] + umat[i__ * 3 + 3] *
		    umat[j * 3 + 3];
/* L105: */
	}
/* L110: */
    }
/* %    WRITE(6,999) USQMAT */

/*     REDUCE AVG OF DIAGONAL TERMS TO ZERO */

    digav = (*aam + *bam + *cam) * third;
/* %    WRITE(6,999) DIGAV */
    *aam -= digav;
    *bam -= digav;
    *cam -= digav;

/*     SETUP COEFFS OF SECULAR EQUATION OF MATRIX WITH TRACE ZERO */

    cc = *fam * *fam + *gam * *gam + *ham * *ham - *aam * *bam - *bam * *cam
	    - *cam * *aam;
    dd = *aam * *bam * *cam + two * (*fam * *gam * *ham) - *aam * *fam * *fam
	    - *bam * *gam * *gam - *cam * *ham * *ham;

/*     THE SECULAR EQN IS Y**3-CC*Y-DD=0  AND DD IS DET(USQMAT) */
/*     REDUCE THIS TO THE FORM COS**3-(3/4)COS- */
/*     (1/4)COS3THETA = 0 */
/*     WITH SOLUTIONS COSTHETA.  SO Y=QQ*COSTHETA */

    if (cc <= eps) {
	goto L115;
    }
    qq = sqrt(forthr * cc);
    cos3th = three * dd / (cc * qq);
    if (fabs(cos3th) > one) {
	cos3th = d_sign(&one, &cos3th);
    }

/*     FUNCTION ARCOS */

    if (cos3th != 0.f) {
	goto L1200;
    }
/* L1100: */
    theta = 1.570796327f;
    goto L1400;
L1200:
    argsq = cos3th * cos3th;
    theta = atan(sqrt(1.f - argsq) / cos3th);
    if (cos3th < 0.f) {
	theta = pi - theta;
    }
L1400:

/*     ROOTS IN ORDER OF SIZE GO 1,2,3 1 LARGEST */

    theta *= third;
    root[0] = qq * cos(theta);
    diff = half * sqrt(three * (qq * qq - root[0] * root[0]));
    root[1] = -root[0] * half + diff;
    root[2] = -root[0] * half - diff;
    goto L120;
L115:

/*     SPECIAL FOR TRIPLY DEGENERATE */

    root[0] = 0.f;
    root[1] = 0.f;
    root[2] = 0.f;
L120:
/*     ADD ON DIGAV AND TAKE SQRT */
    for (j = 1; j <= 3; ++j) {
	rt = root[j - 1] + digav;
	if (rt < eps) {
	    rt = 0.f;
	}
	root[j - 1] = sqrt(rt);
/* L125: */
    }
/* %    WRITE(6,999) ROOT */
/*     IF DETU IS <0 CHANGE SIGN OF ROOT(3) */
    if (isig == -1) {
	root[2] = -root[2];
    }
    *rtsum = root[0] + root[1] + root[2];
/* %    WRITE(6,999) RTSUM */
    return 0;

/*     THIS IS THE FANCY PART */

L200:

/*     FORM USQ = (UT).U    (IN UPPER TRIANGULAR SYMMETRIC STORAGE MODE) */

    for (i__ = 1; i__ <= 3; ++i__) {
	for (j = i__; j <= 3; ++j) {
	    t = 0.f;
	    for (k = 1; k <= 3; ++k) {
		t += umat[k + i__ * 3] * umat[k + j * 3];
/* L210: */
	    }
	    ia = i__ + (j * j - j) / 2;
	    utr[ia - 1] = t;
/* L220: */
	}
    }
/* %    WRITE(6,999) UTR */

/*     CALCULATE EIGENVALUES AND VECTORS */

    eigen_(utr, a, &c__3, &c__0);
    esort_(utr, a, &c__3, &c__0);
/* %    WRITE(6,999) UTR */

    root[0] = utr[0];
    root[1] = utr[2];
    root[2] = utr[5];
/* %    WRITE(6,999) ROOT */
/* %    WRITE(6,999) A */

/*     SET A3 = A1 CROSS A2 */
/*     ROOTS ARE IN ORDER R(1) >= R(2) >= R(3) >= 0 */

    a[6] = a[1] * a[5] - a[2] * a[4];
    a[7] = a[2] * a[3] - a[0] * a[5];
    a[8] = a[0] * a[4] - a[1] * a[3];
/* %    WRITE(6,999) A */

/*     VECTOR SET B=U.A */

    for (i__ = 1; i__ <= 3; ++i__) {
	for (j = 1; j <= 3; ++j) {
	    t = 0.f;
	    for (k = 1; k <= 3; ++k) {
/* L230: */
		t += umat[j + k * 3] * a[k + i__ * 3 - 4];
	    }
	    b[j + i__ * 3 - 4] = t;
/* L240: */
	}
    }

/*      NORMALIZE B1 AND B2 AND CALCULATE B3 = B1 CROSS B2 */

    b1 = sqrt(b[0] * b[0] + b[1] * b[1] + b[2] * b[2]);
    b2 = sqrt(b[3] * b[3] + b[4] * b[4] + b[5] * b[5]);
    for (i__ = 1; i__ <= 3; ++i__) {
	b[i__ - 1] /= b1;
/* L250: */
	b[i__ + 2] /= b2;
    }

/*      CHECK FOR LEFT HANDED ROTATION */

    b13 = b[1] * b[5] - b[2] * b[4];
    b23 = b[2] * b[3] - b[0] * b[5];
    b33 = b[0] * b[4] - b[1] * b[3];

    s = b13 * b[6] + b23 * b[7] + b33 * b[8];
    if (s < 0.f) {
	isig = -1;
    }
    b[6] = b13;
    b[7] = b23;
    b[8] = b33;
/* %    WRITE(6,999) B */

/*     CALCULATE ROTATION MATRIX R */

    for (i__ = 1; i__ <= 3; ++i__) {
	for (j = 1; j <= 3; ++j) {
	    t = 0.f;
	    for (k = 1; k <= 3; ++k) {
/* L260: */
		t += b[i__ + k * 3 - 4] * a[j + k * 3 - 4];
	    }
	    r__[i__ + j * 3] = t;
/* L270: */
	}
    }

/*     RMS ERROR */

    for (i__ = 1; i__ <= 3; ++i__) {
	if (root[i__ - 1] < 0.f) {
	    root[i__ - 1] = 0.f;
	}
	root[i__ - 1] = sqrt(root[i__ - 1]);
/* L280: */
    }

/*     CHANGE SIGN OF EVAL #3 IF LEFT HANDED */

    if (isig < 0) {
	root[2] = -root[2];
    }
    *rtsum = root[2] + root[1] + root[0];
    return 0;
} /* qkfit_ */

#undef usqmat
#undef utr
#undef ham
#undef gam
#undef fam
#undef cam
#undef bam
#undef aam
#undef b
#undef a


/* ---- SUBROUTINE TO COMPUTE EIGENVALUES & EIGENVECTORS OF A REAL SYMMETRIC */
/* ---- MATRIX, STOLEN FROM IBM SSP MANUAL (SEE P165) */
/* ---- DESCRIPTION OF PARAMETERS - */
/* ---- A - ORIGINAL MATRIX STORED COLUMNWISE AS UPPER TRIANGLE ONLY, */
/* ---- I.E. "STORAGE MODE" = 1.  EIGENVALUES ARE WRITTEN INTO DIAGONAL */
/* ---- ELEMENTS OF A  I.E.  A(1)  A(3)  A(6)  FOR A 3*3 MATRIX. */
/* ---- R - RESULTANT MATRIX OF EIGENVECTORS STORED COLUMNWISE IN SAME */
/* ---- ORDER AS EIGENVALUES. */
/* ---- N - ORDER OF MATRICES A & R. */
/* ---- MV = 0 TO COMPUTE EIGENVALUES & EIGENVECTORS. */
/* Subroutine */ int eigen_(double *a, double *r__, int *n,
	int *mv)
{
    /* System generated locals */
    int i__1, i__2;
    double d__1;

    /* Builtin functions */
    double sqrt(double);

    /* Local variables */
    static int i__, j, l, m;
    static double x, y;
    static int ia, ij, il, im, ll, lm, iq, mm, lq, mq, ind, ilq, imq, ilr,
	     imr;
    static double thr, cosx, sinx, cosx2, sinx2, range, anorm, sincs,
	    anrmx;

/* ---- FOR DOUBLE PRECISION, SQRT IN STATEMENTS 40,68,75&78 MUST BE DSQRT, */
/* ---- ABS IN 62 MUST BE DABS AND 1.E-6 IN 5 MUST BE 1.D-12 . */
/* 5	RANGE=1.E-6 */
    /* Parameter adjustments */
    --r__;
    --a;

    /* Function Body */
/* L5: */
    range = 1e-12;
    if (*mv - 1 != 0) {
	goto L10;
    } else {
	goto L25;
    }
L10:
    iq = -(*n);
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	iq += *n;
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    ij = iq + i__;
	    r__[ij] = 0.f;
	    if (i__ - j != 0) {
		goto L20;
	    } else {
		goto L15;
	    }
L15:
	    r__[ij] = 1.f;
L20:
	    ;
	}
    }
/* ---- INITIAL AND FINAL NORMS (ANORM & ANRMX) */
L25:
    anorm = 0.f;
    i__2 = *n;
    for (i__ = 1; i__ <= i__2; ++i__) {
	i__1 = *n;
	for (j = i__; j <= i__1; ++j) {
	    if (i__ - j != 0) {
		goto L30;
	    } else {
		goto L35;
	    }
L30:
	    ia = i__ + (j * j - j) / 2;
/* Computing 2nd power */
	    d__1 = a[ia];
	    anorm += d__1 * d__1;
L35:
	    ;
	}
    }
    if (anorm <= 0.) {
	goto L165;
    } else {
	goto L40;
    }
/* 40	ANORM=SQRT(2.*ANORM) */
L40:
    anorm = sqrt(anorm * 2.f);
    anrmx = anorm * range / *n;
/* ---- INITIALIZE INDICATORS AND COMPUTE THRESHOLD */
    ind = 0;
    thr = anorm;
L45:
    thr /= *n;
L50:
    l = 1;
L55:
    m = l + 1;
/* ---- COMPUTE SIN & COS */
L60:
    mq = (m * m - m) / 2;
    lq = (l * l - l) / 2;
    lm = l + mq;
/* L62: */
    if ((d__1 = a[lm], fabs(d__1)) - thr >= 0.) {
	goto L65;
    } else {
	goto L130;
    }
L65:
    ind = 1;
    ll = l + lq;
    mm = m + mq;
    x = (a[ll] - a[mm]) * .5f;
/* 68	Y=-A(LM)/SQRT(A(LM)**2+X*X) */
/* L68: */
/* Computing 2nd power */
    d__1 = a[lm];
    y = -a[lm] / sqrt(d__1 * d__1 + x * x);
    if (x >= 0.) {
	goto L75;
    } else {
	goto L70;
    }
L70:
    y = -y;
/* 75	SINX=Y/SQRT(2.*(1.+(SQRT(1.-Y*Y)))) */
L75:
    sinx = y / sqrt((sqrt(1.f - y * y) + 1.f) * 2.f);
/* Computing 2nd power */
    d__1 = sinx;
    sinx2 = d__1 * d__1;
/* 78	COSX=SQRT(1.-SINX2) */
/* L78: */
    cosx = sqrt(1.f - sinx2);
/* Computing 2nd power */
    d__1 = cosx;
    cosx2 = d__1 * d__1;
    sincs = sinx * cosx;
/* ---- ROTATE L & M COLUMNS */
    ilq = *n * (l - 1);
    imq = *n * (m - 1);
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	iq = (i__ * i__ - i__) / 2;
	if (i__ - l != 0) {
	    goto L80;
	} else {
	    goto L115;
	}
L80:
	if ((i__2 = i__ - m) < 0) {
	    goto L85;
	} else if (i__2 == 0) {
	    goto L115;
	} else {
	    goto L90;
	}
L85:
	im = i__ + mq;
	goto L95;
L90:
	im = m + iq;
L95:
	if (i__ - l >= 0) {
	    goto L105;
	} else {
	    goto L100;
	}
L100:
	il = i__ + lq;
	goto L110;
L105:
	il = l + iq;
L110:
	x = a[il] * cosx - a[im] * sinx;
	a[im] = a[il] * sinx + a[im] * cosx;
	a[il] = x;
L115:
	if (*mv - 1 != 0) {
	    goto L120;
	} else {
	    goto L125;
	}
L120:
	ilr = ilq + i__;
	imr = imq + i__;
	x = r__[ilr] * cosx - r__[imr] * sinx;
	r__[imr] = r__[ilr] * sinx + r__[imr] * cosx;
	r__[ilr] = x;
L125:
	;
    }
    x = a[lm] * 2.f * sincs;
    y = a[ll] * cosx2 + a[mm] * sinx2 - x;
    x = a[ll] * sinx2 + a[mm] * cosx2 + x;
    a[lm] = (a[ll] - a[mm]) * sincs + a[lm] * (cosx2 - sinx2);
    a[ll] = y;
    a[mm] = x;
/* ---- TESTS FOR COMPLETION */
/* ---- TEST FOR M = LAST COLUMN */
L130:
    if (m - *n != 0) {
	goto L135;
    } else {
	goto L140;
    }
L135:
    ++m;
    goto L60;
/* ---- TEST FOR L =PENULTIMATE COLUMN */
L140:
    if (l - (*n - 1) != 0) {
	goto L145;
    } else {
	goto L150;
    }
L145:
    ++l;
    goto L55;
L150:
    if (ind - 1 != 0) {
	goto L160;
    } else {
	goto L155;
    }
L155:
    ind = 0;
    goto L50;
/* ---- COMPARE THRESHOLD WITH FINAL NORM */
L160:
    if (thr - anrmx <= 0.) {
	goto L165;
    } else {
	goto L45;
    }
L165:
    return 0;
/* ---- SORT EIGENVALUES AND EIGENVECTORS IN DESCENDING ORDER OF EIGENVALUES */
} /* eigen_ */

/* Subroutine */ int esort_(double *a, double *r__, int *n,
	int *mv)
{
    /* System generated locals */
    int i__1, i__2, i__3;

    /* Local variables */
    static int i__, j, k;
    static double x;
    static int ll, iq, jq, mm, ilr, imr;

    /* Parameter adjustments */
    --r__;
    --a;

    /* Function Body */
    iq = -(*n);
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	iq += *n;
	ll = i__ + (i__ * i__ - i__) / 2;
	jq = *n * (i__ - 2);
	i__2 = *n;
	for (j = i__; j <= i__2; ++j) {
	    jq += *n;
	    mm = j + (j * j - j) / 2;
	    if (a[ll] - a[mm] >= 0.) {
		goto L185;
	    } else {
		goto L170;
	    }
L170:
	    x = a[ll];
	    a[ll] = a[mm];
	    a[mm] = x;
	    if (*mv - 1 != 0) {
		goto L175;
	    } else {
		goto L185;
	    }
L175:
	    i__3 = *n;
	    for (k = 1; k <= i__3; ++k) {
		ilr = iq + k;
		imr = jq + k;
		x = r__[ilr];
		r__[ilr] = r__[imr];
/* L180: */
		r__[imr] = x;
	    }
L185:
	    ;
	}
    }
    return 0;
} /* esort_ */

