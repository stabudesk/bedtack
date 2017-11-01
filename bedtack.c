/* bgmergbl: bedgraph merge bed lines: the first bed refers to an over all bed file and the second bed will be treated on a line by line basis
 * essentially this programs summarises the first bed file according to the lines in the second */
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<unistd.h> // required for optopt, opterr and optarg.

#ifdef DBG
#define GBUF 4
#define WBUF 4
#else
#define GBUF 32
#define WBUF 32
#endif

#define NUMBUCKETS 20

// the following is the way we cut out columns hat have nothing in them.
#define MXCOL2VIEW 4

#define CONDREALLOC(x, b, c, a, t); \
	if((x)>=((b)-1)) { \
		(b) += (c); \
		(a)=realloc((a), (b)*sizeof(t)); \
		memset(((a)+(b)-(c)), 0, (c)*sizeof(t)); \
	}

typedef unsigned char boole;

typedef struct /* ia_t integer array type, includes iab the buffer */
{
	int *a;
	unsigned b /* int array buf */, z /* int array size*/;
} ia_t;

typedef struct  /* opt_t, a struct for the options */
{
	boole dflg; /* details / information only */
	boole sflg; /* split outout in two files */
	char *istr; /* first bedgraph file, the target of the filtering by the second */
	char *fstr; /* the name of the second bedgraph file */
	char *ustr; /* the name of a file with the list of elements to be unified */
	char *pstr; /* depth file name */
} opt_t;

typedef struct /* i4_t */
{
	int sc; /* number of same chromosome (occurences?) */
	float mc; /* min signal value */
	int b1i; /* index of the 1st bgr_t, which satisfies the conditions */
	int lgbi; /* last good bgr_t index */
} i4_t; /* 4 vals of some sort? */

typedef struct /* bgr_t */
{
	char *n;
	size_t nsz; /* size of the name r ID field */
	long c[2]; /* coords: 1) start 2) end */
	float co; /* signal value */
} bgr_t; /* bedgraph row type */

typedef struct /* bgr_t2 */
{
	char *n;
	size_t nsz; /* size of the name r ID field */
	long c[2]; /* coords: 1) start 2) end */
	char *f; /* f for feature .. 4th col */
	size_t fsz; /* size of the feature field*/
} bgr_t2; /* bedgraph row type 2i. column is the feature */

typedef struct /* words_t: file with only single words per line */
{
	char *n;
	size_t nsz; /* size of the name r ID field */
} words_t; /* bedgraph row type */

typedef struct /* dpf_t : depth file type ... just chr name, pos and read quant */
{
	char *n;
	size_t nsz; /* size of the name r ID field */
	long p; /* position */
	int d; /* depth reading */
} dpf_t;

typedef struct /* wseq_t */
{
	size_t *wln;
	size_t wsbuf;
	size_t quan;
	size_t lbuf; /* a buffer for the number of lines */
	size_t numl; /* number of lines, i.e. rows */
	size_t *wpla; /* words per line array: the number of words on each line */
} wseq_t;

wseq_t *create_wseq_t(size_t initsz)
{
	wseq_t *words=malloc(sizeof(wseq_t));
	words->wsbuf = initsz;
	words->quan = initsz;
	words->wln=calloc(words->wsbuf, sizeof(size_t));
	words->lbuf=WBUF;
	words->numl=0;
	words->wpla=calloc(words->lbuf, sizeof(size_t));
	return words;
}

int catchopts(opt_t *opts, int oargc, char **oargv)
{
	int c;
	opterr = 0;

	while ((c = getopt (oargc, oargv, "dsi:f:u:p:")) != -1)
		switch (c) {
			case 'd':
				opts->dflg = 1;
				break;
			case 's':
				opts->sflg = 1;
				break;
			case 'i':
				opts->istr = optarg;
				break;
			case 'f':
				opts->fstr = optarg;
				break;
			case 'u': /* unify certain bed2 elements into one file */
				opts->ustr = optarg;
				break;
			case 'p': /* depth file */
				opts->pstr = optarg;
				break;
			case '?':
				fprintf (stderr, "Unknown option character `\\x%x'.\n", optopt);
				return 1;
			default:
				fprintf (stderr, "Wrong arguments. Please launch without arguments to see help file.\n");
				exit(EXIT_FAILURE);
		}
	return 0;
}

void free_wseq(wseq_t *wa)
{
	free(wa->wln);
	free(wa->wpla);
	free(wa);
}

words_t *processwordf(char *fname, int *m, int *n)
{
	/* In order to make no assumptions, the file is treated as lines containing the same amount of words each,
	 * except for lines starting with #, which are ignored (i.e. comments). These words are checked to make sure they contain only floating number-type
	 * characters [0123456789+-.] only, one string variable is continually written over and copied into a growing floating point array each time */

	/* declarations */
	FILE *fp=fopen(fname,"r");
	int i;
	size_t couc /*count chars per line */, couw=0 /* count words */, oldcouw = 0;
	int c;
	boole inword=0;
	wseq_t *wa=create_wseq_t(GBUF);
	size_t bwbuf=WBUF;
	char *bufword=calloc(bwbuf, sizeof(char)); /* this is the string we'll keep overwriting. */

	words_t *bedword=malloc(GBUF*sizeof(words_t));

	while( (c=fgetc(fp)) != EOF) { /* grab a char */
		if( (c== '\n') | (c == ' ') | (c == '\t') | (c=='#')) { /* word closing events */
			if( inword==1) { /* first word closing event */
				wa->wln[couw]=couc;
				bufword[couc++]='\0';
				bufword = realloc(bufword, couc*sizeof(char)); /* normalize */
				/* for the struct, we want to know if it's the first word in a line, like so: */
				if(couw==oldcouw) {
					bedword[wa->numl].n=malloc(couc*sizeof(char));
					bedword[wa->numl].nsz=couc;
					strcpy(bedword[wa->numl].n, bufword);
				}
				couc=0;
				couw++;
			}
			if(c=='#') { /* comment case */
				while( (c=fgetc(fp)) != '\n') ;
				continue;
			} else if(c=='\n') { /* end of a line */
				if(wa->numl == wa->lbuf-1) { /* enought space in our current array? */
					wa->lbuf += WBUF;
					wa->wpla=realloc(wa->wpla, wa->lbuf*sizeof(size_t));
					bedword=realloc(bedword, wa->lbuf*sizeof(words_t));
					memset(wa->wpla+(wa->lbuf-WBUF), 0, WBUF*sizeof(size_t));
				}
				wa->wpla[wa->numl] = couw-oldcouw; /* number of words in current line */
				if(couw-oldcouw >4) {
					printf("Error, each row cannot exceed 4 words: revise your input file\n"); 
					/* need to release all memory too */
					free_wseq(wa);
					exit(EXIT_FAILURE);
				}
				oldcouw=couw; /* restart words per line count */
				wa->numl++; /* brand new line coming up */
			}
			inword=0;
		} else if(inword==0) { /* deal with first character of new word, + and - also allowed */
			if(couw == wa->wsbuf-1) {
				wa->wsbuf += GBUF;
				wa->wln=realloc(wa->wln, wa->wsbuf*sizeof(size_t));
				bedword=realloc(bedword, wa->wsbuf*sizeof(words_t));
				for(i=wa->wsbuf-GBUF;i<wa->wsbuf;++i)
					wa->wln[i]=0;
			}
			couc=0;
			bwbuf=WBUF;
			bufword=realloc(bufword, bwbuf*sizeof(char)); /* don't bother with memset, it's not necessary */
			bufword[couc++]=c; /* no need to check here, it's the first character */
			inword=1;
		} else {
			if(couc == bwbuf-1) { /* the -1 so that we can always add and extra (say 0) when we want */
				bwbuf += WBUF;
				bufword = realloc(bufword, bwbuf*sizeof(char));
			}
			bufword[couc++]=c;
		}

	} /* end of big for statement */
	fclose(fp);
	free(bufword);

	/* normalization stage */
	wa->quan=couw;
	wa->wln = realloc(wa->wln, wa->quan*sizeof(size_t)); /* normalize */
	bedword = realloc(bedword, wa->quan*sizeof(words_t)); /* normalize */
	wa->wpla= realloc(wa->wpla, wa->numl*sizeof(size_t));

	*m= wa->numl;
	int k=wa->wpla[0];
	for(i=1;i<wa->numl;++i)
		if(k != wa->wpla[i])
			printf("Warning: Numcols is not uniform at %i words per line on all lines. This file has one with %zu.\n", k, wa->wpla[i]); 
	*n= k; 
	free_wseq(wa);

	return bedword;
}

bgr_t *processinpf(char *fname, int *m, int *n)
{
	/* In order to make no assumptions, the file is treated as lines containing the same amount of words each,
	 * except for lines starting with #, which are ignored (i.e. comments). These words are checked to make sure they contain only floating number-type
	 * characters [0123456789+-.] only, one string variable is continually written over and copied into a growing floating point array each time */

	/* declarations */
	FILE *fp=fopen(fname,"r");
	int i;
	size_t couc /*count chars per line */, couw=0 /* count words */, oldcouw = 0;
	int c;
	boole inword=0;
	wseq_t *wa=create_wseq_t(GBUF);
	size_t bwbuf=WBUF;
	char *bufword=calloc(bwbuf, sizeof(char)); /* this is the string we'll keep overwriting. */

	bgr_t *bgrow=malloc(GBUF*sizeof(bgr_t));

	while( (c=fgetc(fp)) != EOF) { /* grab a char */
		if( (c== '\n') | (c == ' ') | (c == '\t') | (c=='#')) { /* word closing events */
			if( inword==1) { /* first word closing event */
				wa->wln[couw]=couc;
				bufword[couc++]='\0';
				bufword = realloc(bufword, couc*sizeof(char)); /* normalize */
				/* for the struct, we want to know if it's the first word in a line, like so: */
				if(couw==oldcouw) {
					bgrow[wa->numl].n=malloc(couc*sizeof(char));
					bgrow[wa->numl].nsz=couc;
					strcpy(bgrow[wa->numl].n, bufword);
				} else if((couw-oldcouw)<3) /* it's not the first word, and it's 1st and second col */
					bgrow[wa->numl].c[couw-oldcouw-1]=atol(bufword);
				else if( (couw-oldcouw)==3) { // assume float
					bgrow[wa->numl].co=atof(bufword);
				}
				couc=0;
				couw++;
			}
			if(c=='#') { /* comment case */
				while( (c=fgetc(fp)) != '\n') ;
				continue;
			} else if(c=='\n') { /* end of a line */
				if(wa->numl == wa->lbuf-1) { /* enought space in our current array? */
					wa->lbuf += WBUF;
					wa->wpla=realloc(wa->wpla, wa->lbuf*sizeof(size_t));
					bgrow=realloc(bgrow, wa->lbuf*sizeof(bgr_t));
					memset(wa->wpla+(wa->lbuf-WBUF), 0, WBUF*sizeof(size_t));
				}
				wa->wpla[wa->numl] = couw-oldcouw; /* number of words in current line */
				if(couw-oldcouw >4) {
					printf("Error, each row cannot exceed 4 words: revise your input file\n"); 
					/* need to release all memory too */
					free_wseq(wa);
					exit(EXIT_FAILURE);
				}
				oldcouw=couw; /* restart words per line count */
				wa->numl++; /* brand new line coming up */
			}
			inword=0;
		} else if(inword==0) { /* deal with first character of new word, + and - also allowed */
			if(couw == wa->wsbuf-1) {
				wa->wsbuf += GBUF;
				wa->wln=realloc(wa->wln, wa->wsbuf*sizeof(size_t));
				bgrow=realloc(bgrow, wa->wsbuf*sizeof(bgr_t));
				for(i=wa->wsbuf-GBUF;i<wa->wsbuf;++i)
					wa->wln[i]=0;
			}
			couc=0;
			bwbuf=WBUF;
			bufword=realloc(bufword, bwbuf*sizeof(char)); /* don't bother with memset, it's not necessary */
			bufword[couc++]=c; /* no need to check here, it's the first character */
			inword=1;
		} else {
			if(couc == bwbuf-1) { /* the -1 so that we can always add and extra (say 0) when we want */
				bwbuf += WBUF;
				bufword = realloc(bufword, bwbuf*sizeof(char));
			}
			bufword[couc++]=c;
		}

	} /* end of big for statement */
	fclose(fp);
	free(bufword);

	/* normalization stage */
	wa->quan=couw;
	wa->wln = realloc(wa->wln, wa->quan*sizeof(size_t)); /* normalize */
	bgrow = realloc(bgrow, wa->quan*sizeof(bgr_t)); /* normalize */
	wa->wpla= realloc(wa->wpla, wa->numl*sizeof(size_t));

	*m= wa->numl;
	int k=wa->wpla[0];
	for(i=1;i<wa->numl;++i)
		if(k != wa->wpla[i])
			printf("Warning: Numcols is not uniform at %i words per line on all lines. This file has one with %zu.\n", k, wa->wpla[i]); 
	*n= k; 
	free_wseq(wa);

	return bgrow;
}

bgr_t2 *processinpf2(char *fname, int *m, int *n) /*fourth column is string, other columns to be ignored */
{
	/* In order to make no assumptions, the file is treated as lines containing the same amount of words each,
	 * except for lines starting with #, which are ignored (i.e. comments). These words are checked to make sure they contain only floating number-type
	 * characters [0123456789+-.] only, one string variable is continually written over and copied into a growing floating point array each time */

	/* declarations */
	FILE *fp=fopen(fname,"r");
	int i;
	size_t couc /*count chars per line */, couw=0 /* count words */, oldcouw = 0;
	int c;
	boole inword=0;
	wseq_t *wa=create_wseq_t(GBUF);
	size_t bwbuf=WBUF;
	char *bufword=calloc(bwbuf, sizeof(char)); /* this is the string we'll keep overwriting. */

	bgr_t2 *bgrow=malloc(GBUF*sizeof(bgr_t2));

	while( (c=fgetc(fp)) != EOF) { /* grab a char */
		if( (c== '\n') | (c == ' ') | (c == '\t') | (c=='#')) { /* word closing events */
			if( inword==1) { /* first word closing event */
				wa->wln[couw]=couc;
				bufword[couc++]='\0';
				bufword = realloc(bufword, couc*sizeof(char)); /* normalize */
				/* for the struct, we want to know if it's the first word in a line, like so: */
				if(couw==oldcouw) {
					bgrow[wa->numl].n=malloc(couc*sizeof(char));
					bgrow[wa->numl].nsz=couc;
					strcpy(bgrow[wa->numl].n, bufword);
				} else if((couw-oldcouw)<3) /* it's not the first word, and it's 1st and second col */
					bgrow[wa->numl].c[couw-oldcouw-1]=atol(bufword);
				else if( (couw-oldcouw)==3) { // assume float
					bgrow[wa->numl].f=malloc(couc*sizeof(char));
					bgrow[wa->numl].fsz=couc;
					strcpy(bgrow[wa->numl].f, bufword);
				}
				couc=0;
				couw++;
			}
			if(c=='#') { /* comment case */
				while( (c=fgetc(fp)) != '\n') ;
				continue;
			} else if(c=='\n') { /* end of a line */
				if(wa->numl == wa->lbuf-1) { /* enought space in our current array? */
					wa->lbuf += WBUF;
					wa->wpla=realloc(wa->wpla, wa->lbuf*sizeof(size_t));
					bgrow=realloc(bgrow, wa->lbuf*sizeof(bgr_t));
					memset(wa->wpla+(wa->lbuf-WBUF), 0, WBUF*sizeof(size_t));
				}
				wa->wpla[wa->numl] = couw-oldcouw; /* number of words in current line */
				oldcouw=couw; /* restart words per line count */
				wa->numl++; /* brand new line coming up */
			}
			inword=0;
		} else if(inword==0) { /* deal with first character of new word, + and - also allowed */
			if(couw == wa->wsbuf-1) {
				wa->wsbuf += GBUF;
				wa->wln=realloc(wa->wln, wa->wsbuf*sizeof(size_t));
				bgrow=realloc(bgrow, wa->wsbuf*sizeof(bgr_t2));
				for(i=wa->wsbuf-GBUF;i<wa->wsbuf;++i)
					wa->wln[i]=0;
			}
			couc=0;
			bwbuf=WBUF;
			bufword=realloc(bufword, bwbuf*sizeof(char)); /* don't bother with memset, it's not necessary */
			bufword[couc++]=c; /* no need to check here, it's the first character */
			inword=1;
		} else {
			if(couc == bwbuf-1) { /* the -1 so that we can always add and extra (say 0) when we want */
				bwbuf += WBUF;
				bufword = realloc(bufword, bwbuf*sizeof(char));
			}
			bufword[couc++]=c;
		}

	} /* end of big for statement */
	fclose(fp);
	free(bufword);

	/* normalization stage */
	wa->quan=couw;
	wa->wln = realloc(wa->wln, wa->quan*sizeof(size_t)); /* normalize */
	bgrow = realloc(bgrow, wa->quan*sizeof(bgr_t2)); /* normalize */
	wa->wpla= realloc(wa->wpla, wa->numl*sizeof(size_t));

	*m= wa->numl;
	int k=wa->wpla[0];
	for(i=1;i<wa->numl;++i)
		if(k != wa->wpla[i])
			printf("Warning: Numcols is not uniform at %i words per line on all lines. This file has one with %zu.\n", k, wa->wpla[i]); 
	*n= k; 
	free_wseq(wa);

	return bgrow;
}

dpf_t *processdpf(char *fname, int *m, int *n) /*fourth column is string, other columns to be ignored */
{
	/* In order to make no assumptions, the file is treated as lines containing the same amount of words each,
	 * except for lines starting with #, which are ignored (i.e. comments). These words are checked to make sure they contain only floating number-type
	 * characters [0123456789+-.] only, one string variable is continually written over and copied into a growing floating point array each time */

	/* declarations */
	FILE *fp=fopen(fname,"r");
	int i;
	size_t couc /*count chars per line */, couw=0 /* count words */, oldcouw = 0;
	int c;
	boole inword=0;
	wseq_t *wa=create_wseq_t(GBUF);
	size_t bwbuf=WBUF;
	char *bufword=calloc(bwbuf, sizeof(char)); /* this is the string we'll keep overwriting. */

	dpf_t *dpf=malloc(GBUF*sizeof(dpf_t));

	while( (c=fgetc(fp)) != EOF) { /* grab a char */
		if( (c== '\n') | (c == ' ') | (c == '\t') | (c=='#')) { /* word closing events */
			if( inword==1) { /* first word closing event */
				wa->wln[couw]=couc;
				bufword[couc++]='\0';
				bufword = realloc(bufword, couc*sizeof(char)); /* normalize */
				/* for the struct, we want to know if it's the first word in a line, like so: */
				if(couw==oldcouw) {
					dpf[wa->numl].n=malloc(couc*sizeof(char));
					dpf[wa->numl].nsz=couc;
					strcpy(dpf[wa->numl].n, bufword);
				} else if((couw-oldcouw)==1) /* it's not the first word, and it's 1st and second col */
					dpf[wa->numl].p=atol(bufword);
				else if((couw-oldcouw)==2) /* it's not the first word, and it's 1st and second col */
					dpf[wa->numl].d=atoi(bufword);
				couc=0;
				couw++;
			}
			if(c=='#') { /* comment case */
				while( (c=fgetc(fp)) != '\n') ;
				continue;
			} else if(c=='\n') { /* end of a line */
				if(wa->numl == wa->lbuf-1) { /* enought space in our current array? */
					wa->lbuf += WBUF;
					wa->wpla=realloc(wa->wpla, wa->lbuf*sizeof(size_t));
					dpf=realloc(dpf, wa->lbuf*sizeof(bgr_t));
					memset(wa->wpla+(wa->lbuf-WBUF), 0, WBUF*sizeof(size_t));
				}
				wa->wpla[wa->numl] = couw-oldcouw; /* number of words in current line */
				oldcouw=couw; /* restart words per line count */
				wa->numl++; /* brand new line coming up */
			}
			inword=0;
		} else if(inword==0) { /* deal with first character of new word, + and - also allowed */
			if(couw == wa->wsbuf-1) {
				wa->wsbuf += GBUF;
				wa->wln=realloc(wa->wln, wa->wsbuf*sizeof(size_t));
				dpf=realloc(dpf, wa->wsbuf*sizeof(dpf_t));
				for(i=wa->wsbuf-GBUF;i<wa->wsbuf;++i)
					wa->wln[i]=0;
			}
			couc=0;
			bwbuf=WBUF;
			bufword=realloc(bufword, bwbuf*sizeof(char)); /* don't bother with memset, it's not necessary */
			bufword[couc++]=c; /* no need to check here, it's the first character */
			inword=1;
		} else {
			if(couc == bwbuf-1) { /* the -1 so that we can always add and extra (say 0) when we want */
				bwbuf += WBUF;
				bufword = realloc(bufword, bwbuf*sizeof(char));
			}
			bufword[couc++]=c;
		}

	} /* end of big for statement */
	fclose(fp);
	free(bufword);

	/* normalization stage */
	wa->quan=couw;
	wa->wln = realloc(wa->wln, wa->quan*sizeof(size_t)); /* normalize */
	dpf = realloc(dpf, wa->quan*sizeof(dpf_t)); /* normalize */
	wa->wpla= realloc(wa->wpla, wa->numl*sizeof(size_t));

	*m= wa->numl;
	int k=wa->wpla[0];
	for(i=1;i<wa->numl;++i)
		if(k != wa->wpla[i])
			printf("Warning: Numcols is not uniform at %i words per line on all lines. This file has one with %zu.\n", k, wa->wpla[i]); 
	*n= k; 
	free_wseq(wa);

	return dpf;
}

void prtbd2ia(bgr_t2 *bed2, int n, ia_t *ia)
{
	int i, j;
	for(i=0;i<ia->z;++i) {
		for(j=0;j<n;++j) {
			if(j==0)
				printf("%s ", bed2[ia->a[i]].n);
			else if(j==3)
				printf("%s ", bed2[ia->a[i]].f);
			else
				printf("%li ", bed2[ia->a[i]].c[j-1]);
			}
			printf("\n"); 
	}
	return;
}

void bed2in2(char *bed2fn, bgr_t2 *bed2, int m, int n, ia_t *ia) // split into 2 files
{
	int i, j, k=0;
	size_t lfn=strlen(bed2fn);
	char *outfn1=calloc(4+lfn ,sizeof(char));
	char *outfn2=calloc(4+lfn, sizeof(char));
	int rootsz=(int)(strchr(bed2fn, '.')-bed2fn);
	sprintf(outfn1, "%.*s_p1.bed", rootsz, bed2fn);
	sprintf(outfn2, "%.*s_p2.bed", rootsz, bed2fn);
	FILE *of1=fopen(outfn1, "w");
	FILE *of2=fopen(outfn2, "w");
	printf("bgr_t is %i rows by %i columns and is as follows:\n", m, n); 
	for(i=0;i<m;++i) {
		if(i==ia->a[k]){
			for(j=0;j<n;++j) {
				if(j==0)
					fprintf(of2, "%s\t", bed2[i].n);
				else if(j==3)
					fprintf(of2, "%s\n", bed2[i].f);
				else
					fprintf(of2, "%li\t", bed2[i].c[j-1]);
			}
			k++;
		} else {
			for(j=0;j<n;++j) {
				if(j==0)
					fprintf(of1, "%s\t", bed2[i].n);
				else if(j==3)
					fprintf(of1, "%s\n", bed2[i].f);
				else
					fprintf(of1, "%li\t", bed2[i].c[j-1]);
			}
		}
	}
	fclose(of1);
	fclose(of2);
	free(outfn1);
	free(outfn2);
	return;
}

void prtobed(bgr_t *bgrow, int m, int n, float minsig) // print over bed ... a value that is over a certain signal
{
	int i, j;
	printf("bgr_t is %i rows by %i columns and is as follows:\n", m, n); 
	for(i=0;i<m;++i) {
		if(bgrow[i].co >= minsig) {
			for(j=0;j<n;++j) {
				if(j==0)
					printf("%s ", bgrow[i].n);
				else if(j==3)
					printf("%2.6f ", bgrow[i].co);
				else
					printf("%li ", bgrow[i].c[j-1]);
			}
			printf("\n"); 
		}
	}
	return;
}

int *hist_co(bgr_t *bgrow, int m, float mxco, float mnco, int numbuckets)
{
	int i, j;
	float step=(mxco-mnco)/(float)numbuckets;
	float *bucketlimarr=malloc((numbuckets-1)*sizeof(float));
	int *bucketarr=calloc(numbuckets, sizeof(int));
	bucketlimarr[0]=step+mnco;
	for(i=1;i<numbuckets-1;++i) 
		bucketlimarr[i]=bucketlimarr[i-1]+step;

	for(i=0;i<m;++i)
		if(bgrow[i].co>=bucketlimarr[numbuckets-2]) {
			bucketarr[numbuckets-1]++;
			continue;
		} else {
			for(j=0;j<numbuckets-1;++j)
				if(bgrow[i].co < bucketlimarr[j]) {
					bucketarr[j]++;
					break;
				}
		}
	free(bucketlimarr);
	return bucketarr;
}

void prthist(char *histname, int *bucketarr, int numbuckets, int m, float mxco, float mnco)
{
	int i;
	printf("%s value %d-bin hstgrm for: %-24.24s (totels=%04i):\n", histname, numbuckets, histname, m); 
	printf("minval=%4.6f<-", mnco); 
	for(i=0;i<numbuckets;++i) 
		printf("| %i ", bucketarr[i]);
	printf("|->maxval=%4.6f\n", mxco); 
	return;
}

void prtdets(bgr_t *bgrow, int m, int n, char *label)
{
	int i;
	float mxco=.0, mnco=10e20;
	printf("bgr_t is %i rows by %i columns and is as follows:\n", m, n); 
	for(i=0;i<m;++i) {
		if(bgrow[i].co > mxco)
			mxco=bgrow[i].co;
		if(bgrow[i].co < mnco)
			mnco = bgrow[i].co;
	}
	int *hco=hist_co(bgrow, m, mxco, mnco, NUMBUCKETS);
	prthist(label, hco, NUMBUCKETS, m, mxco, mnco);
	free(hco);
	return;
}

void prtbed2s(bgr_t2 *bed2, int m, int n, words_t *bedword, int m3, int n3, char *label)
{
	/* TODO what you want is a copy of the data structure */
	int i, j, k;
	boole foundifeat;
	printf("Separated feature file %s is %i rows by %i columns and is as follows:\n", label, m, n); 
	for(i=0;i<m;++i) {
		foundifeat=0;
		for(k=0;k<m3;++k) {
			if(!strcmp(bedword[k].n, bed2[i].f) ) {
				foundifeat=1;
				for(j=0;j<n;++j) {
					if(j==0)
						printf("%s ", bed2[i].n);
					else if(j==3)
						printf("%s ", bed2[i].f);
					else
						printf("%li ", bed2[i].c[j-1]);
				}
				printf("\n"); 
			}
			if(foundifeat)
				break;
		}
	}
	return;
}

ia_t *gensplbdx(bgr_t2 *bed2, int m, int n, words_t *bedword, int m3, int n3) /* generate split bed index */
{
	/* TODO what you want is a copy of the data structure:
	 * NOPE! what you want is an array of indices */
	int i, k;
	ia_t *ia=calloc(1, sizeof(ia_t));
	ia->b=GBUF;
	ia->a=calloc(ia->b, sizeof(int));
	boole foundifeat;
	for(i=0;i<m;++i) {
		foundifeat=0;
		for(k=0;k<m3;++k) {
			if(!strcmp(bedword[k].n, bed2[i].f) ) {
				foundifeat=1;
				CONDREALLOC(ia->z, ia->b, GBUF, ia->a, int);
				ia->a[ia->z]=i;
				ia->z++;
			}
			if(foundifeat)
				break;
		}
	}
	ia->a=realloc(ia->a, ia->z*sizeof(int)); /*normalize */
	return ia;
}

void prtbed2f(bgr_t2 *bgrow, int m, int n, char *label) /* just print the feature names iof the f-bed (bed2) file */
{
	int i;
	printf("Feature file %s is %i rows by %i columns and is as follows:\n", label, m, n); 
	for(i=0;i<m;++i)
		printf("%s\n", bgrow[i].f);
	printf("\n"); 
	return;
}

void prtmbed(bgr_t **bgra, i4_t *dca, int dcasz, int n) /* the 2D version */
{
	int i, j;
	for(i=0;i<dcasz;++i) {
		for(j=0;j<dca[i].sc;++j) { // we're cycling though all of them, though we're really only interested in the first and last.
			if(j==0) { 
				printf("%s ", bgra[i][j].n);
				printf("%li ", bgra[i][j].c[0]);
			}
			if(j==dca[i].sc-1) { // note this cannot be an else if, because if only one line j = 0 = dca[i]-1.
				printf("%li ", bgra[i][j].c[1]);
				printf("%2.6f ", dca[i].mc);
			}
		}
		printf("\n"); 
	}
	return;
}

void m2beds(bgr_t *bgrow, bgr_t2 *bed2, int m2, int m) /* match up 2 beds */
{
	/* TODO: there could be an issue with intensity l;ines that span the end of one region and the start of another
	 * Need to look into that. this will only introduce a small error though.
	 */
	int i, j;
	int reghits; /* hits for region: number of lines in bed1 which coincide with a region in bed2 */
	int cloci; /* as opposed to hit, catch the number of loci */
	int rangecov=0;
	double assoctval=0;
	int istarthere=0, catchingi=0;
	boole caught;
	for(j=0;j<m2;++j) {
		caught=0;
		reghits=0;
		cloci=0;
		assoctval=0;
		for(i=istarthere;i<m;++i) {
			if( !(strcmp(bgrow[i].n, bed2[j].n)) & (bgrow[i].c[0] >= bed2[j].c[0]) & (bgrow[i].c[1] <= bed2[j].c[1]) ) {
				reghits++;
				rangecov=bgrow[i].c[1] - bgrow[i].c[0]; // range covered by this hit
				cloci+=rangecov;
				assoctval+=rangecov * bgrow[i].co;
				catchingi=i;
				caught=1;
			} else if (caught) { // will catch first untruth after a series of truths.
				caught=2;
				break; // bed1 is ordered so we can forget about trying to match anymore.
			}
		}
		if(caught==2)
			istarthere=catchingi+1;
		printf("Bed2idx %i / name %s / size %li got %i hits from bed1 , being %i loci and total assoc (prob .intensty) val %4.2f\n", j, bed2[j].f, bed2[j].c[1]-bed2[j].c[0], reghits, cloci, assoctval);
		if(istarthere >= m)
			break;
	}
	return;
}

void mbedp(dpf_t *dpf, bgr_t2 *bed2, int m2, int m) /* match up 2 beds */
{
	int i, j, min, max;
	int reghits; /* hits for region: number of lines in bed1 which coincide with a region in bed2 */
	int cloci; /* as opposed to hit, catch the number of loci */
	long assoctval=0;
	int istarthere=0, catchingi=0;
	boole caught;
	for(j=0;j<m2;++j) {
		caught=0;
		reghits=0;
		cloci=0;
		assoctval=0;
		min=9999999;
		max=0;
		for(i=istarthere;i<m;++i) {
			if( !(strcmp(dpf[i].n, bed2[j].n)) & (dpf[i].p >= bed2[j].c[0]) & (dpf[i].p < bed2[j].c[1]) ) {
				reghits++;
				cloci++;
				assoctval+= dpf[i].d;
				catchingi=i;
				if(dpf[i].d<min)
					min=dpf[i].d;
				if(dpf[i].d>max)
					max=dpf[i].d;
				caught=1;
			} else if (caught) { // will catch first untruth after a series of truths.
				caught=2;
				break; // bed1 is ordered so we can forget about trying to match anymore.
			}
		}
		if(caught==2)
			istarthere=catchingi+1;
#ifdef DBG
		printf("Bed2idx %i / name %s / size %li got %i hits from dpf , being %i loci and accumulated depth val of %lu\n", j, bed2[j].f, bed2[j].c[1]-bed2[j].c[0], reghits, cloci, assoctval);
#else
		printf("%s\t%li\t%li\t%s\t%i\t%i\t%li\t%4.4f\n", bed2[j].n, bed2[j].c[0], bed2[j].c[1], bed2[j].f, min, max, assoctval, (float)assoctval/cloci);

#endif
		if(istarthere >= m)
			break;
	}
	return;
}

i4_t *difca(bgr_t *bgrow, int m, int *dcasz, float minsig) /* An temmpt to merge bgraph quickly, no hope */
{
	int i, goodi=0 /* the last i at which minsig was satisfied */;
	boole seenminsig=0;
	/* how many different chromosomes are there? the dc (different chromsosome array */
	int dcbf=GBUF, dci=0;
	i4_t *dca=calloc(dcbf, sizeof(i4_t));
	char *tstr=NULL;
	/* find first bgrow element which is over the minimum coverage */
	for(i=0;i<m;++i)
		if(bgrow[i].co >= minsig) {
			tstr=malloc(bgrow[i].nsz*sizeof(char)); /* tmp string */
			strcpy(tstr, bgrow[i].n);
			dca[dci].sc++;
			dca[dci].mc=bgrow[i].co;
			dca[dci].b1i=i;
			dca[dci].lgbi=i;
			seenminsig=1;
			goodi=i;
			break;
		}
	if(!seenminsig) {
		printf("Error. No bedgraph element was able to satisfy the minimum signal value that was specified: abandoning ship.\n");
		exit(EXIT_FAILURE);
	}

	for(i=goodi+1;i<m;++i) {
		/* the same now means same name and contiguous */
		if( (!strcmp(tstr, bgrow[i].n)) & (bgrow[i].c[0] == bgrow[dca[dci].lgbi].c[1]) & (bgrow[i].co >= minsig) ) {
			dca[dci].sc++;
			dca[dci].lgbi=i;
			if(bgrow[i].co<dca[dci].mc)
				dca[dci].mc=bgrow[i].co;
		} else if (bgrow[i].co >= minsig) {
			CONDREALLOC(dci, dcbf, GBUF, dca, i4_t);
			dci++;
			dca[dci].sc++;
			dca[dci].mc=bgrow[i].co;
			dca[dci].b1i=i;
			dca[dci].lgbi=i;
			/* new string could be differnt length*/
			tstr=realloc(tstr, bgrow[i].nsz*sizeof(char)); /* tmp string */
			strcpy(tstr, bgrow[i].n);
		}
	}
	dca=realloc(dca, (dci+1)*sizeof(i4_t));
#ifdef DBG
	printf("Num of different chromcontigs=%i. How many of each? Let's see:\n", dci+1); 
	printf("dcbf=%i. 4-tupe is sc/mc/b1i/lgbi\n", dcbf); 
	for(i=0;i<=dci;++i) 
		printf("%i/%i/%i/%i ",dca[i].sc, dca[i].mc, dca[i].b1i, dca[i].lgbi); 
	printf("\n"); 
#endif
	*dcasz=dci+1;
	free(tstr);
	return dca;
}

void prtusage()
{
	printf("bgmergbl: this takes a bedgraph file, specified by -i, probably the bedgraph from a MACS2 intensity signal,\n");
	printf("and another bedgraph file, specified by -f, and merges the first into lines defined by the second.\n");
	printf("Before filtering however, please run with the -d (details) option. This will showi a rough spread of the values,\n");
	printf("so you can run a second time choosing filtering value (-f) more easily.\n");
	return;
}

int main(int argc, char *argv[])
{
	/* argument accounting */
	if(argc == 1) {
		prtusage();
		exit(EXIT_FAILURE);
	}
	int i, m, n, m2, n2, m3, n3, m4, n4;
	opt_t opts={0};
	catchopts(&opts, argc, argv);

	bgr_t *bgrow=NULL; /* usually macs signal */
	bgr_t2 *bed2=NULL; /* usually bed file from gff */
	words_t *bedword=NULL; /* usually feature names of interest */
	dpf_t *dpf=NULL; /* usually feature names of interest */
	if(opts.istr)
		bgrow=processinpf(opts.istr, &m, &n);
	if(opts.fstr)
		bed2=processinpf2(opts.fstr, &m2, &n2);
	if(opts.ustr)
		bedword=processwordf(opts.ustr, &m3, &n3);
	if(opts.pstr)
		dpf=processdpf(opts.pstr, &m4, &n4);

	if((opts.dflg) && (opts.istr)) {
		prtdets(bgrow, m, n, "Target bedgraph (1st) file");
		goto final;
	}
	if((opts.dflg) && (opts.fstr)) {
		prtbed2f(bed2, m2, n2, "Feature (bed2) file feature names");
		goto final;
	}
	// prtbed2(bed2, m2, MXCOL2VIEW);
	if((opts.istr) && (opts.fstr))
		m2beds(bgrow, bed2, m2, m);
	if((opts.ustr) && (opts.fstr) && (!opts.sflg)) {
		printf("bedwords:\n"); 
		for(i=0;i<m3;++i)
			printf("%s\n", bedword[i].n);
	}
	if((opts.pstr) && (opts.fstr) )
		mbedp(dpf, bed2, m2, m4);

	// if((opts.ustr) && (opts.fstr) && opts.sflg)
	// 	prtbed2s(bed2, m2, MXCOL2VIEW, bedword, m3, n3, "bed2 features that are in interesting-feature-file");

	ia_t *ia=NULL;
	if((opts.ustr) && (opts.fstr) && opts.sflg) {
		ia=gensplbdx(bed2, m2, n2, bedword, m3, n3);
		bed2in2(opts.fstr, bed2, m2, n2, ia);
	}

final:
	if(opts.pstr) {
		for(i=0;i<m4;++i)
			free(dpf[i].n);
		free(dpf);
	}
	if(opts.istr) {
		for(i=0;i<m;++i)
			free(bgrow[i].n);
		free(bgrow);
	}
	if(opts.fstr) {
		for(i=0;i<m2;++i) {
			free(bed2[i].n);
			free(bed2[i].f);
		}
		free(bed2);
	}
	if(opts.ustr) {
		for(i=0;i<m3;++i)
			free(bedword[i].n);
		free(bedword);
	}

	return 0;
}
