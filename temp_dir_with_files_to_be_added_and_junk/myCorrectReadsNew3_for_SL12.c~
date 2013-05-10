/* correct for read errors and put in the matrix result vector
 * input
 * % dbseqs - a list of compressed input sequence strings (1000000 x 4) of the matrix (uint32)
 * % invec - the compressed input (noisy) reads (XX x 4) (uint32)
 * % innum - the number of times each read in invec appears (XX x 1) (double)
 * % output
 * % resvec - the read error corrected (1000000) count vector
 * %
 * % note : 0,1,2,3 = A,C,G,T - Orzuk compression
 */

#include "stdio.h"
#include "math.h"
#include "mex.h"
#include "matrix.h"
#include "uthash.h"
#include <inttypes.h>

#define SEQLENGTH   8
#define MAXLISTSIZE 1000000

#define USE_64

struct SeqData {
    int seqpos;
    uint64_t seq[SEQLENGTH]; /* need to set the right length */
    
    UT_hash_handle hh; /* makes this structure hashable */
};

struct SeqData *hTable=NULL;

/* The gateway routine */
void mexFunction(int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[]) {
    
    mxArray *num_minutes, *str;
    /*unsigned int *dbseqs, *invec;*/
    uint64_t  *dbseqs, *invec;
    
    int numOfSeqs;
    int seqLen;
    int numReads;
    double *outDat;
    struct SeqData *cseqdat;
    struct SeqData *tseqdat,*testseqdat;
    int a, b;
    /*unsigned int creadseq[SEQLENGTH]; */
    /* unsigned int oldnuc; */
    uint64_t creadseq[SEQLENGTH];
    uint64_t oldnuc,newnuc,cnuc,oldcreadseq;
    uint64_t three;
    
    int cpos;
    int found0, found1;
    int *ppos;
    int clpos;
    unsigned char tmpbuf[255];
    unsigned char tmpstr[255];
    /* indices for non unique database kmers*/
    int *nonUnique;
    
    double totfreq;
    struct seqData *current_user, *tmp;
    int cpack;
    double alpha, p;
    int ccpos;
    double *innum;

    /* use three to maybe prevent uint64 problem */
    three=3;

    /* prevent trash in kmers shorter that MAXSEQLEN */
    for (b=0;b<SEQLENGTH;b++) {
        creadseq[b]=0;
    }
    
    dbseqs = (uint64_t *)  mxGetPr(prhs[0]);
    
    numOfSeqs = mxGetM(prhs[0]);
    seqLen = mxGetN(prhs[0]);
    if (seqLen!=SEQLENGTH) {
        printf("Different seq lens from #define!!!\n");
        #ifndef USE_64
                return;
        #endif
                /*seqLen=SEQLENGTH;*/
    }
    invec = (uint64_t *)  mxGetPr(prhs[1]);
    
    numReads = mxGetM(prhs[1]);
    if (mxGetN(prhs[1])!=seqLen) {
        printf("Different seq lens!!!\n");
        #ifndef USE_64
                return;
        #endif
    }
    
    plhs[0] = mxCreateDoubleMatrix(numOfSeqs, 1, mxREAL);
    outDat=mxGetPr(plhs[0]);
    
    if ((tseqdat=(struct SeqData *) calloc(numOfSeqs, sizeof(struct SeqData)))==NULL) {
        printf("Can't allocate seqdata");
        return;
    }
    
    innum = (double *)  mxGetPr(prhs[2]);
    if (mxGetM(prhs[2])!=numReads) {
        printf("innum not same size as number of reads!!!\n");
        return;
    }
    if (mxGetN(prhs[2])!=1) {
        printf("innunm should be a vector and not a matrix!!!\n");
        return;
    }
    
    if ((nonUnique=(int *) calloc(numOfSeqs, sizeof(int)))==0) {
        printf("Can't allocate nonUnique\n");
        return;
    }
    

    /* add database sequences to hash */
    for (a=0;a<numOfSeqs;a++) {
        outDat[a]=0;
        cseqdat=&(tseqdat[a]);
        cseqdat->seqpos=a;
        
        for (b=0;b<seqLen;b++) {
            cseqdat->seq[b]=dbseqs[a+b*numOfSeqs];
            creadseq[b]=dbseqs[a+b*numOfSeqs];
        }
        /* look if kmer already in the database (non-unique) */
        HASH_FIND( hh, hTable, creadseq, SEQLENGTH*sizeof(uint64_t), testseqdat);
        if (testseqdat) {
            printf("non unique %d already in %d\n",a,testseqdat->seqpos);
            /* k-mer non unique - store it's similar k-mer position */
            nonUnique[a]=testseqdat->seqpos;
        } else 
        {
            /* k-mer is unique so add to hash */
            HASH_ADD( hh, hTable, seq, SEQLENGTH*sizeof(uint64_t), cseqdat );
            nonUnique[a]=-1;
        }
    }
    
    
    if ((ppos=(int *) calloc(numOfSeqs, sizeof(int)))==0) {
        printf("Can't allocate ppos\n");
        return;
    }
    

    /* now go over all the reads */
    for (a=0;a<numReads;a++) {
        for (b=0;b<seqLen;b++) {
            creadseq[b]=invec[a+b*numReads];
        }
        found0=0;
        found1=0;
        clpos=0;
        
        /* if original sequence is found, don't modify frequencies of other sequences */
        HASH_FIND(hh, hTable, creadseq, SEQLENGTH*sizeof(uint64_t), cseqdat);
        if (cseqdat) {
/*            printf("found orig %d position %d\n",a,cseqdat->seqpos);*/
            found0++;
            /* add the observed number of reads */
            outDat[cseqdat->seqpos]+=innum[a];
        } else {
            /* original sequence is not found - find all neighbours */
            ccpos=sizeof(uint64_t);
            ccpos=0;
            for (cpack=0;cpack<seqLen;cpack++) {
                oldcreadseq=creadseq[cpack];
                for (cpos=0;cpos<sizeof(uint64_t)*4;cpos++) {
                    ccpos++;
                    oldnuc=creadseq[cpack];
                    oldnuc=(oldnuc>>(cpos*2));
                    newnuc=oldnuc&three;
                    for (cnuc=0;cnuc<4;cnuc++) {
                        if (cnuc!=newnuc) {
                            creadseq[cpack]=(creadseq[cpack]&~(three<<(cpos*2)))|(cnuc<<cpos*2);
                            HASH_FIND( hh, hTable, creadseq, SEQLENGTH*sizeof(uint64_t), cseqdat);
                            if (cseqdat) {
                                ppos[clpos]=cseqdat->seqpos;
                                clpos++;
                                if (clpos>numOfSeqs) {
                                    printf("Error clpos=%d numofseqs=%d\n",clpos,numOfSeqs);
                                    clpos--;
                                }
                            }
                        }
                    }
                    creadseq[cpack]=oldcreadseq;
                }
            }
            /* add 1 to frequency of all neighboring seqs found */
            if (clpos>0) {
                for (b=0;b<clpos;b++) {
                    outDat[ppos[b]]+=innum[a];
                }
            }
        }
    }
    /* now copy results for non-unique sequences*/
    for (a=0;a<numOfSeqs;a++) {
        if (nonUnique[a]>-1) {
/*            printf("list non unique %d is similar to %d\n",a,nonUnique[a]);*/
            outDat[a]=outDat[nonUnique[a]];
        }
    }

    free(ppos);
    free(nonUnique);
    
    HASH_CLEAR(hh, hTable);
    free(tseqdat);
}
