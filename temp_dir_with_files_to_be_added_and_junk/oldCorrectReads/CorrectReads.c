/* correct for read errors and put in the matrix result vector
 * input
 * % dbseqs - a list of compressed input sequence strings (1000000 x 4) of the matrix (uint32)
 * % invec - the compressed input (noisy) reads (? x 4) (uint32)
 * % errorMat - the original ZZ model error matrix
 * % E0 - the error prob. at position 0
 * % E50 - the error prob. at position 50
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

#define SEQLENGTH   4
#define MAXLISTSIZE 1000000

// comment this if u pass the data packed as uint32 table
// uncomment if the data is packed as uint64
#define USE_64

struct SeqData {
    int seqpos;
    unsigned int seq[SEQLENGTH]; /* need to set the right length */
    UT_hash_handle hh; /* makes this structure hashable */
};

struct SeqData *hTable=NULL;

/* The gateway routine */
void mexFunction(int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[]) {
    
    mxArray *num_minutes, *str;
    unsigned int *dbseqs, *invec;
    int numOfSeqs;
    int seqLen;
    int numReads;
    double *outDat;
    struct SeqData *cseqdat;
    struct SeqData *tseqdat;
    int a, b;
    unsigned int creadseq[SEQLENGTH];
    unsigned int oldnuc;
    int cpos;
    int cnuc;
    int found0, found1;
    int *ppos;
    double *pfreq;
    int clpos;
    unsigned char tmpbuf[255];
    unsigned char tmpstr[255];
    int tmpseq[SEQLENGTH];
    double totfreq;
    struct seqData *current_user, *tmp;
    int cpack;
    double alpha,p;
    int ccpos;
    double *errormat;
    double *E0,*E50;
    
        
    dbseqs = (unsigned int *)  mxGetPr(prhs[0]);
    numOfSeqs = mxGetM(prhs[0]);
    seqLen = mxGetN(prhs[0]);
    if (seqLen!=SEQLENGTH) {
        printf("Different seq lens from #define!!!\n");
 #ifndef USE_64
         return;
 #endif
        seqLen=SEQLENGTH;
    }
    invec = (unsigned int *)  mxGetPr(prhs[1]);
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
    
    errormat = (double *)  mxGetPr(prhs[2]);
    if (mxGetN(prhs[2])!=4) {
        printf("Bad error mat!!!\n");
        return;
    }
    if (mxGetM(prhs[2])!=4) {
        printf("Bad error mat!!!\n");
        return;
    }
    E0 = (double *)  mxGetPr(prhs[3]);
    E50 = (double *)  mxGetPr(prhs[4]);
//    printf("E0-%f E50-%f\n",*E0,*E50);
    alpha=log((*E50)/(*E0))/49;
    
  
    for (a=0;a<numOfSeqs;a++) {
        for (b=0;b<seqLen;b++) {
            tmpseq[b]=dbseqs[a+b*numOfSeqs];
        }
        /*        if ((cseqdat=(struct SeqData *)malloc(sizeof(struct SeqData)))==NULL) {
            printf("Out of mem\n");
            return;
        }*/
        cseqdat=&(tseqdat[a]);
        cseqdat->seqpos=a;
 //       printf("Added seq %d\n",a);
        for (b=0;b<seqLen;b++) {
            #ifndef USE_64
            cseqdat->seq[b]=dbseqs[a+b*numOfSeqs];
            #endif
            #ifdef USE_64
            cseqdat->seq[b]=dbseqs[a*2+(b/2)*numOfSeqs*2+(b%2)];
            #endif
            for (cpos=0;cpos<sizeof(int)*4;cpos++) {
                oldnuc=(cseqdat->seq[b]&(3<<(cpos*2)))>>(cpos*2);
//                printf("%d",oldnuc);
            }
        }
//        printf("\n");
        HASH_ADD( hh, hTable, seq, SEQLENGTH*sizeof(int), cseqdat );
    }
    
    if ((ppos=(int *) calloc(numOfSeqs, sizeof(int)))==0) {
        printf("Can't allocate ppos\n");
        return;
    }
    if ((pfreq=(double *) calloc(numOfSeqs, sizeof(double)))==0) {
        printf("Can't allocate pfreq\n");
        return;
    }
    
    for (a=0;a<numReads;a++) {
//            printf("\n seq %d:\n",a);
        for (b=0;b<seqLen;b++) {
            #ifndef USE_64
            creadseq[b]=invec[a+b*numReads];
            #endif
            #ifdef USE_64
//                    printf("%d\n",a*2+(b/2)*numReads*2+(b%2));
            creadseq[b]=invec[a*2+(b/2)*numReads*2+(b%2)];
            #endif
//            creadseq[b]=invec[a+b*numReads];
        }
        found0=0;
        found1=0;
        clpos=0;
//         HASH_FIND(hh, hTable, creadseq, SEQLENGTH*sizeof(int), cseqdat);
//         if (cseqdat) {
//             found0++;
//             pfreq[clpos]=1;
//             ppos[clpos]=cseqdat->seqpos;
//             clpos++;
//         }
        ccpos=0;
        for (cpack=0;cpack<seqLen;cpack++) {
            for (cpos=0;cpos<sizeof(int)*4;cpos++) {
                ccpos++;
                oldnuc=(creadseq[cpack]&(3<<(cpos*2)))>>(cpos*2);
//                printf("%d",oldnuc);
                for (cnuc=0;cnuc<4;cnuc++) {
                    if (cnuc!=oldnuc) {
                        creadseq[cpack]=(creadseq[cpack]&~(3<<(cpos*2)))|(cnuc<<cpos*2);
                        HASH_FIND( hh, hTable, creadseq, SEQLENGTH*sizeof(int), cseqdat);
                        if (cseqdat) {
                            found1++;
//                            printf("ccpos=%d alpha=%f\n",ccpos,alpha);
                            p=(*E0)*exp(alpha*(ccpos-1));
                            pfreq[clpos]=errormat[cnuc+oldnuc*4]*p;
                            ppos[clpos]=cseqdat->seqpos;
                            clpos++;
                        }
                    }
                }
                creadseq[cpack]=(creadseq[cpack]&~(3<<(cpos*2)))|(oldnuc<<cpos*2);
            }
        }
        totfreq=0;
        if (clpos>0) {
            for (b=0;b<clpos;b++) {
//                printf("similar - %d %f\n",ppos[b],pfreq[b]);
                totfreq+=pfreq[b];
            }
            if (totfreq>0) {
                for (b=0;b<clpos;b++) {
//                    outDat[ppos[b]]+=pfreq[b]/totfreq;
                    outDat[ppos[b]]+=pfreq[b];
                }
            } else {
                printf("Error-freqsum=0\n");
            }
        }
        HASH_FIND(hh, hTable, creadseq, SEQLENGTH*sizeof(int), cseqdat);
        if (cseqdat) {
            found0++;
            outDat[cseqdat->seqpos]+=1-totfreq;
//            printf("\n found myself\n");
//            pfreq[clpos]=1-totfreq;
//            ppos[clpos]=cseqdat->seqpos;
//            clpos++;
        } 
//        else {
//            printf("I don't exist!\n");
//        }
            

    }
    free(ppos);
    free(pfreq);
    
    HASH_CLEAR(hh, hTable);
    free(tseqdat);
}
