#include <ctype.h>
#include <dirent.h>
#include <errno.h>
#include <fcntl.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <pwd.h>
#include <setjmp.h>
#include <signal.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>


#define PRM_TPL_3TRI                     1
#define PRM_TPL_4QUAD                    2
#define PRM_TPL_6TRI                    15
#define PRM_TPL_8QUAD                   31
#define PRM_TPL_9QUAD                   32 

#define PRM_TPL_4TET                     3
#define PRM_TPL_5TET                    28
#define PRM_TPL_5PYRAMID                 4
#define PRM_TPL_6WEDGE                   5
#define PRM_TPL_8BRICK                   6
#define PRM_TPL_10TET                   16
#define PRM_TPL_13PYRAMID               33
#define PRM_TPL_14PYRAMID               34
#define PRM_TPL_15WEDGE                 35
#define PRM_TPL_18WEDGE                 36
#define PRM_TPL_20BRICK                 37
#define PRM_TPL_27BRICK                 38 

#define PRM_TPL_4TET_BND                 7
#define PRM_TPL_5TET_BND                29
#define PRM_TPL_5PYRAMID_BND_P           8
#define PRM_TPL_5PYRAMID_BND_Q           9
#define PRM_TPL_6WEDGE_BND_P            10
#define PRM_TPL_6WEDGE_BND_Q            11
#define PRM_TPL_8BRICK_BND              12
#define PRM_TPL_10TET_BND               17

#define PRM_TPL_MIXED                   13
#define PRM_TPL_MIXED_BND               14

#define PRM_TPL_4TET_EDGE               18
#define PRM_TPL_5TET_EDGE               30
#define PRM_TPL_5PYRAMID_EDGE_P1        19
#define PRM_TPL_5PYRAMID_EDGE_P2        20
#define PRM_TPL_5PYRAMID_EDGE_P3        21
#define PRM_TPL_5PYRAMID_EDGE_Q         22
#define PRM_TPL_6WEDGE_EDGE_P           23
#define PRM_TPL_6WEDGE_EDGE_Q1          24
#define PRM_TPL_6WEDGE_EDGE_Q2          25
#define PRM_TPL_8BRICK_EDGE             26
#define PRM_TPL_10TET_EDGE              27

typedef int     Integer ;               /* integer type                 */
typedef double  Real ;                  /* real    type                 */
typedef float   Float ;                 /* float   type                 */
typedef char*   String ;                /* string  type                 */
typedef void*   Data ;                  /* data    type                 */
typedef void    Void ;                  /* void    type                 */
typedef size_t  Size ;                  /* size    type                 */
typedef long long       Integer64 ;     /* integer type 64 bit          */
typedef unsigned int    UInteger ;      /* unsigned integer type        */

#ifdef TRUE
#undef TRUE
#endif
#ifdef FALSE
#undef FALSE
#endif

typedef enum {
    FALSE = 0 ,
    TRUE  = 1
} Bool ;                                /* boolean type                 */

#define newMem(T,n)     malloc( sizeof(T)*n ) 

/*===========================================================================
 *
 * "MtxHeap":  Heap data structure
 *
 *===========================================================================
 */
typedef struct _MtxHeap {
    Integer*	heap ;			/* heap pointer			*/
    Integer*	invMap ;		/* heap inverse map		*/
    Real*	wght ;			/* heap values			*/
    Integer	endId ;                 /* size of wght                 */
    Integer	maxFlag ;		/* sort max to min?		*/
    Integer	nItems ;		/* No. items in the heap	*/
    Integer	maxItems ;              /* length of the heap array     */
} MtxHeap ;

typedef struct _MtxHeap* MtxHeapHd ;

/*===========================================================================
 *
 * Useful heap macros
 *
 *===========================================================================
 */
#define	mtxHeapDim(H)		(H)->nItems
#define	mtxHeapHead(H)		(H)->heap[0]
#define	mtxIndxHeap(H,I)	(H)->invMap[I]
#define	mtxPopHeapItem(H,I)	mtxPopHeap(H,mtxIndxHeap(H,I))
#define	mtxAdjHeapItem(H,I)	mtxAdjHeap(H,mtxIndxHeap(H,I))
#define mtxSwapInt(A,B)         { Integer I ; I = A ; A = B ; B = I ; }
#define mtxSwapReal(A,B)        { Real    I ; I = A ; A = B ; B = I ; }
#define mtxDot2Vec3(A,B)        ((A)[0]*(B)[0]+(A)[1]*(B)[1]+(A)[2]*(B)[2])
#define mtxDaxpyVec3(A,B,S)                                             \
                ((B)[0]+=(S)*(A)[0],(B)[1]+=(S)*(A)[1],(B)[2]+=(S)*(A)[2])
#define mtxScaleVec3(A,S)       ((A)[0]*=(S),(A)[1]*=(S),(A)[2]*=(S))


Void		mtxAddHeap( 		MtxHeapHd       heapHd,
            				Real            value     	);
MtxHeapHd	mtxNewHeap(		Real*		wght,	
					Integer		maxFlag,
					Integer		nItems,
					Integer		maxItems	) ;
Void		mtxInitHeap(		MtxHeapHd	heapHd,
					Integer*	heap,
					Integer*	invMap,
					Real*		wght,	
					Integer		maxFlag,
					Integer		nItems		) ;
Integer		mtxPopHeap(		MtxHeapHd	heapHd,	
					Integer		index		) ;
Void		mtxAdjHeap(		MtxHeapHd	heapHd,	
					Integer		index		) ;
Void		mtxFreeHeap(		MtxHeapHd	heapHd		) ;

Void            mtxOrth1Vec3(           Real*           src,
                                        Real*           dst             ) ;
Void            mtxCross(               Real*           src1,
                                        Real*           src2,
                                        Real*           dst             ) ;

Void            mtxZeroInt(             Integer*        a,
                                        Integer         n               ) ;

#define MTH_MAX_INT             INT_MAX
#define MTH_INT_MAX             INT_MAX

#define MTH_MAX_INT64           LONG_MAX
#define MTH_INT64_MAX           LONG_MAX

#define MTH_MAX_REAL            DBL_MAX
#define MTH_REAL_MAX            DBL_MAX

#define MTH_SQRT_REAL_MAX       sqrt((double)DBL_MAX)
#define MTH_SQRT_MAX_REAL       sqrt((double)DBL_MAX)
#define MTH_FTRT_REAL_MAX       sqrt((double)MTH_SQRT_REAL_MAX)
#define MTH_FTRT_MAX_REAL       sqrt((double)MTH_SQRT_MAX_REAL)

#define MTH_EPS                 DBL_EPSILON
#define MTH_SQRT_EPS            sqrt((double)DBL_EPSILON)

#define MTH_PI                  3.1415926535897931E+00
#define MTH_E                   2.7182818284590452E+00
#define MTH_SQRT_2              1.4142135623730951E+00
#define MTH_INV_SQRT_2          7.0710678118654752E-01
#define MTH_SQRT_3              1.7320508075688773E+00
#define MTH_INV_SQRT_3          5.7735026918962577E-01

#define MTH_1K                  1024
#define MTH_1MEG                1048576

#define mthMax(A,B)             (((A) > (B)) ? (A) :  (B))
#define mthMin(A,B)             (((A) < (B)) ? (A) :  (B))
#define mthAbs(A)               (((A) > 0  ) ? (A) : -(A))
#define mthPow2(A)              ((A)*(A))

#define mthAcos(A)              acos((double) (A))
#define mthAsin(A)              asin((double) (A))
#define mthAtan(A)              atan((double) (A))
#define mthAtan2(A,B)           atan2((double) (A), (double) (B))
#define mthCeil(A)              ceil((double) (A))
#define mthCos(A)               cos((double) (A))
#define mthCosh(A)              cosh((double) (A))
#define mthExp(A)               exp((double) (A))
#define mthFloor(A)             floor((double) (A))
#define mthFmod(A,B)            fmod((double) (A), (double) (B))
#define mthLog(A)               log((double) (A))
#define mthLog10(A)             log10((double) (A))
#define mthLog2(A)              (log((double) (A))/log((double)2))
#define mthMod(A,B)             fmod((double) (A), (double) (B))
#define mthPow(A,B)             pow((double) (A), (double) (B))
#define mthSin(A)               sin((double) (A))
#define mthSinh(A)              sinh((double) (A))
#define mthSqrt(A)              sqrt((double) (A))



typedef struct _MemBlk*		MemBlkHd ;
typedef struct _MemFree*	MemFreeHd ;

/*===========================================================================
 *
 * "memFree": free structure
 *
 *===========================================================================
 */
typedef struct _MemFree {
    MemFreeHd	next ;			/* next free item		*/
} MemFree ;

/*===========================================================================
 *
 * "memBlk": data structure
 *
 *===========================================================================
 */
typedef struct _MemBlk {
    MemFreeHd	free ;			/* free list			*/
    Data	aTmp ;			/* temporary location		*/
    Data*	bTmp ;			/* temporary location		*/
    size_t	ind ;			/* index to the next item	*/
    size_t	itemTop ;		/* top element of the blk	*/
    Integer	itemDim ;		/* dimension of items		*/
    Integer	itemSize ;		/* size of the items		*/
    size_t	blkSize ;		/* size of the blk		*/
    size_t	blkDim ;		/* dimension of the blk		*/
    Data*	data ;			/* array of data		*/
} MemBlk ;

/*===========================================================================
 *
 * "memBlkNew":
 *
 *===========================================================================
 */
#define	memBlkNew(T,S)		memBlkNewF(sizeof(T),S)

/*===========================================================================
 *
 * "memBlkGet": Get a new items from block allocation
 *
 *===========================================================================
 */
#define	memBlkGet(T,B)							\
    (T*)((B)->free							\
	?( (B)->aTmp=(B)->free,						\
	   (B)->free=(B)->free->next,					\
	   (B)->aTmp )							\
	:(B)->ind 							\
	?( (B)->ind-=(B)->itemDim,					\
	   (B)->data+(B)->ind )						\
	:( (B)->bTmp=newMem(Data,(B)->blkDim),			        \
	   (B)->bTmp[(B)->blkDim-1]=(B)->data,				\
	   (B)->data=(B)->bTmp,						\
	   (B)->ind=(B)->itemTop,					\
	   (B)->data+(B)->ind ) )

/*===========================================================================
 *
 * "memBlkDel":  Give back an item to the memory block
 *
 *===========================================================================
 */
#define	memBlkDel(B,D)							\
    ((MemFreeHd)(D))->next = (B)->free ; (B)->free = (MemFreeHd)(D)


MemBlkHd        memBlkNewF(             Integer         itemSize,
                                        size_t         blkSize         ) ;
Void            memBlkFree(             MemBlkHd        memBlkHd        ) ;



#define fMtxCrdBnd fmtxcrdbnd_


Void            fMtxCrdBnd(             Real*           crd,
                                        Real*           bnd,
                                        Integer*        dim             ) ;

