#include <cstdlib>
#include <cstdio>

#define SYS_BUFF_SIZE                   256
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


extern "C" {

struct _Ste ;
typedef struct _Ste* SteHd ;


SteHd		steNew(			Real*		crd,
					Integer*	cnn,
					Integer		tplType,
					Integer		nElems,		
					Bool		cpuOpt		) ;
Void		steFree(		SteHd		steHd		) ;
Integer		steFindElem(		SteHd		steHd,
					Real*		pCrd,
					Real*		lCrd,
					Real*		dist		) ;
Integer		steFindElemNext(	SteHd		steHd,
					Real*		pCrd,
					Real*		lCrd,
					Real*		dist		) ;
Integer		steFindElemNextHeap(    SteHd           steHd ,
                        		Real*           pCrd ,
                        		Real*           lCrd ,
                        		Real*           dist     	) ;

}
