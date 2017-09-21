/****************************************************************************
**  Copyright (c) 1994-2017 ALTAIR Engineering, Inc.
**  All rights reserved.
**  This source code is confidential and may not be disclosed.
****************************************************************************/

/*===========================================================================
**
** "ste.h":  Split Tree Element
**
** Original: Farzin Shakib (Jan 98)
**===========================================================================
*/

#ifndef	__STE_H__
#define	__STE_H__

/*===========================================================================
 *
 * Include needed files
 *
 *===========================================================================
 */

#define STE_B_BND_REALS		6             
#define STE_E_BND_REALS		18

typedef struct _Bbnd {
    Real	box[STE_B_BND_REALS] ;	/* basic B-Box in xyz		*/
} Bbnd ;
typedef Bbnd*	BbndHd ;

typedef struct _Ebnd {
    Real        box[STE_E_BND_REALS] ;	/* enhanced B-Boxes		*/
} Ebnd ;
typedef Ebnd*   EbndHd ;

/*===========================================================================
 *
 * "SteBox":  Ste Binary Tree Set
 *
 *===========================================================================
 */
typedef struct _SteBox* SteBoxHd ;

typedef struct _SteBox {
    Integer	nRefs ;			/* No. references		*/
    Integer*	refs ;			/* references of this set	*/
    Real*	box  ;                  /* bnd boxes of this set	*/
    SteBoxHd	child1 ;		/* child 1			*/
    SteBoxHd	child2 ;		/* child 2			*/
} SteBox ;
    
/*===========================================================================
 *
 * "Ste":  Ste main structure
 *
 *===========================================================================
 */
typedef struct _Ste {
    Real*	crd ;			/* coordinates			*/
    Integer*	cnn ;			/* connectivity			*/
    Integer	nElems ;		/* No. elements			*/
    Integer	tplType ;		/* topology type		*/
    Integer	nElemNodes ;		/* No. element nodes		*/
    Integer	nCuts ;			/* No. elements cut		*/
    Integer*	refs ;			/* reference			*/
    Real*	elmBnds ;		/* bounds for elements		*/
    Bool	cpuOpt ;		/* option for CPU optimization	*/
    SteBoxHd	rootHd ;		/* root				*/
    MemBlkHd	memHd ;			/* memory handle for boxHd	*/
    MemBlkHd	mbsHd ;			/* memory handle of boxHd->box	*/
    MtxHeapHd	heapHd ;		/* heap for tree traversing	*/
    Real*	wght ;			/* heap keys ( distances )	*/
    SteBoxHd*	boxHds ;                /* array of SteBoxHds		*/
    Integer	nBdrs ;			/* 6 or 18, No. real numbers	*/
					/* defining 1 or 3 bnd boxes	*/
} Ste ;

typedef struct _Ste* SteHd ;

/*===========================================================================
 *
 * Fortran names
 *
 *===========================================================================
 */
#define fSteGetDistBrick fstegetdistbrick_
#define fSteGetDistTet fstegetdisttet_
#define fSteGetDistPyramid fstegetdistpyramid_
#define fSteGetDistWedge fstegetdistwedge_
#define fSteSplitBnds fstesplitbnds_
#define fSteElmBnds fsteelmbnds_
#define fSteUpdBoxBnds fsteupdboxbnds_
#define fSteUpdBoxEBnds fsteupdboxebnds_
#define fSteCutPlaneBlk fstecutplaneblk_
#define fStePntRotX fstepntrotx_
#define fStePntRotY fstepntroty_
#define fStePntRotZ fstepntrotz_

/*===========================================================================
 *
 * Function prototype
 *
 *===========================================================================
 */
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
SteBoxHd	steSplitTree(		SteHd		steHd,
					Integer*	refs,
					Integer		nRefs		) ;
Void		steGetDist(		SteHd		steHd,
					SteBoxHd	steBoxHd,
					Real*		pCrd,
					Real*		lCrd,
					Integer*	id,
					Real*		dist		) ;
Integer		steCutPlane(		SteHd		steHd,
					Real*		cpCrds, 
					Integer*	tags		) ;
Void		steCutPlaneBox(		SteHd		steHd,
					SteBoxHd	boxHd,
					Integer*	tags,
					Real*		cpPnts,
					Real*		cpBnd,
					Real*		cpDirs,
					Real*		cpBndRot,
					Real		lxRot,
					Real		lyRot,
					Integer*	count		) ;
Real		steCalcBoxDist(		Real*           bnds,
                               		Real*           pCrd,
					Bool		cpuOpt          ) ;
Integer		steFindElemNextHeap(    SteHd           steHd ,
                        		Real*           pCrd ,
                        		Real*           lCrd ,
                        		Real*           dist     	) ;
/*===========================================================================
 *
 * Fortran functions
 *
 *===========================================================================
 */
Void		fStePntRotX( 		Real*		cpDirs, 
					Real*		orig, 
					Real*		pnt,
					Real*		rotPnt 		) ;
Void            fStePntRotY(            Real*           cpDirs,
                                        Real*           orig,
                                        Real*           pnt,
                                        Real*           rotPnt          ) ;
Void            fStePntRotZ(            Real*           cpDirs,
                                        Real*           orig,
                                        Real*           pnt,
                                        Real*           rotPnt          ) ;
Void		fSteElmBnds(		Real*		elmBnds,
					Real*		crd,
					Integer*	cnn,
                                        Integer*        nElemNodes, 
					Integer*	nRefs,
					Integer*	nBdrs		) ;
Void		fSteUpdBoxBnds(	        Real*		box,
					Real*		elmBnds,
					Integer*	refs,
					Integer*	nRefs,
					Integer*	nBdrs		) ;
Void            fSteUpdBoxEBnds(        Real*           box,
                                        Real*           elmBnds,
                                        Integer*        refs,
                                        Integer*        nRefs,
                                        Integer*        nBdrs           ) ;
Void		fSteSplitBnds(	        Real*           elmBnds,	
					Integer*	refs,
					Integer*	nRefs,
					Integer*	nRef1s,
					Integer*	nBdrs		) ;
Void		fSteGetDistTet(		Real*		pCrd,
					Real*		crd,
					Integer*	cnn,
					Integer*	refs,
					Real*		lCrd,
					Integer*	nRefs,
					Integer*	id,
					Real*		dist		) ;
Void		fSteGetDistWedge(	Real*		pCrd,
					Real*		crd,
					Integer*	cnn,
					Integer*	refs,
					Real*		lCrd,
					Integer*	nRefs,
					Integer*	id,
					Real*		dist		) ;
Void		fSteGetDistPyramid(	Real*		pCrd,
					Real*		crd,
					Integer*	cnn,
					Integer*	refs,
					Real*		lCrd,
					Integer*	nRefs,
					Integer*	id,
					Real*		dist		) ;
Void		fSteGetDistBrick(	Real*		pCrd,
					Real*		crd,
					Integer*	cnn,
					Integer*	refs,
					Real*		lCrd,
					Integer*	nRefs,
					Integer*	id,
					Real*		dist		) ;
Void		fSteCutPlaneBlk(	Real*		cpDirs, 	
					Real*		cpBndRot,	
					Integer*	tags,		
					Real*		crd,	
					Integer*	cnn,	
					Integer*	refs,	
					Integer*	nRefs,	
					Integer*	nElemNodes,		
					Integer*	nCuts		) ;

/*===========================================================================
 *
 * End of the file
 *
 *===========================================================================
 */
#endif	/* __STE_H__ */
