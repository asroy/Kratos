/****************************************************************************
**  Copyright (c) 1994-2017 ALTAIR Engineering, Inc.
**  All rights reserved.
**  This source code is confidential and may not be disclosed.
****************************************************************************/

/*===========================================================================
**
** "steDim.h":  Split Tree Element
**
** Original: Farzin Shakib (Jan 98)
**===========================================================================
*/

#ifndef	__STE_DIM_H__
#define	__STE_DIM_H__

/*===========================================================================
 *
 * Macros
 *
 *===========================================================================
 */
#define	STE_MIN_SIZE	16		/* size of the leaf nodes	*/
/*===========================================================================
 *
 * The STE_HEAP_SIZE is approximate, and for large models, this number may
 * not be sufficient.   Testing to date has indicated this is OK, but
 * if there are problems, simply increase STE_HEAP_SIZE.
 *
 *===========================================================================
 */
#define	STE_HEAP_SIZE   131072          /* size of the heap array 	*/
/*===========================================================================
 *
 * End of the file
 *
 *===========================================================================
 */
#endif	/* __STE_DIM_H__ */
