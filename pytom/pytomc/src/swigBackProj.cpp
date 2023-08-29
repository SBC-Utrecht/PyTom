/*
 * swigBackProj.cpp
 *
 *  Created on: Dec 7, 2010
 *      Author: hrabe
 */

/*=================================================================
 *
 * tom_backproj3d.c	Performs a 3D Backprojection
 *
 *
 * The syntax is:
 *
 *		tom_backproj3d(IN,PHI,PSI,DIMENSIONS[3],OFFSET[3])
 *
 *  Electron Tomography toolbox of the
 *  Max-Planck-Institute for Biochemistry
 *  Dept. Molecular Structural Biology
 *  82152 Martinsried, Germany
 *  http://www.biochem.mpg.de
 *
 * Last changes: 04.07.2002
 * By: Stephan Nickell
 * Revision: 1.00 by
 *
 *=================================================================*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mex.h"


//#include <swigVolume.hpp>
//#include <swigBackProj.hpp>


#define PI 3.141592653589793238

/* Input Arguments */
#define	IN	prhs[1]
#define	PHI	prhs[2]
#define	THE	prhs[3]
#define	OFFS	prhs[4]
#define	VOL	prhs[0]



static void tom_backproj3d(float *O,float *I,size_t Ox_max,size_t Oy_max,size_t Oz_max,
		float Phi,float The,size_t Ix_max,size_t Iy_max,
		int x_offset,int y_offset,int z_offset
)
{
	long int DRX, DRY,DRZ;
	float XE,YE,ZE;
	int DPX,DPY;
	float XEP,YEP;
	float T[4][4];
	int INDR;
	int IZO,IZ,IZM;
	float RADD;
	long int IY,IX;
	float RX,RY,XD,YD;
	long int INDA;
	float XIN1,XIN2,YIN;
	int DXR,DYR,DZR,DXP,DYP,LL,MM;
	float RXADD, RYADD;

	INDA=0;
	/*  OBJEKTDATEN:  */
	DRX=Ox_max;
	DRY=Oy_max;
	DRZ=Oz_max;
	IZO=1;
	IZM=DRZ;
	XE = 1.0;
	YE = 1.0;
	ZE = 1.0;
	/*  PROJEKTIONSDATEN: */
	DPX = Ix_max;
	DPY = Iy_max;
	XEP = 1.0;
	YEP = 1.0;

	/*  Prepare rotation matrix */
	Phi=(float)Phi*(float)PI/(float)180.0;    /*radians*/
	The=(float)The*(float)PI/(float)180.0;    /*radians*/
	T[1][1] = (float)cos(The)*(float)cos(Phi)*(float)XE/(float)XEP;
	T[2][1] = (float)cos(The)*(float)sin(Phi)*(float)YE/(float)XEP;
	T[3][1] = -(float)sin(The)*(float)ZE/(float)XEP;
	T[1][2] = -(float)sin(Phi)*(float)XE/(float)YEP;
	T[2][2] = (float)cos(Phi)*(float)YE/(float)YEP;
	/*C  Define centre of the reconstruction body*/
	DXR=DRX/2+1 - x_offset;
	DYR=DRY/2+1 - y_offset;
	DZR=DRZ/2+1 - z_offset;
	DXP=DPX/2+1;
	DYP=DPY/2+1;

	/*C
	C  Main loop over all pixels of the object*/
	INDR=0;
	INDA=0;
	XIN1=0.0;
	XIN2=0.0;
	YIN=0.0;

	for (IZ=1;IZ<=IZM;IZ++){
		RADD = (IZ-DZR)*T[3][1];

		for (IY=1;IY<=DRY;IY++){
			RXADD = (IY-DYR)*T[2][1] + RADD;
			RYADD = (IY-DYR)*T[2][2];

			for (IX=1;IX<=DRX;IX++){

				/*C  Calculate point (RX,RY) in the projection plane corresponding to
				C  the object point (IX,IY,IZ)*/
				RX = (IX-DXR)*T[1][1] + RXADD + DXP;
				RY = (IX-DXR)*T[1][2] + RYADD + DYP;

				LL = (int)(RX);
				MM = (int)(RY);
				XD = (RX-LL)/XEP;
				YD = (RY-MM)/YEP;
				/*C  In case of projection point outside frame of projection do nothing.
				IF ((LL.LT.1).OR.(MM.LT.1).OR.(LL.GT.DPX).OR.(MM.GT.DPY)) */

				if (LL < 1 || MM < 1 || LL > DPX || MM > DPY){
					INDR = INDR + 1;
				}

				else{

					/*  Perform bilinear interpolation */
					INDA = (MM-1)*DPX + LL-1; /*-1 ???*/

					/*  LL last column of projection array?*/
					if (LL==DPX) {XD=0.0;}
					XIN1 = I[INDA] + (I[INDA+1]-I[INDA])*XD;
					/*  MM last row of projection array?*/
					if (MM==DPY)
					{YD=0.0; XIN2=0.0;}
					else
					{
						XIN2 = I[INDA+DPX] + (I[INDA+DPX+1]-I[INDA+DPX]) *XD;
					}
					YIN = XIN1 + (XIN2-XIN1) *YD;
					O[INDR] = O[INDR] + YIN;
					INDR = INDR + 1;
				}
			}
		}
	}
}


void mexFunction( int nlhs, mxArray *plhs[],
		int nrhs, const mxArray*prhs[] )

{
	float *vol;
	double *p_offs;
	const mwSize *pp_dims;
	float *p_in;
	float a1,a2;
	size_t m,n;
	size_t dima[3];
	int offa[3];
	size_t ndim;

	/* Check for proper number of arguments */

	if (nrhs != 5) {
		mexErrMsgTxt("Five input arguments required.\n Syntax: out=tom_backproj3dc(vol,in,Phi,The,Offset[])");
	} else if (nlhs > 1) {
		mexErrMsgTxt("Too many output arguments.");
	}
	ndim = mxGetNumberOfDimensions(prhs[1]);
	if (ndim != 2) {
		mexErrMsgTxt("Two dimensional image as input requiered.\n");
	}
	if ( mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) )
		mexErrMsgTxt("Input volume must be single.\n");

	/* Check the dimensions of IN. */
	/* Get the number of dimensions in the input argument. */
	m = mxGetM(IN);
	n = mxGetN(IN);

	a1 = mxGetScalar(PHI);
	a2 = mxGetScalar(THE);

	pp_dims=mxGetDimensions(prhs[0]);
	dima[0]=pp_dims[0];
	dima[1]=pp_dims[1];
	dima[2]=pp_dims[2];
	p_offs = mxGetData(OFFS);
	offa[0]=p_offs[0];
	offa[1]=p_offs[1];
	offa[2]=p_offs[2];

	/* Assign pointers to the various parameters */
	p_in = mxGetData(IN);
	vol = mxGetData(VOL);
	/* Do the actual computations in a subroutine */
	tom_backproj3d(vol,p_in,dima[0],dima[1],dima[2],a1,a2,m,n,offa[0],offa[1],offa[2]);
	return;

}


/* proj=emread('proj');proj=single(proj.Value); */
/* vol=emread('vol');vol=single(vol.Value); */

/*  i=single(zeros(128,128,128));tom_backproj3d(i,single(proj),10,12,[0 0 0]); Data.Value=double(vol)-double(i);max(max(max(Data.Value))) */



//namespace swigTom{
//
///****************************************************************************//**
// * 	\brief Wrapper for the back projection function
// * 	\param[in]
// * 	Wrapper for the back projection function. It performs a back projection of 2D projection images to return a 3D volume
// *******************************************************************************/
//template<typename T,typename TSCALE_SHIFT>
//void backProject(swigVolume<T,TSCALE_SHIFT>& src,
//		 	 	 swigVolume<T,TSCALE_SHIFT>& dst,
//		 	 	 swigVolume<T,TSCALE_SHIFT>& phi,
//		 	 	 swigVolume<T,TSCALE_SHIFT>& theta,
//		 	 	 swigVolume<T,TSCALE_SHIFT>& offset){
//
//
//	size_t vol_x, vol_y, vol_z;
//	float phi_, the_;
//	size_t img_x, img_y;
//	int off1, off2, off3;
//	int img_num;
//
//	vol_x = dst.size_x();
//	vol_y = dst.size_y();
//	vol_z = dst.size_z();
//
//	img_x = src.size_x();
//	img_y = src.size_y();
//	img_num = src.size_z();
//
//	off1 = offset.getV(0, 0, 0);
//	off2 = offset.getV(0, 0, 1);
//	off3 = offset.getV(0, 0, 2);
//
//	float vol_buffer[vol_x*vol_y*vol_z];
//	float img_buffer[img_x*img_y];
//
//	for(int i = 0; i < vol_x*vol_y*vol_z; i++){
//
//		vol_buffer[i] = 0.0;
//	}
//
//	for(int i = 0; i < img_num; i++){
//
//		for(int j = 0; j < vol_x; j++){
//			for(int k = 0; k < vol_y; k++){
//
//				img_buffer[(k*vol_x) + j] = src.getV(j, k, i);
//			}
//		}
//
//		phi_ = phi.getV(0, 0, i);
//		the_ = theta.getV(0, 0, i);
//
//		tom_backproj3d(vol_buffer, img_buffer, vol_x, vol_y, vol_z, phi_, the_, img_x, img_y, off1, off2, off3);
//	}
//
//	for(int i = 0; i < vol_x; i++){
//		for(int j = 0; j < vol_y; j++){
//			for(int k = 0; k < vol_z; k++){
//
//				dst.setV(vol_buffer[((j*vol_x) + i)+(k*vol_x*vol_y)], i, j, k);
//			}
//		}
//	}
//
//
//}
//
//}

