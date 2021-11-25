/****************************************************************************//**
 * \file swigVolumeFnc.cpp
 * \brief Source file for swigVolume functions
 * \author  Thomas Hrabe
 * \version 0.2
 * \date    1.12.2008
 **************************************************************************/


#include <swigVolume.hpp>
#include <swigVolumeFnc.hpp>
#include <tom/transf/transform.hpp>
#include <tom/volume_fcn.hpp>
#include <tom/io/io.hpp>

#include <iostream>




namespace swigTom{

/****************************************************************************//**
 * 	\brief	Reads a EM file from disk and stores it in a swigVolume.
 * 	\param[in] fileName Name of the file
 * 	The wrapped tom::read_from_em function throws exceptions that are caught outside readEM.
 *******************************************************************************/
template<typename T,typename TSCALE_SHIFT>
swigVolume<T,TSCALE_SHIFT> read(std::string fileName) {
	tom::Volume<T> *p;
	tom::io::read_from_file(p,fileName.c_str(),tom::io::UNKNOWN,NULL,NULL,NULL);
	std::auto_ptr<tom::Volume<T> > pp(p);
	return swigTom::swigVolume<T,TSCALE_SHIFT>(*(pp.get()));
}

/****************************************************************************//**
 * 	\brief	Reads a EM file from disk and stores it in a swigVolume.
 * 	\param[in] fileName Name of the file
 * 	The wrapped tom::read_from_em function throws exceptions that are caught outside readEM.
 *******************************************************************************/
template<typename T,typename TSCALE_SHIFT>
swigVolume<T,TSCALE_SHIFT> read(std::string fileName,std::size_t subregion1,std::size_t subregion2,std::size_t subregion3,
													 std::size_t subregion4,std::size_t subregion5,std::size_t subregion6,
													 std::size_t sampling1,std::size_t sampling2,std::size_t sampling3,
													 std::size_t binning1,std::size_t binning2,std::size_t binning3) {


	uint32_t* subregion = NULL;
	std::auto_ptr<uint32_t> subregionP(subregion);
	if(subregion1 >=0 && subregion2 >=0 && subregion3 >=0 && ( subregion4 >0 || subregion5 >0 || subregion6 >0) ){
		subregion = (uint32_t*)malloc(sizeof(uint32_t)*6);
		subregion[0] = (uint32_t)subregion1;
		subregion[1] = (uint32_t)subregion2;
		subregion[2] = (uint32_t)subregion3;
		subregion[3] = (uint32_t)subregion4;
		subregion[4] = (uint32_t)subregion5;
		subregion[5] = (uint32_t)subregion6;
	}

	uint32_t* sampling = NULL;
	std::auto_ptr<uint32_t> samplingP(sampling);
	if(sampling1 >0 || sampling2 >0 || sampling3 >0 ){
		sampling= (uint32_t*)malloc(sizeof(uint32_t)*3);
		sampling[0] = (uint32_t)sampling1;
		sampling[1] = (uint32_t)sampling2;
		sampling[2] = (uint32_t)sampling3;
	}

	uint32_t* binning = NULL;
	std::auto_ptr<uint32_t> binningP(binning);
	if(binning1 >0 || binning2 >0 || binning3 >0 ){
		binning = (uint32_t*)malloc(sizeof(uint32_t)*3);
		binning[0] = (uint32_t)binning1;
		binning[1] = (uint32_t)binning2;
		binning[2] = (uint32_t)binning3;
	}

	tom::Volume<T> *p;
	tom::io::read_from_file(p,fileName.c_str(),tom::io::UNKNOWN,subregion,sampling,binning);
	std::auto_ptr<tom::Volume<T> > pp(p);
	return swigTom::swigVolume<T,TSCALE_SHIFT>(*(pp.get()));
}

/****************************************************************************//**
 * 	\brief Performs correlation of v1 with conj(v2), where v2 is masked with v3
 * 	\param[in,out] v1 Search volume and result volume
 * 	\param[in] v2 The particle searched
 * 	\param[in] v3 A multiplication mask
 *
 *******************************************************************************/
template<typename T,typename TSCALE_SHIFT>
void conj_mult(swigTom::swigVolume<T,TSCALE_SHIFT>& v1,const swigTom::swigVolume<T,TSCALE_SHIFT>& v2){
	tom::element_wise_conj_multiply(v1,v2);
}

/****************************************************************************//**
 * 	\brief Wrapper function for filling the current volume with a sphere
 * 	\param[in,out] vol The volume where the sphere is stored
 * 	\param[in] radius
 * 	\param[in] sigma
 * 	\param[in] max_radius
 *
 *******************************************************************************/
template<typename T,typename TSCALE_SHIFT>
void initSphere(swigVolume<T,TSCALE_SHIFT>& vol ,float radius,float sigma, float max_radius, float centerx, float centery, float centerz){
	tom::init_sphere(vol,radius,sigma,max_radius, centerx, centery, centerz);
}

/****************************************************************************//**
 * 	\brief Wrapper for normalisation under mask
 * 	\param[in,out] vol The volume to be normalised
 * 	\param[in] mask The mask
 * 	\param[in] is_boolean_mask
 * Wrapper for normalisation under a mask. The normalisation is always performed to mean==0, std==1
 *******************************************************************************/
template<typename T,typename TSCALE_SHIFT>
void maskNorm(swigVolume<T,TSCALE_SHIFT>& vol, swigVolume<T,TSCALE_SHIFT>& mask,bool is_boolean_mask){
	T variance;
	tom::norm_mask(vol,mask,tom::norm::NORM_STDDEV_1,&variance,is_boolean_mask);
}


/****************************************************************************//**
 * 	\brief Wrapper for the volume transformation function
 * 	\param[in]
 * 	Wrapper for the volume transformation function. It performs rotation and translation defined in the transformation matrix.
 *******************************************************************************/
template<typename T,typename TSCALE_SHIFT>
void transform(swigVolume<T,TSCALE_SHIFT>& src,swigVolume<T,TSCALE_SHIFT>& dst,
			   double phi,double psi,double theta,double centerX,double centerY,double centerZ,
			   double preX,double preY,double preZ,double postX,double postY,double postZ){

	phi = tom::math::deg2rad(phi);
	psi = tom::math::deg2rad(psi);
	theta = tom::math::deg2rad(theta);

	tom::transf::rotate(src,dst,phi,psi,theta,centerX,centerY,centerZ,preX,preY,preZ,postX,postY,postZ);
}


/****************************************************************************//**
 * 	\brief Wrapper for the volume transformSpline function
 * 	\param[in]
 * 	Wrapper for the volume transformSpline function. It performs general transformation according to the given transformation matrix.
 *******************************************************************************/
template<typename T,typename TSCALE_SHIFT>
void general_transform(swigVolume<T,TSCALE_SHIFT>& src,swigVolume<T,TSCALE_SHIFT>& dst, swigVolume<T,TSCALE_SHIFT>& mtx){

	tom::transf::transformSpline(src,dst,mtx);
}


/****************************************************************************//**
 * 	\brief Wrapper for the volume transformation function
 * 	\param[in]
 * 	Wrapper for the volume transformation function. It performs a rotation around the volumes size/2+1 without any shifting
 *******************************************************************************/
template<typename T,typename TSCALE_SHIFT>
void rotate(swigVolume<T,TSCALE_SHIFT>& src,swigVolume<T,TSCALE_SHIFT>& dst,double phi,double psi,double theta){

	phi = tom::math::deg2rad(phi);
	psi = tom::math::deg2rad(psi);
	theta = tom::math::deg2rad(theta);

	int z = src.sizeZ();

	//if 2D image, rotate inplane
	if(z == 1)
		z=0;
	else
		z = z/2;

	tom::transf::rotate(src,dst,phi,psi,theta,src.sizeX()/2,src.sizeY()/2,z,0,0,0,0,0,0);

}

/****************************************************************************//**
* 	\brief Wrapper for the volume transformation function
* 	\param[in]
* 	Wrapper for the volume transformation function. It performs rotation and translation defined in the transformation matrix.
*******************************************************************************/
template<typename T,typename TSCALE_SHIFT>
void transformCubic(swigVolume<T,TSCALE_SHIFT>& src,swigVolume<T,TSCALE_SHIFT>& dst,
			   double phi,double psi,double theta,double centerX,double centerY,double centerZ,
			   double preX,double preY,double preZ,double postX,double postY,double postZ){

	phi = tom::math::deg2rad(phi);
	psi = tom::math::deg2rad(psi);
	theta = tom::math::deg2rad(theta);

	tom::transf::rotateCubic(src,dst,phi,psi,theta,centerX,centerY,centerZ,preX,preY,preZ,postX,postY,postZ);
}

/****************************************************************************//**
* 	\brief Wrapper for the volume transformation function
* 	\param[in]
* 	Wrapper for the volume transformation function. It performs a rotation around the volumes size/2+1 without any shifting
*******************************************************************************/
template<typename T,typename TSCALE_SHIFT>
void rotateCubic(swigVolume<T,TSCALE_SHIFT>& src,swigVolume<T,TSCALE_SHIFT>& dst,double phi,double psi,double theta){

	phi = tom::math::deg2rad(phi);
	psi = tom::math::deg2rad(psi);
	theta = tom::math::deg2rad(theta);

	int z = src.sizeZ();

	//if 2D image, rotate inplane
	if(z == 1)
	z=0;
	else
	z = z/2;
	tom::transf::rotateCubic(src,dst,phi,psi,theta,src.sizeX()/2,src.sizeY()/2,z,0,0,0,0,0,0);

}

/****************************************************************************//**
* 	\brief Wrapper for the volume transformation function
* 	\param[in]
* 	Wrapper for the volume transformation function. It performs rotation and translation defined in the transformation matrix.
*******************************************************************************/
template<typename T,typename TSCALE_SHIFT>
void transformSpline(swigVolume<T,TSCALE_SHIFT>& src,swigVolume<T,TSCALE_SHIFT>& dst,
			   double phi,double psi,double theta,double centerX,double centerY,double centerZ,
			   double preX,double preY,double preZ,double postX,double postY,double postZ){

	phi = tom::math::deg2rad(phi);
	psi = tom::math::deg2rad(psi);
	theta = tom::math::deg2rad(theta);

	tom::transf::rotateSpline(src,dst,phi,psi,theta,centerX,centerY,centerZ,preX,preY,preZ,postX,postY,postZ);
}

/****************************************************************************//**
* 	\brief Wrapper for the volume transformation function
* 	\param[in]
* 	Wrapper for the volume transformation function. It performs a rotation around the volumes size/2+1 without any shifting
*******************************************************************************/
template<typename T,typename TSCALE_SHIFT>
void rotateSpline(swigVolume<T,TSCALE_SHIFT>& src,swigVolume<T,TSCALE_SHIFT>& dst,double phi,double psi,double theta){

	phi = tom::math::deg2rad(phi);
	psi = tom::math::deg2rad(psi);
	theta = tom::math::deg2rad(theta);

	int z = src.sizeZ();

	//if 2D image, rotate inplane
	if(z == 1)
	z=0;
	else
	z = z/2;
	tom::transf::rotateSpline(src,dst,phi,psi,theta,src.sizeX()/2,src.sizeY()/2,z,0,0,0,0,0,0);

}


/****************************************************************************//**
 * 	\brief Wrapper for the volume shift function
 * 	\param[in]
 * 	Wrapper for the volume shiftfunction. It performs a shift of volume according to preX,preY,preZ
 *******************************************************************************/
template<typename T,typename TSCALE_SHIFT>
void shift(swigVolume<T,TSCALE_SHIFT>& src,swigVolume<T,TSCALE_SHIFT>& dst,double preX,double preY,double preZ){
	tom::transf::rotate(src,dst,0,0,0,0,0,0,preX,preY,preZ,0,0,0);
}



template<typename T,typename TSCALE_SHIFT>
void shiftFourier(swigVolume<T,TSCALE_SHIFT>& src,swigVolume<T,TSCALE_SHIFT>& dst,double shiftX,double shiftY,double shiftZ){
	tom::transf::shiftFourier(src,dst,shiftX,shiftY,shiftZ);
}

/****************************************************************************//**
 * 	\brief Wrapper for the peak function
 * 	\param[in]
 *	Wrapper for the peak function. Returns x,y,z coordinates of the highest pixel in the given volume.
 *******************************************************************************/
template<typename T>
tom::st_idx peak(const swigVolume<T,T>& volume){
	std::vector<tom::st_idx> peaks = tom::peak(volume);

	tom::st_idx peak;

	if(peaks.size() == 1)
		peak = peaks[0];
	else if(peaks.size() > 1)
		peak = peaks[0];
	else
		throw std::runtime_error("No peak found. Is the volume empty?");

	return peak;
}

/****************************************************************************//**
 * 	\brief Wrapper for the peak function
 * 	\param[in]
 *	Wrapper for the peak function. Returns x,y,z coordinates of the highest pixel in the given volume.
 *******************************************************************************/
template<typename T>
tom::st_idx peak(const swigVolume<T,T>& volume,const swigVolume<T,T>& mask){
	std::vector<tom::st_idx> peaks = tom::peak(volume, mask);

	tom::st_idx peak;

	if(peaks.size() == 1)
		peak = peaks[0];
	else if(peaks.size() > 1)
		peak = peaks[0];
	else
		throw std::runtime_error("No peak found. Is the volume empty?");
	return peak;
}

/****************************************************************************//**
 * 	\brief Wrapper for the conjugate function
 * 	\param[in]
 *	Conjugates any complex swigTom::Volume.
 *******************************************************************************/
template<typename T>
void conjugate(swigVolume<std::complex<T>,T>& volume){
	tom::conjugate(volume);
}

/****************************************************************************//**
 * 	\brief Wrapper for the paste function
 * 	\param[in]
 *	Pastes volume into destination at the given positions.
 *******************************************************************************/
template<typename T,typename TSHIFT_SCALE>
void pasteIn(const swigVolume<T,TSHIFT_SCALE>& volume,swigVolume<T,TSHIFT_SCALE>& destination,signed long x,signed long y,signed long z){
	tom::transf::paste<T,T>(volume,destination,x,y,z,NULL);
}

/****************************************************************************//**
 * 	\brief Wrapper for the pasteCenter function
 * 	\param[in]
 *	Pastes volume into the center of destination.
 *******************************************************************************/
template<typename T,typename TSHIFT_SCALE>
void pasteCenter(const swigVolume<T,TSHIFT_SCALE>& volume,swigVolume<T,TSHIFT_SCALE>& destination){
	tom::transf::paste<T,T>(volume,destination,NULL);
}

/****************************************************************************//**
 * 	\brief Wrapper for the power function
 * 	\param[in] volume
 *	Inplace power of volume to exponent
 *******************************************************************************/
template<typename T, typename T2>
void power(swigVolume<T, T2>& volume,const T2 & exponent){
	tom::element_wise_power(volume,exponent);
}

/****************************************************************************//**
 * 	\brief Wrapper for the number_set_voxels function
 * 	\param[in] volume
 *  \param[out] the number of voxels != 0
 *******************************************************************************/
template<typename T>
std::size_t numberSetVoxels(const swigVolume<T,T>& volume){
	return tom::number_set_voxels(volume);
}

/****************************************************************************//**
 * 	\brief Wrapper for the abs function
 * 	\param[in] volume
 *  \param[out] absolute values of volume
 *	Abs values will be returned as copy.
 *******************************************************************************/
template<typename T, typename TSCALE_SHIFT>
swigVolume<T,TSCALE_SHIFT> abs(const swigVolume<T,TSCALE_SHIFT>& volume){

	swigVolume<T,TSCALE_SHIFT> v(volume.getSizeX(),volume.getSizeY(),volume.getSizeZ());
	tom::element_wise_abs(v,volume);

	return v;
}


/****************************************************************************//**
 * 	\brief Returns inner sum of volume
 * 	\param[in] volume
 *  \param[out] Sum of all values
 *******************************************************************************/
template<typename T,typename TSCALE_SHIFT>
T sum(const swigVolume<T,TSCALE_SHIFT>& volume){
	return tom::sum<T, T>(volume);
}

template<typename T,typename TSCALE_SHIFT>
swigVolume<T,TSCALE_SHIFT> projectSum(const swigVolume<T,TSCALE_SHIFT>& volume, int axis){
	if(axis < 0 || axis > 2)
		throw std::runtime_error("Axis must be between 0-2!");

	std::size_t sizeX = volume.getSizeX();
	std::size_t sizeY = volume.getSizeY();
	std::size_t sizeZ = volume.getSizeZ();

	std::size_t x;
	std::size_t y;

	if(axis == 0) {
		x = sizeY;
		y = sizeZ;
	}
	else if(axis == 1) {
		x = sizeX;
		y = sizeZ;
	}
	else {
		x = sizeX;
		y = sizeY;
	}

	tom::Volume<T> res(x,y,1,NULL,NULL);

	for(int i=0; i<x; i++)
		for(int j=0; j<y; j++)
		{
			int sum = 0;

			if(axis == 0) {
				for(int ii=0; ii<sizeX; ii++)
					sum += volume.get(ii,i,j);
			}
			else if(axis == 1) {
				for(int jj=0; jj<sizeY; jj++)
					sum += volume.get(i,jj,j);
			}
			else {
				for(int kk=0; kk<sizeZ; kk++)
					sum += volume.get(i,j,kk);
			}

			res.get(i,j,0) = sum;
		}

	return  swigTom::swigVolume<T,TSCALE_SHIFT>(res);
}

/****************************************************************************//**
 * 	\brief Returns mean of volume
 * 	\param[in] volume
 *  \param[out] Mean of all values
 *******************************************************************************/
template<typename T,typename TSCALE_SHIFT>
double mean(const swigVolume<T,TSCALE_SHIFT>& volume){
	return tom::mean(volume);
}

/****************************************************************************//**
 * 	\brief Returns variance of volume
 * 	\param[in] volume
 *  \param[out] Variance of all values
 *******************************************************************************/
template<typename T,typename TSCALE_SHIFT>
double variance(const swigVolume<T,TSCALE_SHIFT>& volume,bool use_sample_standard_deviation){
	return tom::variance(volume,use_sample_standard_deviation);
}

/****************************************************************************//**
 * 	\brief Returns min of volume
 * 	\param[in] volume
 *  \param[out] Min of all values
 *******************************************************************************/
template<typename T,typename TSCALE_SHIFT>
T min(const swigVolume<T,TSCALE_SHIFT>& volume){
	return tom::min(volume);
}

/****************************************************************************//**
 * 	\brief Returns max of volume
 * 	\param[in] volume
 *  \param[out] Max of all values
 *******************************************************************************/
template<typename T,typename TSCALE_SHIFT>
T max(const swigVolume<T,TSCALE_SHIFT>& volume){
	return tom::max(volume);
}

/****************************************************************************//**
 * 	\brief Replaces values lower or higher than bounds with other values
 * 	\param[in] volume
 *******************************************************************************/
template<typename T,typename TSCALE_SHIFT>
void limit(swigVolume<T,TSCALE_SHIFT>& volume, T lowerBound, T lowerReplacement, T upperBound, T upperReplacement,bool doLower, bool doUpper){
	tom::limit(volume,lowerBound,lowerReplacement,upperBound,upperReplacement,doLower,doUpper);
}

/****************************************************************************//**
 * 	\brief Expands volume from reduced complex to hermitian symetry
 * 	\param[in] volume
 *  \param[out] Full, complex volume
 *******************************************************************************/
template<typename T,typename TSCALE_SHIFT>
swigTom::swigVolume<T,TSCALE_SHIFT> reducedToFull(swigVolume<T,TSCALE_SHIFT>& volume){
	std::size_t x,y,z;
	if(((std::size_t)volume.getFtSizeX()) == 0 && ((std::size_t)volume.getFtSizeY()) == 0 && ((std::size_t)volume.getFtSizeZ()) == 0){
		x = volume.sizeX();
		y = volume.sizeY();
		z = (volume.sizeZ()-1)*2;
		if((z-y == 1 || z-y == -1) && y == x)
			z = y;
	}
	else{
		x = (std::size_t) volume.getFtSizeX();
		y = (std::size_t) volume.getFtSizeY();
		z = (std::size_t) volume.getFtSizeZ();
	}

	tom::Volume<T> resultVolume(x,y,z,NULL,NULL);
	tom::Volume<T> windowVolume(resultVolume,(void*)(&resultVolume.get()),volume.getSizeX(),volume.getSizeY(),volume.getSizeZ(),volume.getStrideX(),volume.getStrideY(),volume.getStrideZ());

	windowVolume.setValues(volume);

	tom::hermitian_symmetry_to_full(resultVolume);

	return  swigTom::swigVolume<T,TSCALE_SHIFT>(resultVolume);

}

/****************************************************************************//**
* 	\brief Expands volume from full complex to half size
* 	\param[in] Full size volume
*  \param[out] Reduced, complex volume
*******************************************************************************/
template<typename T,typename TSCALE_SHIFT>
swigTom::swigVolume<T,TSCALE_SHIFT> fullToReduced(swigVolume<T,TSCALE_SHIFT>& volume){
    
    std::size_t sizeX,sizeY,sizeZ;
    sizeX = volume.getSizeX();
	sizeY = volume.getSizeY();
	// sizeZ = (volume.getSizeZ()+1)/2;
	sizeZ = volume.getSizeZ()/2 + 1;

    swigVolume<T,TSCALE_SHIFT> subvolume(tom::getSubregion(volume,0,0,0,sizeX,sizeY,sizeZ));
    subvolume.setFtSizeX(sizeX);
    subvolume.setFtSizeY(sizeY);
    subvolume.setFtSizeZ(volume.getSizeZ());
    
    return subvolume;
}
    
/****************************************************************************//**
 *	\brief Assigns gaussian noise to each volume voxel
 *	\param[in] volume Volume
 *  \param[in] mean Noise mean
 *  \param[in] sigma Noise deviation
 *******************************************************************************/
template<typename T,typename TSCALE_SHIFT>
void gaussianNoise(swigVolume<T,TSCALE_SHIFT>& volume,double mean, double sigma){
	tom::gaussian_noise(volume,mean,sigma);
}

/****************************************************************************//**
 * 	\brief Divides volume by div
 * 	\param[in] Complex source volume
 *  \param[in] Real volume
 *  \param[return]
 *******************************************************************************/
template<typename T>
swigVolume<std::complex<T>,T > complexDiv(swigVolume<std::complex<T>,T> &volume, const swigVolume<T,T> &div){
	swigVolume<std::complex<T>,T> vol(volume);
	tom::element_wise_div<std::complex<T>,T>(vol,div,std::complex<T>(0,0));
	return vol;
}

/****************************************************************************//**
 * 	\brief Converts a  3d volume to a 1d vector
 *******************************************************************************/
template<typename T>
swigVolume<T,T> vectorize(swigVolume<T,T> &source){
	return swigVolume<T,T>(tom::vectorize(source));
}


/****************************************************************************//**
 * 	\brief Returns a subvolume from source
 *******************************************************************************/
template<typename T,typename TSCALE_SHIFT>
swigVolume<T,TSCALE_SHIFT> getSubregion(swigVolume<T,TSCALE_SHIFT> &source,std::size_t startX,std::size_t startY,std::size_t startZ,std::size_t endX,std::size_t endY,std::size_t endZ){
	return swigVolume<T,TSCALE_SHIFT>(tom::getSubregion(source,startX,startY,startZ,endX,endY,endZ));
}

/****************************************************************************//**
 * 	\brief Puts a small volume (source) into a bigger one (destination)
 *******************************************************************************/
template<typename T,typename TSCALE_SHIFT>
void putSubregion(swigVolume<T,TSCALE_SHIFT> &source, swigVolume<T,TSCALE_SHIFT> &destination,std::size_t positionX,std::size_t positionY,std::size_t positionZ){
	tom::putSubregion(source,destination,positionX,positionY,positionZ);
}

/**
 * @brief Write a subregion into a existing file
 */
template<typename T>
void writeSubregion(swigVolume<T,T> &source, std::string filename, std::size_t positionX, std::size_t positionY, std::size_t positionZ)
{
	uint32_t first_voxel[3]={positionX, positionY, positionZ};
	tom::io::write_to_file_paste(source, filename, first_voxel);
}

/**
 * @brief Used for module localization. Update the result volume with higher value, and corresponding orientation information.
 *
 * @param[out] resV to be updated volume
 * @param[in] newV newly calculated volume
 * @param[out] orientV to be updated orientation volume
 * @param[in] orientIdx newly calculated orientation index
 *
 * @author Yuxiang Chen
 */
template<typename T,typename TSHIFT_SCALE>
void updateResFromIdx(swigVolume<T,TSHIFT_SCALE>& resV,const swigVolume<T,TSHIFT_SCALE>& newV, swigVolume<T,TSHIFT_SCALE>& orientV, const std::size_t orientIdx)
{
	std::size_t sizeX = newV.getSizeX();
	std::size_t sizeY = newV.getSizeY();
	std::size_t sizeZ = newV.getSizeZ();
	if(sizeX > resV.getSizeX() || sizeY > resV.getSizeY() || sizeZ > resV.getSizeZ()
		|| resV.getSizeX() != orientV.getSizeX() || resV.getSizeY() != orientV.getSizeY() || resV.getSizeZ() != orientV.getSizeZ())
		throw std::runtime_error("Volumes sizes not consistent!");

	for(int i=0; i<sizeX; i++)
		for(int j=0; j<sizeY; j++)
			for(int k=0; k<sizeZ; k++)
			{
				if(resV.get(i,j,k) < newV.get(i,j,k))
				{
					resV.get(i,j,k) = newV.get(i,j,k);
					orientV.get(i,j,k) = orientIdx;
				}
			}
}

/**
 * @brief Used for module localization. Update the result volume with higher value, and corresponding orientation information.
 *
 * @param[out] resV to be updated volume
 * @param[in] newV newly calculated volume
 * @param[out] orientV to be updated orientation volume
 * @param[in] neworientV newly calculated orientation volume
 *
 * @author Yuxiang Chen
 */
template<typename T,typename TSHIFT_SCALE>
void updateResFromVol(swigVolume<T,TSHIFT_SCALE>& resV,const swigVolume<T,TSHIFT_SCALE>& newV, swigVolume<T,TSHIFT_SCALE>& orientV, const swigVolume<T,TSHIFT_SCALE>& neworientV)
{
	std::size_t sizeX = resV.getSizeX();
	std::size_t sizeY = resV.getSizeY();
	std::size_t sizeZ = resV.getSizeZ();
	if(sizeX != newV.getSizeX() || sizeY != newV.getSizeY() || sizeZ != newV.getSizeZ()
		|| sizeX != orientV.getSizeX() || sizeY != orientV.getSizeY() || sizeZ != orientV.getSizeZ()
		|| sizeX != neworientV.getSizeX() || sizeY != neworientV.getSizeY() || sizeZ != neworientV.getSizeZ())
		throw std::runtime_error("Volumes have different sizes!");

	for(int i=0; i<sizeX; i++)
		for(int j=0; j<sizeY; j++)
			for(int k=0; k<sizeZ; k++)
			{
				if(resV.get(i,j,k) < newV.get(i,j,k))
				{
					resV.get(i,j,k) = newV.get(i,j,k);
					orientV.get(i,j,k) = neworientV.get(i,j,k);
				}
			}
}

template<typename T,typename TSHIFT_SCALE>
void mirrorVolume(swigVolume<T,TSHIFT_SCALE>& src, swigVolume<T,TSHIFT_SCALE>& des)
{
	std::size_t sizeX = src.getSizeX();
	std::size_t sizeY = src.getSizeY();
	std::size_t sizeZ = src.getSizeZ();
	if(sizeX != des.getSizeX() || sizeY != des.getSizeY() || sizeZ != des.getSizeZ())
		throw std::runtime_error("Volumes have different sizes!");

	// swap according to the center of x
	int l = sizeX-1;
	for(int i=0; i<sizeX; i++)
		for(int j=0; j<sizeY; j++)
			for(int k=0; k<sizeZ; k++)
			{
				des.get(l-i,j,k) = src.get(i,j,k);
			}
}

template<typename T,typename TSHIFT_SCALE>
void rescale(const swigVolume<T,TSHIFT_SCALE>& source,swigVolume<T,TSHIFT_SCALE>& destination){
	tom::transf::rescale(source,destination);
}

template<typename T,typename TSHIFT_SCALE>
void rescaleCubic(const swigVolume<T,TSHIFT_SCALE>& source,swigVolume<T,TSHIFT_SCALE>& destination){
	tom::transf::rescale(source,destination,tom::transf::InterpolTriCubic<T>(0));
}

template<typename T,typename TSHIFT_SCALE>
void rescaleSpline(const swigVolume<T,TSHIFT_SCALE>& source,swigVolume<T,TSHIFT_SCALE>& destination){
	tom::transf::rescale(source,destination,tom::transf::InterpolCubicSpline<T>(0));
}

template<typename T,typename TSHIFT_SCALE>
T interpolate(const swigVolume<T,TSHIFT_SCALE>& source,const float x,const float y,const float z){
	const T defaultValue = 0;
	tom::transf::InterpolTriLinear<T> interpolationObject = tom::transf::InterpolTriLinear<T>(defaultValue);

	interpolationObject.setvolume(source);

	return interpolationObject.interpolate(x,y,z);
}

template<typename T,typename TSHIFT_SCALE>
T interpolateCubic(const swigVolume<T,TSHIFT_SCALE>& source,const float x,const float y,const float z){
	const T defaultValue = 0;
	tom::transf::InterpolTriCubic<T> interpolationObject = tom::transf::InterpolTriCubic<T>(defaultValue);

	interpolationObject.setvolume(source);

	return interpolationObject.interpolate(x,y,z);
}

template<typename T,typename TSHIFT_SCALE>
T interpolateSpline(const swigVolume<T,TSHIFT_SCALE>& source,const float x,const float y,const float z){
	const T defaultValue = 0;
	tom::transf::InterpolCubicSpline<T> interpolationObject = tom::transf::InterpolCubicSpline<T>(defaultValue);

	interpolationObject.setvolume(source);

	return interpolationObject.interpolate(x,y,z);
}

template<typename T,typename TSHIFT_SCALE>
T interpolateFourierSpline(const swigVolume<T,TSHIFT_SCALE>& source,const float x,const float y,const float z){
	const T defaultValue = 0;
	tom::transf::InterpolFourierSpline<T> interpolationObject = tom::transf::InterpolFourierSpline<T>(defaultValue);

	interpolationObject.setvolume(source);

	return interpolationObject.interpolate(x,y,z);
}

/****************************************************************************//**
 * 	\brief Wrapper for the back projection function
 * 	\param[in]
 * 	Wrapper for the back projection function. It performs a back projection of 2D projection images to return a 3D volume
 *******************************************************************************/


template<typename T,typename TSCALE_SHIFT>
void backProject(swigVolume<T,TSCALE_SHIFT>& src,
		 	 	 swigVolume<T,TSCALE_SHIFT>& dst,
		 	 	 swigVolume<T,TSCALE_SHIFT>& phi,
		 	 	 swigVolume<T,TSCALE_SHIFT>& theta,
		 	 	 swigVolume<T,TSCALE_SHIFT>& offsetPaticle ,
		 	 	 swigVolume<T,TSCALE_SHIFT>& offsetProjections){


	int vol_x, vol_y, vol_z;
	float phi_, the_;
	int img_x, img_y;
	int img_num;
	float off1, off2, off3;

	vol_x = dst.sizeX();
	vol_y = dst.sizeY();
	vol_z = dst.sizeZ();

	img_x = src.sizeX();
	img_y = src.sizeY();
	img_num = src.sizeZ();

	off1 = offsetPaticle.getV(0, 0, 0);
	off2 = offsetPaticle.getV(1, 0, 0);
	off3 = offsetPaticle.getV(2, 0, 0);


	dst.setAll(0.0);

	for(int i = 0; i < img_num; i++){

		phi_ = phi.getV(0, 0, i);
		the_ = theta.getV(0, 0, i);

		swigVolume<T,TSCALE_SHIFT> slice = getSubregion(src, 0, 0, i, img_x, img_y, 1);
		//std::cout << "HALLO" << std::endl;
		tom::back_project(slice, dst,
						  vol_x, vol_y, vol_z,
						  phi_, the_, 0,
						  img_x, img_y,
						  off1, off2, off3,offsetProjections.getV(0,0,i),offsetProjections.getV(0,1,i));
	}

}

template<typename T,typename TSCALE_SHIFT>
void backProject(swigVolume<T,TSCALE_SHIFT>& src,
				 swigVolume<T,TSCALE_SHIFT>& dst,
				 swigVolume<T,TSCALE_SHIFT>& phi,
				 swigVolume<T,TSCALE_SHIFT>& theta,
				 swigVolume<T,TSCALE_SHIFT>& psi,
				 swigVolume<T,TSCALE_SHIFT>& offsetPaticle ,
				 swigVolume<T,TSCALE_SHIFT>& offsetProjections){
		
		
		int vol_x, vol_y, vol_z;
		float phi_, the_, psi_;
		int img_x, img_y;
		int img_num;
		
		vol_x = dst.sizeX();
		vol_y = dst.sizeY();
		vol_z = dst.sizeZ();
		
		img_x = src.sizeX();
		img_y = src.sizeY();
		img_num = src.sizeZ();
		
		dst.setAll(0.0);
		
		for(int i = 0; i < img_num; i++){
			
			phi_ = phi.getV(0, 0, i);
			the_ = theta.getV(0, 0, i);
			psi_ = psi.getV(0, 0, i);
			swigVolume<T,TSCALE_SHIFT> slice = getSubregion(src, 0, 0, i, img_x, img_y, 1);

			tom::back_project(slice, dst,
							  vol_x, vol_y, vol_z,
							  phi_, the_, psi_,
							  img_x, img_y,
							  offsetPaticle.getV(0,i,0), offsetPaticle.getV(1,i,0),offsetPaticle.getV(2,i,0),
							  offsetProjections.getV(0,0,i),offsetProjections.getV(0,1,i));
		}
		
}

template<typename T, typename TSCALE_SHIFT>
swigVolume<std::complex<T>,TSCALE_SHIFT> complexRealMult(const swigVolume<std::complex<T>,TSCALE_SHIFT> &vol,
														 const swigVolume<T,TSCALE_SHIFT> &otherVol){
	swigVolume<std::complex<T>,TSCALE_SHIFT> newVol(vol);
	tom::element_wise_multiply(newVol, otherVol);
	return newVol;
}

template<typename T, typename TSCALE_SHIFT>
swigVolume<T,TSCALE_SHIFT> real(const swigVolume<std::complex<T>,TSCALE_SHIFT> &vol){
        
    swigVolume<T,TSCALE_SHIFT> realVol = swigVolume<T,TSCALE_SHIFT>(vol.sizeX(),vol.sizeY(),vol.sizeZ());
        
    std::complex<T> val=0;
        
    for(std::size_t x = 0; x<vol.sizeX();x++){
        for(std::size_t y = 0; y<vol.sizeY();y++){
            for(std::size_t z = 0; z<vol.sizeZ();z++){
                val = vol.get(x,y,z);
                realVol(val.real(),x,y,z);
            }
        }
    }
        
    return realVol;
        
}

template<typename T, typename TSCALE_SHIFT>
swigVolume<T,TSCALE_SHIFT> imag(const swigVolume<std::complex<T>,TSCALE_SHIFT> &vol){
        
    swigVolume<T,TSCALE_SHIFT> imagVol = swigVolume<T,TSCALE_SHIFT>(vol.sizeX(),vol.sizeY(),vol.sizeZ());
        
    std::complex<T> val=0;
        
    for(std::size_t x = 0; x<vol.sizeX();x++){
        for(std::size_t y = 0; y<vol.sizeY();y++){
            for(std::size_t z = 0; z<vol.sizeZ();z++){
                val = vol.get(x,y,z);
                imagVol(val.imag(),x,y,z);
            }
        }
    }
        
    return imagVol;
        
}
    
template<typename T, typename TSCALE_SHIFT>
swigVolume<std::complex<T>,TSCALE_SHIFT> mergeRealImag(const swigVolume<T,TSCALE_SHIFT> &real,const swigVolume<T,TSCALE_SHIFT> &imag){
    
    if(real.getSizeX() != imag.getSizeX() || real.getSizeY() != imag.getSizeY() || real.getSizeZ() != imag.getSizeZ()){
        throw std::runtime_error("Volumes in mergeRealImag must have same size!");
    }
    
    swigVolume<std::complex<T>,TSCALE_SHIFT> complexVol = swigVolume<std::complex<T>,TSCALE_SHIFT>(real.sizeX(),real.sizeY(),real.sizeZ());
    
    for(std::size_t x = 0; x<real.sizeX();x++){
        for(std::size_t y = 0; y<real.sizeY();y++){
            for(std::size_t z = 0; z<real.sizeZ();z++){
                std::complex<T> val(real.get(x,y,z), imag.get(x,y,z));
                complexVol(val,x,y,z);
            }
        }
    }
    
    return complexVol;
}


template<typename T, typename TSCALE_SHIFT>
swigVolume<T,TSCALE_SHIFT> volTOsf(const swigVolume<T,TSCALE_SHIFT> &vol, const double r, const int b, const float m_x, const float m_y, const float m_z){
    
    swigVolume<T,TSCALE_SHIFT> res = swigVolume<T,TSCALE_SHIFT>(4*b*b, 1, 1);
    
    double the = 0.0;
    double phi = 0.0;
    float x = 0.0;
    float y = 0.0;
    float z = 0.0;
    T c;

    for(int jj=0; jj < 2*b; jj++)
    {
		for(int kk=0; kk < 2*b; kk++)
		{
			the = M_PI * (2*jj+1) / (4*b);
			phi = M_PI * kk / b;
			x = r*cos(phi)*sin(the);
			y = r*sin(phi)*sin(the);
			z = r*cos(the);
			c = interpolateSpline(vol, x+m_x, y+m_y, z+m_z);
			res(c, jj*2*b+kk, 0,0);
		}
    }
    
    return res;
}

template<typename T, typename TSCALE_SHIFT>
swigVolume<std::complex<T>,TSCALE_SHIFT> fvolTOsf(const swigVolume<std::complex<T>,TSCALE_SHIFT> &vol, const double r, const int b){
	// assume the zero frequency of the input volume has already been shifted in the center
	// construct the result complex volume
    swigVolume<std::complex<T>,TSCALE_SHIFT> res = swigVolume<std::complex<T>,TSCALE_SHIFT>(4*b*b, 1, 1);
    
    double the = 0.0;
    double phi = 0.0;
    float x = 0.0;
    float y = 0.0;
    float z = 0.0;
    int m_x = vol.getSizeX()/2;
    int m_y = vol.getSizeY()/2;
    int m_z = vol.getSizeZ()/2;

    swigVolume<T,TSCALE_SHIFT> fvol = real(vol);
    swigVolume<T,TSCALE_SHIFT> ivol = imag(vol);
    T re;
    T im;

    for(int jj=0; jj < b; jj++) // only loop half sphere
    {
		for(int kk=0; kk < 2*b; kk++)
		{
			the = M_PI * (2*jj+1) / (4*b);
			phi = M_PI * kk / b;
			x = r*cos(phi)*sin(the);
			y = r*sin(phi)*sin(the);
			z = r*cos(the);
			re = interpolateSpline(fvol, x+m_x, y+m_y, z+m_z);
			im = interpolateSpline(ivol, x+m_x, y+m_y, z+m_z);

			res(std::complex<T>(re, im), jj*2*b+kk, 0,0);
			res(std::complex<T>(re,-im), (2*b-jj-1)*2*b+(kk<b?kk+b:kk-b), 0,0); // set the other half
		}
    }
        
    return res; 
}


// template<typename T, typename TSCALE_SHIFT>
// swigVolume<T,TSCALE_SHIFT> fvolGetSF(const swigVolume<std::complex<T>,TSCALE_SHIFT> &vol, const int r, const int b){
// 	// assume the zero frequency of the input volume has already been shifted in the center
	
//     swigVolume<T,TSCALE_SHIFT> res = swigVolume<T,TSCALE_SHIFT>(4*b*b, 1, 1);
    
//     double the = 0.0;
//     double phi = 0.0;
//     double x = 0.0;
//     double y = 0.0;
//     double z = 0.0;
//     int m_x = vol.getSizeX()/2;
//     int m_y = vol.getSizeY()/2;
//     int m_z = vol.getSizeZ()/2;

//     swigVolume<std::complex<T>,TSCALE_SHIFT> tmp = abs(vol);
//     swigVolume<T,TSCALE_SHIFT> amp_vol = real(tmp);

//     for(int jj=0; jj < b; jj++) // only loop half sphere
//     {
// 		for(int kk=0; kk < 2*b; kk++)
// 		{
// 			the = M_PI * (2*jj+1) / (4*b);
// 			phi = M_PI * kk / b;
// 			x = r*cos(phi)*sin(the);
// 			y = r*sin(phi)*sin(the);
// 			z = r*cos(the);
// 			T amp = interpolateSpline(amp_vol, x+m_x, y+m_y, z+m_z);

// 			res(amp, jj*2*b+kk, 0,0);
// 			res(amp, (2*b-jj-1)*2*b+(kk<b?kk+b:kk-b), 0,0); // set the other half
// 		}
//     }
        
//     return res; 
// }


template<typename T, typename TSCALE_SHIFT>
swigVolume<std::complex<T>,TSCALE_SHIFT> sfFourierShift(const swigVolume<std::complex<T>,TSCALE_SHIFT> &vol, const double r, const int b, int sizex, int sizey, int sizez, double shiftX,double shiftY,double shiftZ){
	// assume the zero frequency of the input volume has already been shifted in the center
	// validate the input vol
	if(vol.getSizeX()!=4*b*b || vol.getSizeY()!=1 || vol.getSizeZ()!=1)
		throw std::runtime_error("Volumes size not right!"); 

	// construct the result complex volume
    swigVolume<std::complex<T>,TSCALE_SHIFT> res = swigVolume<std::complex<T>,TSCALE_SHIFT>(4*b*b, 1, 1);
    
    double the = 0.0;
    double phi = 0.0;
    float x = 0.0;
    float y = 0.0;
    float z = 0.0;
    // int m_x = sizex/2; // the center of fourier volume
    // int m_y = sizey/2;
    // int m_z = sizez/2;

    std::complex<T> c;

    for(int jj=0; jj < b; jj++) // only loop half sphere
    {
		for(int kk=0; kk < 2*b; kk++)
		{
			the = M_PI * (2*jj+1) / (4*b);
			phi = M_PI * kk / b;
			x = r*cos(phi)*sin(the);
			y = r*sin(phi)*sin(the);
			z = r*cos(the);

			std::complex<T> shift = std::complex<T>(cos(2*M_PI/sizex*x*shiftX), -sin(2*M_PI/sizex*x*shiftX)) *
							  		std::complex<T>(cos(2*M_PI/sizey*y*shiftY), -sin(2*M_PI/sizey*y*shiftY)) *
							  		std::complex<T>(cos(2*M_PI/sizez*z*shiftZ), -sin(2*M_PI/sizez*z*shiftZ));

			c = vol.get(jj*2*b+kk, 0, 0) * shift;

			res(c, jj*2*b+kk, 0,0);
			res(std::conj(c), (2*b-jj-1)*2*b+(kk<b?kk+b:kk-b), 0,0); // set the other half
		}
    }
        
    return res; 
}


/*******************************************************************************
* 	\brief Wrapper for the volume rotation function in fourier space with spline interpolation
* 	\param[in]
* 	Wrapper for the volume rotate function. It performs a rotation around the volumes size/2+1 without any shifting
*******************************************************************************/
template<typename T,typename TSCALE_SHIFT>
swigVolume<T,TSCALE_SHIFT> rotateSplineInFourier(swigVolume<T,TSCALE_SHIFT>& src,double phi,double psi,double theta){
    
    phi = tom::math::deg2rad(phi);
    psi = tom::math::deg2rad(psi);
    theta = tom::math::deg2rad(theta);

    swigTom::swigVolume<T,TSCALE_SHIFT> full = reducedToFull(src);
    
    tom::fftshift(full,false);
    swigVolume<T,TSCALE_SHIFT> dst(full.sizeX(),full.sizeY(),full.sizeZ());

	std::size_t z = full.sizeZ();
    
    //if 2D image, rotate inplane
    if(z == 1)
        z=0;
    else
        z = z/2;

    tom::transf::rotateFourierSpline(full,dst,phi,psi,theta,(double)full.sizeX()/2,(double)full.sizeY()/2,(double)z,(T)0);
    tom::fftshift(dst,true);

    return fullToReduced(dst);

}    
    
/*******************************************************************************
* 	\brief Wrapper for the volume rotation function in fourier space with spline interpolation and phase shift for shift
* 	\param[in]
* 	Wrapper for the volume rotate function. It performs a rotation around the volumes size/2+1 without any shifting
*******************************************************************************/
template<typename T,typename TSCALE_SHIFT>
swigVolume<T,TSCALE_SHIFT> transformFourierSpline(swigVolume<T,TSCALE_SHIFT>& src,double phi,double psi,double theta,double shiftX,double shiftY,double shiftZ){

    phi = tom::math::deg2rad(phi);
    psi = tom::math::deg2rad(psi);
    theta = tom::math::deg2rad(theta);

    swigTom::swigVolume<T,TSCALE_SHIFT> full = reducedToFull(src);

    tom::fftshift(full,false);
    swigVolume<T,TSCALE_SHIFT> dst(full.sizeX(),full.sizeY(),full.sizeZ());

	std::size_t z = full.sizeZ();

    //if 2D image, rotate inplane
    if(z == 1)
        z=0;
    else
        z = z/2;

    tom::transf::transformFourierSpline(full,dst,phi,psi,theta,(double)full.sizeX()/2,(double)full.sizeY()/2,(double)z, shiftX, shiftY, shiftZ,(T)0);
    tom::fftshift(dst,true);

    return fullToReduced(dst);

}

}
