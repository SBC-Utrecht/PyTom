/****************************************************************************//**
 * \file swigVolume.hpp
 * \brief The header file for the class swigTom::swigVolume.
 * \author  Thomas Hrabe
 * \version 0.2
 * \date    1.12.2008
 *******************************************************************************/

#ifndef SWIGVOLUME_HPP_
#define SWIGVOLUME_HPP_


#include <tom/volume.hpp>
#include <fftw3.h>
#include <tom/volume_fcn.hpp>
#include <tom/io/io.h>

namespace swigTom{
/**********************************************************************//**
 *	\brief Volume class used for wraping tom::Volume
 *
 *
 *	swigVolume wraps tom::Volume<T> for interfacing with SWIG.
 * 	Methods defined here will be visible in the interfaced language such as Python.
 *	There will not be much of a documentation for class methods because they basically wrap baseclass methods.
 *	Template parameters have become complicated because of the predefined template functions / members of the parent classes.
 */
template<typename T, typename TSCALE_SHIFT>
class swigVolume : public tom::Volume<T>{

private:
	float ftSizeX;
	float ftSizeY;
	float ftSizeZ;

public:

	//swigVolume(std::size_t sizex,std::size_t sizey,std::size_t sizez): tom::Volume<T>(sizex,sizey,sizez,&fftw_malloc,&fftw_free){ this->ftSizeX = 0;this->ftSizeY = 0;this->ftSizeZ = 0; };
	swigVolume(std::size_t sizex,std::size_t sizey,std::size_t sizez): tom::Volume<T>(sizex,sizey,sizez,NULL,NULL){ this->ftSizeX = 0;this->ftSizeY = 0;this->ftSizeZ = 0; };
	swigVolume(const tom::Volume<T>& v) : tom::Volume<T>(v){ this->ftSizeX = 0;this->ftSizeY = 0;this->ftSizeZ = 0; };
	swigVolume(const swigVolume<T,TSCALE_SHIFT>& v) : tom::Volume<T>(v){ this->ftSizeX = v.getFtSizeX();this->ftSizeY = v.getFtSizeY();this->ftSizeZ = v.getFtSizeZ(); };
	
	swigVolume(T *data, std::size_t sizex, std::size_t sizey, std::size_t sizez, std::size_t stridex, std::size_t stridey, std::size_t stridez) : tom::Volume<T>(data, sizex, sizey, sizez, stridex, stridey, stridez, false, NULL)
	{
		this->ftSizeX = 0;this->ftSizeY = 0;this->ftSizeZ = 0;
	};

	~swigVolume(){};


	/**
	 *	\brief Wrapper function
	 */
	std::size_t sizeX() const{
		return this->getSizeX();
	};
	/**
	 *	\brief Wrapper function
	 */
	std::size_t sizeY() const{
		return this->getSizeY();
	};
	/**
	 *	\brief Wrapper function
	 */
	std::size_t sizeZ() const{
		return this->getSizeZ();
	};

	float getFtSizeX() const{
		return this->ftSizeX;
	};
	/**
	 *	\brief Wrapper function
	 */
	float getFtSizeY() const{
		return this->ftSizeY;
	};
	/**
	 *	\brief Wrapper function
	 */
	float getFtSizeZ() const{
		return this->ftSizeZ;
	};

	void setFtSizeX(float sizeX){
		this->ftSizeX=sizeX;
	};
	/**
	*	\brief Wrapper function
	*/
	void setFtSizeY(float sizeY){
		this->ftSizeY=sizeY;
	};
	/**
	*	\brief Wrapper function
	*/
	void setFtSizeZ(float sizeZ){
		this->ftSizeZ=sizeZ;
	};
	/**
	 *	\brief Wrapper function
	 */

	std::size_t strideX() const{
		return this->getStrideX();
	};
	/**
	 *	\brief Wrapper function
	 */
	std::size_t strideY() const{
		return this->getStrideY();
	};
	/**
	 *	\brief Wrapper function
	 */
	std::size_t strideZ() const{
		return this->getStrideZ();
	};
	/**
	 *	\brief Wrapper function
	 */
	void write(std::string fileName){
		tom::io::write_to_em(*this,fileName,NULL);
	};

	void write(std::string fileName,std::string fileType){

		if(fileType == "mrc")
			tom::io::write_to_mrc(*this,fileName,NULL);
		else if(fileType == "ccp4")
			tom::io::write_to_ccp4(*this,fileName,NULL);
		else if(fileType == "em")
			tom::io::write_to_em(*this,fileName,NULL);
		else{
			std::cout << std::endl << "Filetype unknown, saving data as EM to " << fileName << std::endl;
			tom::io::write_to_em(*this,fileName,NULL);
		}
	};


	/**
	 *	\brief Wrapper function
	 */
	void info (const std::string &name) const{
		this->printInfo(name);
	}
	/**
	 *	\brief Wrapper function
	 */
	std::size_t numelem() const{
		return this->numel();
	}
	/**
	 *	\brief Wrapper function
	 */
	bool equalsTo(const swigTom::swigVolume<T,TSCALE_SHIFT> &v) const{
		return (*this) == v;
	}
	/**
	 *	\brief Wrapper function
	 */
	T getV(std::size_t x,std::size_t y, std::size_t z){
		return this->get(x,y,z);
	}
	/**
	 *	\brief Wrapper function
	 */
	void setV(T val,std::size_t x,std::size_t y, std::size_t z){
		this->get(x,y,z) = val;
	}
	/**
	 *	\brief Wrapper function
	 */
	void setAll(T val){
		this->setValues(val);
	}
	/**
	 *	\brief Wrapper function
	 */
	void copyVolume(const swigTom::swigVolume<T,TSCALE_SHIFT>& v2){
		this->setValues(v2);
		this->ftSizeX = v2.getFtSizeX();
		this->ftSizeY = v2.getFtSizeY();
		this->ftSizeZ = v2.getFtSizeZ();
	}
	/**
	 *	\brief Wrapper function
	 */
	void shiftscale(const TSCALE_SHIFT & shiftV, const TSCALE_SHIFT & scaleV){
		this->shift_scale(shiftV,scaleV);
	}

	/**
	 *	\brief Wrapper function
	 */
	std::size_t dims(){
		if(this->getSizeZ() >0)
			return 3;
		else if(this->getSizeY() >0)
			return 2;
		else
			return 1;
	}

	swigVolume<T,TSCALE_SHIFT> operator+(const swigVolume<T,TSCALE_SHIFT> &otherVol) const;
	swigVolume<T,TSCALE_SHIFT> operator+(const TSCALE_SHIFT &value) const;
	swigVolume<T,TSCALE_SHIFT> operator*(const swigVolume<T,TSCALE_SHIFT> &otherVol) const;
	swigVolume<T,TSCALE_SHIFT> operator*(const TSCALE_SHIFT &value) const;
	swigVolume<T,TSCALE_SHIFT> operator-(const swigVolume<T,TSCALE_SHIFT> &otherVol) const;
	swigVolume<T,TSCALE_SHIFT> operator-(const TSCALE_SHIFT &value) const;
	swigVolume<T,TSCALE_SHIFT> operator/(const swigVolume<T,TSCALE_SHIFT> &otherVol) const;
	swigVolume<T,TSCALE_SHIFT> operator/(const TSCALE_SHIFT &value) const;
	const T operator()(const std::size_t &x,const std::size_t &y,const std::size_t &z);
	void operator()(const T &value,const std::size_t &x,const std::size_t &y,const std::size_t &z);
};


// The following is the class for EM header
/*
class EMHeader {

	tom_io_em_header raw_data;

public:

	EMHeader() {
		this->raw_data.machine = 6; // set the machine code to PC
	};
	
	~EMHeader() {
		delete raw_data;
	};
	
	void set_voltage(long voltage)
	{
		this->raw_data.emdata[0] = voltage;
	}
	
};
*/


}


#endif
