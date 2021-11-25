/****************************************************************************//**
 * \file swigVolume.hpp
 * \brief The header file for the class swigTom::swigFreqWeight.
 * \author  Thomas Hrabe
 * \version 0.3
 * \date    30.7.2010
 *******************************************************************************/

#ifndef SWIGSINGLEAXISWEDGE_HPP_
#define SWIGSINGLEAXISWEDGE_HPP_

#include <tom/FreqWeight.hpp>
#include <swigVolume.hpp>
#include <fftw3.h>
namespace swigTom{

template<typename T>
class swigFreqWeight{

public:
	swigFreqWeight(double angle, double cutoff_radius,const std::size_t x,const std::size_t y,const std::size_t z){
		isAllPass = false;
		//bloody radians conversion
		angle = tom::math::deg2rad(angle);
		if(angle != 0.0){
			tom::FreqWeight_SingleAxisWedge<T>* fw = new tom::FreqWeight_SingleAxisWedge<T>(x,y,z,angle,cutoff_radius);
			this->fwPointer.reset(fw);
		}
		else{
			tom::FreqWeight_AllPass<T>* fw = new tom::FreqWeight_AllPass<T>(x,y,z);
			this->fwPointer.reset(fw);
			this->isAllPass = true;
		}

	};

	swigFreqWeight(double angle1, double angle2, double cutoff_radius,const std::size_t x,const std::size_t y,const std::size_t z,double smooth){
		isAllPass = false;
		//bloody radians conversion
		angle1 = tom::math::deg2rad(angle1);
		angle2 = tom::math::deg2rad(angle2);

		if(angle1 > 0.0 && angle2 > 0.0){
			tom::FreqWeight_SmoothedAsymmetricSingleAxisWedge<T>* fw = new tom::FreqWeight_SmoothedAsymmetricSingleAxisWedge<T>(x,y,z,angle1,angle2,cutoff_radius,smooth);
			this->fwPointer.reset(fw);
		}
		else{
			tom::FreqWeight_AllPass<T>* fw = new tom::FreqWeight_AllPass<T>(x,y,z);
			this->fwPointer.reset(fw);
			this->isAllPass = true;
		}

	};

	swigFreqWeight(swigVolume<T,T>& vol){
		isAllPass = false;
		tom::FreqWeight_Volume<T>* fw = new tom::FreqWeight_Volume<T>(vol,true);
		this->fwPointer.reset(fw);
	};

	swigFreqWeight(const std::string &filename, const std::size_t x,const std::size_t y,const std::size_t z){
		isAllPass = false;
		tom::FreqWeight_EmFile<T>* fw= new tom::FreqWeight_EmFile<T>(filename,x,y,z,true,true);
		this->fwPointer.reset(fw);
	};

	swigFreqWeight(const float lowestFrequency, const float highestFrequency,const std::size_t x,const std::size_t y,const std::size_t z,const T smooth){
		isAllPass = false;
		tom::FreqWeight_Bandpass<T>* fw = new tom::FreqWeight_Bandpass<T>(lowestFrequency,highestFrequency,x,y,z,smooth);
		this->fwPointer.reset(fw);
	}

	void rotate(double phi, double psi, double theta){
		phi = tom::math::deg2rad(phi);
		psi = tom::math::deg2rad(psi);
		theta = tom::math::deg2rad(theta);
		if(!this->isAllPass)
			this->fwPointer.get()->rotate(phi,psi,theta);
	};

	bool apply(swigVolume<std::complex<T>,T>& vsrc){
		if(!this->isAllPass){
			tom::FreqWeight<T>* fwP = this->fwPointer.get();
			return fwP->weight(vsrc);
		}
		else
			return false;
	};

	swigVolume<T,T> getWeightVolume(const bool reducedComplex){
		/*
		 * reducedComplex = if true, the volume will be scaled down to x,y,z/2+1.
		 * Thus, you must always provide the ORIGINAL size of the volume you are about to filter.
		 *
		 * */
		tom::FreqWeight<T>* fwP = this->fwPointer.get();
		tom::Volume<T>* v = const_cast<tom::Volume<T> * >(fwP->get_weight(reducedComplex));
		swigVolume<T,T> s(*(v));
		return s;
	};



private:
	std::auto_ptr<tom::FreqWeight<T> > fwPointer;
	bool isAllPass;
};

}


#endif /* SWIGSINGLEAXISWEDGE_HPP_ */
