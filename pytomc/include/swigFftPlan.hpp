/****************************************************************************//**
 * \file swigFftPlan.hpp
 * \brief The header file for the class swigTom::swigFftPlan.
 * \author  Thomas Hrabe
 * \version 0.2
 * \date    3.12.2008
 *******************************************************************************/


#ifndef SWIGFFTPLAN_HPP_
#define SWIGFFTPLAN_HPP_

#include <swigVolume.hpp>
#include <tom/fftw/fftw_plan.hpp>
#include <fftw3.h>

namespace swigTom{
/**********************************************************************//**
 *	\brief Plan class used for wraping tom::fftw::plan<T>
 *
 *	This class wraps tom::fftw::plan<T> to save either a r2c or c2r plan for two swigVolumes<T>.
 *	swigFftPlan<T> will stick with the two volumes and is only applicable to those.Furthermore, it stores pointers to the volumes internally.
 *	Make sure the plan is destroyed when one of each volumes is destroyed because the program will yield a memory leak otherwise.
 *
 */
template<typename T,typename TSCALE_SHIFT>
class swigFftPlan{
public:
		/**
		 *	\brief Wrapper constructor
		 */
		swigFftPlan(swigVolume<T,TSCALE_SHIFT> &v,swigVolume<std::complex<T>,TSCALE_SHIFT > &cv){
			tom::Volume<T>* pv = new tom::Volume<T>(&v.get(), v.sizeX(), v.sizeY(), v.sizeZ(), v.strideX(), v.strideY(), v.strideZ(), false, NULL);
			this->timeDomainP.reset(pv);
			tom::Volume<std::complex<T> >* pcv = new tom::Volume<std::complex<T> >(&cv.get(), cv.sizeX(), cv.sizeY(), cv.sizeZ(), cv.strideX(), cv.strideY(), cv.strideZ(), false, NULL);
			this->frequencyDomainP.reset(pcv);
			tom::fftw::Plan<T>* pp = new tom::fftw::Plan<T>(*this->timeDomainP.get(),*this->frequencyDomainP.get(), FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);
			this->r2cPlanP.reset(pp);
			this->isInverse = false;
		};
		/**
		 *	\brief Wrapper constructor
		 */
		swigFftPlan(swigVolume<std::complex<T>,TSCALE_SHIFT > &cv,swigVolume<T,TSCALE_SHIFT> &v){
			tom::Volume<T>* pv = new tom::Volume<T>(&v.get(), v.sizeX(), v.sizeY(), v.sizeZ(), v.strideX(), v.strideY(), v.strideZ(), false, NULL);
			this->timeDomainP.reset(pv);
			tom::Volume<std::complex<T> >* pcv = new tom::Volume<std::complex<T> >(&cv.get(), cv.sizeX(), cv.sizeY(), cv.sizeZ(), cv.strideX(), cv.strideY(), cv.strideZ(), false, NULL);
			this->frequencyDomainP.reset(pcv);
			tom::fftw::Plan<T>* pp = new tom::fftw::Plan<T>(*this->frequencyDomainP.get(),*this->timeDomainP.get(), FFTW_ESTIMATE);
			this->c2rPlanP.reset(pp);
			this->isInverse = true;
		};
		/**
		 *	\brief Wrapper function
		 */
		void transform(){
			if(!this->isInverse)
				this->r2cPlanP.get()->execute(*this->timeDomainP.get(),*this->frequencyDomainP.get());
			else
				this->c2rPlanP.get()->execute(*this->frequencyDomainP.get(),*this->timeDomainP.get());
		};

private:
		bool isInverse;
		std::auto_ptr<tom::Volume<T> > timeDomainP;
		std::auto_ptr<tom::Volume<std::complex<T> > > frequencyDomainP;
		std::auto_ptr<tom::fftw::Plan<T> > r2cPlanP;
		std::auto_ptr<tom::fftw::Plan<T> > c2rPlanP;
};
}

#endif /* SWIGFFTPLAN_HPP_ */
