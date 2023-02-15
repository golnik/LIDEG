/*
 * File:   uBLASBoostOdeint.hpp
 * Author: Nikolay Golubev
 *
 * Created on September 25, 2015, 6:35 PM
 */

#ifndef UBLASBOOSTODEINT_HPP
#define	UBLASBOOSTODEINT_HPP

#include <boost/numeric/ublas/matrix.hpp>

#include <boost/numeric/odeint.hpp>

#include <boost/numeric/odeint/algebra/vector_space_algebra.hpp>

namespace boost{
namespace numeric{
namespace odeint{
    template<typename T>
    struct vector_space_norm_inf<boost::numeric::ublas::matrix<T>>{
	typedef typename T::value_type result_type;
	typedef boost::numeric::ublas::matrix<T> state_t;

	result_type operator()(const state_t& x) const{
	    result_type maxElement=0.;

	    auto comparator=[](const T& x1, const T& x2){
		return std::abs(x1)<std::abs(x2);
	    };

	    std::vector<T> elements;
	    for(auto it1=x.begin1(); it1!=x.end1(); it1++){
		for(auto it2=it1.begin(); it2!=it1.end(); it2++){
		    elements.push_back(*it2);
		}
	    }

	    auto result=std::max_element(elements.begin(),elements.end(),comparator);
	    maxElement=std::abs(*result);

	    return maxElement;
	}
    };

    template<typename T>
    struct vector_space_norm_inf<boost::numeric::ublas::vector<T>>{
	typedef typename T::value_type result_type;
	typedef boost::numeric::ublas::vector<T> state_t;

	result_type operator()(const state_t& x) const{
	    result_type maxElement=0.;

	    auto comparator=[](const T& x1, const T& x2){
		return std::abs(x1)<std::abs(x2);
	    };

	    auto result=std::max_element(x.begin(),x.end(),comparator);
	    maxElement=std::abs(*result);

	    return maxElement;
	}
    };
}
}
}

#endif	/* UBLASBOOSTODEINT_HPP */

