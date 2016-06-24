/*
 * rxnClassifier.hh
 *
 *  Created on: Feb 13, 2011
 *      Author: Leonard Harris
 */

#ifndef RXNCLASSIFIER_HH_
#define RXNCLASSIFIER_HH_

#include "../../std_include.hh"
#include "../../model/reaction.hh"

namespace network3{

	class RxnClassifier{
	public:
		enum{EXACT_STOCHASTIC,POISSON,LANGEVIN,DETERMINISTIC};
		RxnClassifier();
		RxnClassifier(const RxnClassifier& rc);
		virtual ~RxnClassifier();
		virtual void classifyRxns(vector<int>& classif, double tau, bool reclassify_all) = 0;
//		virtual RxnClassifier* clone() const{ return new RxnClassifier(*this); }
		void forceClassifications(int force){ this->force = force; }
	protected:
		int force;
	};
}
#endif /* RXNCLASSIFIER_HH_ */
