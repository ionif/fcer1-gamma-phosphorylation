/*
 * eRungeKutta_TC_FG_rbPL.cpp
 *
 *  Created on: May 10, 2011
 *      Author: Leonard Harris
 */

#include "eRungeKutta_EXTRA.hh"
/*
eRungeKuttaRB_TC_FG_PL::eRungeKuttaRB_TC_FG_PL(){
	if (MoMMA::debug)
		cout << "eRungeKuttaRB_TC_FG_PL constructor called." << endl;
}
*/
eRungeKutta_TC_FG_rbPL::eRungeKutta_TC_FG_rbPL(ButcherTableau bt, double eps, double p, double pp, double q, double w,
		vector<SimpleSpecies*>& sp, vector<Reaction*>& rxn, Preleap_TC& ptc) : eRungeKutta_FG(bt,sp,rxn), p(p),
		pp(pp), q(q), w(w), preCalc(true), ptc(ptc), rxn(rxn){
	if (debug)
		cout << "eRungeKutta_TC_FG_rbPL constructor called." << endl;
	// Error check
	if (this->pp < this->p){
		cout << "Error in eRungeKutta_TC_FG_rbPL constructor: ";
		cout << "pp must be >= p; you have pp = " << this->pp << ", p = " << this->p << endl;
		exit(1);
	}
	if (this->q < 1.0){
		cout << "Error in eRungeKutta_TC_FG_rbPL constructor: ";
		cout << "q must be >= 1.0; your q = " << this->q << endl;
		exit(1);
	}
	if (this->w <= 0.0 || this->w >= 1.0){
		cout << "Error in eRungeKutta_TC_FG_rbPL constructor: ";
		cout << "w must be > 0.0 and < 1.0; your w = " << this->w << endl;
		exit(1);
	}
	this->ch = new RBChecker(eps,this->rxn);
	this->bc = new BinomialCorrector_RK(p,this->rxn);
	// Add rxns
	for (unsigned int v=0;v < this->rxn.size();v++){
		this->addRxn();
	}
}

eRungeKutta_TC_FG_rbPL::eRungeKutta_TC_FG_rbPL(ButcherTableau bt, double eps, double p, double pp, double q, double w,
		vector<SimpleSpecies*>& sp, vector<Reaction*>& rxn, Preleap_TC& ptc, bool round) : eRungeKutta_FG(bt,sp,rxn,round),
		p(p), pp(pp), q(q), w(w), preCalc(true), ptc(ptc), rxn(rxn){
	if (debug)
		cout << "eRungeKutta_TC_FG_rbPL constructor called." << endl;
	// Error check
	if (this->pp < this->p){
		cout << "Error in eRungeKutta_TC_FG_rbPL constructor: ";
		cout << "pp must be >= p; you have pp = " << this->pp << ", p = " << this->p << endl;
		exit(1);
	}
	if (this->q < 1.0){
		cout << "Error in eRungeKutta_TC_FG_rbPL constructor: ";
		cout << "q must be >= 1.0; your q = " << this->q << endl;
		exit(1);
	}
	if (this->w <= 0.0 || this->w >= 1.0){
		cout << "Error in eRungeKutta_TC_FG_rbPL constructor: ";
		cout << "w must be > 0.0 and < 1.0; your w = " << this->w << endl;
		exit(1);
	}
	this->ch = new RBChecker(eps,this->rxn);
	this->bc = new BinomialCorrector_RK(p,this->rxn);
	// Add rxns
	for (unsigned int v=0;v < this->rxn.size();v++){
		this->addRxn();
	}
}

eRungeKutta_TC_FG_rbPL::eRungeKutta_TC_FG_rbPL(const eRungeKutta_TC_FG_rbPL& tc_fg_pl) : eRungeKutta_FG(tc_fg_pl),
		p(tc_fg_pl.p), pp(tc_fg_pl.pp), q(tc_fg_pl.q), w(tc_fg_pl.w), preCalc(true), ptc(tc_fg_pl.ptc),
		rxn(tc_fg_pl.rxn){
	if (debug)
		cout << "eRungeKutta_TC_FG_rbPL copy constructor called." << endl;
	this->ch = new RBChecker(*tc_fg_pl.ch);
	this->bc = new BinomialCorrector_RK(*tc_fg_pl.bc);
	// Add rxns
	for (unsigned int v=0;v < this->rxn.size();v++){
		this->addRxn();
	}
}

eRungeKutta_TC_FG_rbPL::~eRungeKutta_TC_FG_rbPL(){
	if (debug)
		cout << "eRungeKutta_TC_FG_rbPL destructor called." << endl;
	delete this->ch;
	delete this->bc;
	for (unsigned int v=0;v < this->oldPop.size();v++){
		delete[] this->oldPop[v];
		delete[] this->projPop[v];
	}
}

void eRungeKutta_TC_FG_rbPL::getNewTau(double& tau){
	// Check for new rxns
	while (this->oldPop.size() != this->rxn.size() && this->projPop.size() != this->rxn.size()){
		this->addRxn();
	}
	// Get new tau
	if (this->preCalc){
		this->ptc.getNewTau(tau);
		this->preCalc = false;
	}
	else{
		if (this->substantially){ // Step was substantially accepted, increase tau
			tau *= this->q;
		}
		else{ // Step was barely accepted, reduce tau by a little bit
			tau *= this->pp;
		}
	}
	// Perform pre-check
	bool ok = false;
	while (!ok){
		//
		// Calculate a_eff[]
		this->aCalc->calc_aEff(tau);
		double mean_dX[this->aCalc->X_eff.size()];
		double sdev_dX[this->aCalc->X_eff.size()];
		//
		// Calculate projected species population changes
		for (unsigned int j=0;j < this->aCalc->X_eff.size();j++){
			mean_dX[j] = 0.0;
			sdev_dX[j] = 0.0;
			double z_vj;
			unsigned int R_v;
			for (unsigned int v=0;v < this->aCalc->spInRxn[j].size();v++){
				z_vj = this->aCalc->stoich[j][v];
				R_v = this->aCalc->spInRxn[j][v];
				mean_dX[j] += z_vj*this->aCalc->a_eff[R_v];
				sdev_dX[j] += z_vj*z_vj*this->aCalc->a_eff[R_v];
			}
			mean_dX[j] *= tau;
			sdev_dX[j] *= tau;
			sdev_dX[j] = sqrt(sdev_dX[j]);
			if (mean_dX[j] < 0.0){ // If the mean is negative, make the sdev negative
				sdev_dX[j] = -sdev_dX[j];
			}
		}
		// Calculate elements of projPop[][]
		unsigned int S_j;
		for (unsigned int v=0;v < this->projPop.size();v++){
			for (unsigned int j=0;j < this->rxn[v]->rateSpecies.size();j++){
				S_j = this->aCalc->rateSp[v][j];
				this->projPop[v][j] = this->oldPop[v][j] + mean_dX[S_j] + sdev_dX[S_j];
			}
		}
		// Check against current rates
		ok = this->ch->check(1.0,this->aCalc->a_eff,this->projPop,false);
		if (!ok){
			tau *= this->p; // Reduce
		}
	}
}

void eRungeKutta_TC_FG_rbPL::fireRxns(vector<double>& k, vector<int>& classif, double tau){
	// a_eff[] elements have already been calculated in getNewTau()
	this->fg->fireRxns(k,classif,tau,this->aCalc->a_eff);
}

bool eRungeKutta_TC_FG_rbPL::check(){
	// Check for new rxns
	while (this->oldPop.size() != this->rxn.size() && this->projPop.size() != this->rxn.size()){
		this->addRxn();
	}
//	cout << "**Checking**" << endl;
	bool ok;
	this->substantially = this->ch->check(this->w,this->aCalc->a_eff,this->oldPop,true);
	if (this->substantially){
		ok = true;
//		cout << "**Substantially accepted**" << endl;
	}
	else{
		ok = this->ch->check(1.0,this->aCalc->a_eff,this->oldPop,true);
//		if (ok) cout << "**Barely accepted**" << endl;
//		else cout << "**Rejected**" << endl;
	}
	return ok;
}

void eRungeKutta_TC_FG_rbPL::update(){
	// Update oldPop[][]
	for (unsigned int v=0;v < this->oldPop.size();v++){
		for (unsigned int j=0;j < this->rxn[v]->rateSpecies.size();j++){
			this->oldPop[v][j] = this->rxn[v]->rateSpecies[j]->population;
		}
	}
	// Just in case
	if (this->oldPop.size() != this->rxn.size()){
		cout << "Error in eRungeKutta_TC_FG_rbPL::update(): Sizes of 'oldPop' and 'rxn' vectors not equal. "
			 << "Shouldn't happen. Exiting." << endl;
		exit(1);
	}
	if (this->projPop.size() != this->rxn.size()){
		cout << "Error in eRungeKutta_TC_FG_rbPL::update(): Sizes of 'projPop' and 'rxn' vectors not equal. "
			 << "Shouldn't happen. Exiting." << endl;
		exit(1);
	}
}

void eRungeKutta_TC_FG_rbPL::addRxn(){
	if (this->oldPop.size() < this->rxn.size()  && this->projPop.size() < this->rxn.size()
			&& this->oldPop.size() == this->projPop.size()){
		unsigned int u = this->oldPop.size();
		this->oldPop.push_back(new double[this->rxn[u]->rateSpecies.size()]);
		this->projPop.push_back(new double[this->rxn[u]->rateSpecies.size()]);
		for (unsigned int j=0;j < this->rxn[u]->rateSpecies.size();j++){
			this->oldPop[u][j] = this->rxn[u]->rateSpecies[j]->population;
			this->projPop[u][j] = 0.0; // Initialize
		}
	}
	else{
		cout << "Error in eRungeKutta_TC_FG_rbPL::addRxn(): No rxns to add (oldPop.size = " << this->oldPop.size()
			 << ", projPop.size = " << this->projPop.size() << ", rxn.size = " << this->rxn.size()
			 << "). Shouldn't happen. Exiting." << endl;
		exit(1);
	}
}
