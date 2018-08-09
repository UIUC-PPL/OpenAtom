// compute the short range Hartree and e-N energies and ion forces
// EeNshort = \sum_J \sum_f wf rho_J(rf) \sum_<K!=J> eQ_K (erfc(alpha|rf-R_KJ|))/(|rf-R_KJ|) +
//						EeNselfshort
// EHarshort = \sum_J \sum_f wf rho_J(rf) \sum_<K!=J>\sum_f' e^2/2 rho_K(rf') (erfc(alpha|rf-rf'-R_JK|))/(|rf-rf'-R_JK|)  +
//						EHarselfshort
// This routine gets passed rho_J(rf), rho_K(rf) for all nearest neighbors of J, only a subset of J arrives (myAtmStart to 
//		myAtmEnd)
