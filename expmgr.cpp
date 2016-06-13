//	expmgr.cpp

//	Kamil A. Grajski
//	2-November-2014
//

#include "expmgr.h"
#include "ioutil.h"
#include "ran01.h"

	//////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////
	//
	//  START: ExpManager class definitions.
	//
	//////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////

ExpManager::ExpManager ( NM & nm ) {
	

			//
			//	SYNDACTYLY.  This section implements control experiments. 
			//

			//
			//	Initial random network.
			//
	nm.InitializeNetworkAndWorkspaceSparse ();
	nm.SaveNetworkToFileSparse("SysInit", nm.numBaselineCycles, 0);
	cout << "InitializeNetworkAndWorkspace.  Done." << endl;

			//
			//	Do a system spin-up to shake out transients and do an initial RF map.
			//
	nm.SystemSpinUpSparse();
	nm.SaveNetworkToFileSparse ( "Base", nm.numBaselineCycles, 0 ); cout << "SystemSpinUp Done." << endl;
	nm.ReceptiveFieldMapSparseStreamed ( "Base.RFMap.Raw", nm.numBaselineCycles, 0 ); cout << "ReceptiveFieldMap.  Done." << endl;
	nm.SaveRFExpResToFile ( "Base.RFMap", nm.numBaselineCycles, 0 ); cout << "ReceptiveFieldMap Files Written.  Done." << endl;

			//
			//	Part I: Baseline refinement.
			//
	for ( long iCycle = 1; iCycle <= nm.numBaselineCycles; ++iCycle ) {		
		nm.BaselineRefinementSparse2 ( "Base.RawPatch", nm.numBaselineCycles, iCycle );
		if ( ( nm.dumpRFMapRawToFile && ( (iCycle % nm.doRFMapStepSize) == 0 )) || ( iCycle == 1 ) ) {
			nm.SaveNetworkToFileSparse ( "Base", nm.numBaselineCycles, iCycle );
			nm.ReceptiveFieldMapSparseStreamed ( "Base.RFMap.Raw", nm.numBaselineCycles, iCycle );
			nm.SaveRFExpResToFile ( "Base.RFMap", nm.numBaselineCycles, iCycle ); cout << "BR & RF Map & Files Written.  Cycle: " << iCycle << " Done." << endl;
		} // if ( nm.doRFMapStepSize && ( (iCycle % nm.doRFMapStepSize) == 0 ) ) {		
		nm.wBeta *= nm.wBetaDecay; nm.wBetaE0 *= nm.wBetaDecay; nm.wBetaEE *= nm.wBetaDecay; nm.wBetaEI *= nm.wBetaDecay; nm.wBetaIE *= nm.wBetaDecay;
	} // 	for ( long iCycle = 1; iCycle <= nm.numBaselineCycles; ++iCycle ) {

			//
			//	Part II: Syndactyly refinement.
			//
	nm.wBeta = nm.wBetaCopy; nm.wBetaE0 = nm.wBetaE0Copy; nm.wBetaEE = nm.wBetaEECopy; nm.wBetaEI = nm.wBetaEICopy; nm.wBetaIE = nm.wBetaIECopy;
	for ( long iCycle = 1; iCycle <= nm.numBaselineCycles; ++iCycle ) {
		nm.SyndactylySparse2 ( "SyndactExp.RawPatch", nm.numBaselineCycles, iCycle );
		if ( ( nm.dumpRFMapRawToFile && ( (iCycle % nm.doRFMapStepSize) == 0 )) || ( iCycle == 1 ) ) {
			nm.SaveNetworkToFileSparse ( "SyndactExp", nm.numBaselineCycles, iCycle );
			nm.ReceptiveFieldMapSparseStreamed ( "SyndactExp.RFMap.Raw", nm.numBaselineCycles, iCycle );
			nm.SaveRFExpResToFile ( "SyndactExp.RFMap", nm.numBaselineCycles, iCycle ); cout << "SYN & RF Map & Files Written.  Cycle: " << iCycle << " Done." << endl;
		} // if (nm.dumpRFMapRawToFile && ((iCycle % nm.doRFMapStepSize) == 0)) {
		nm.wBeta *= nm.wBetaDecay; nm.wBetaE0 *= nm.wBetaDecay; nm.wBetaEE *= nm.wBetaDecay; nm.wBetaEI *= nm.wBetaDecay; nm.wBetaIE *= nm.wBetaDecay;
	} // 	for ( long iCycle = 1; iCycle <= nm.numSyndactCycles; ++iCycle ) {

			//
			//	Part III Syndactyly Release experiment.
			//
	nm.wBeta = nm.wBetaCopy; nm.wBetaE0 = nm.wBetaE0Copy; nm.wBetaEE = nm.wBetaEECopy; nm.wBetaEI = nm.wBetaEICopy; nm.wBetaIE = nm.wBetaIECopy;
	for ( long iCycle = 1; iCycle <= nm.numBaselineCycles; ++iCycle ) {															
		nm.BaselineRefinementSparse2 ( "SyndactReleaseExp.RawPatch", nm.numBaselineCycles, iCycle );
		if ( ( nm.dumpRFMapRawToFile && ( (iCycle % nm.doRFMapStepSize) == 0 )) || ( iCycle == 1 ) ) {
			nm.SaveNetworkToFileSparse ( "SyndactReleaseExp", nm.numBaselineCycles, iCycle );
			nm.ReceptiveFieldMapSparseStreamed ( "SyndactReleaseExp.RFMap.Raw", nm.numBaselineCycles, iCycle );
			nm.SaveRFExpResToFile ( "SyndactReleaseExp.RFMap", nm.numBaselineCycles, iCycle ); cout << "SYN_REL & RF Map & Files Written.  Cycle: " << iCycle << " Done." << endl;
		} // if (nm.dumpRFMapRawToFile && ((iCycle % nm.doRFMapStepSize) == 0)) {
		nm.wBeta *= nm.wBetaDecay; nm.wBetaE0 *= nm.wBetaDecay; nm.wBetaEE *= nm.wBetaDecay; nm.wBetaEI *= nm.wBetaDecay; nm.wBetaIE *= nm.wBetaDecay;
	} // 	for ( long iCycle = 1; iCycle <= nm.numBaselineCycles; ++iCycle ) {
			
			//
			//	Part IV: Syndactyly control.
			//
	nm.LoadNetworkFromFileSparse ( "Base", nm.numBaselineCycles, nm.numBaselineCycles );
	nm.wBeta = nm.wBetaCopy; nm.wBetaE0 = nm.wBetaE0Copy; nm.wBetaEE = nm.wBetaEECopy; nm.wBetaEI = nm.wBetaEICopy; nm.wBetaIE = nm.wBetaIECopy;
	for ( long iCycle = 1; iCycle <= nm.numBaselineCycles; ++iCycle ) {															
		nm.BaselineRefinementSparse2 ( "SyndactCtl.RawPatch", nm.numBaselineCycles, iCycle );
		if ( ( nm.dumpRFMapRawToFile && ( (iCycle % nm.doRFMapStepSize) == 0 )) || ( iCycle == 1 ) ) {
			nm.SaveNetworkToFileSparse ( "SyndactCtl", nm.numBaselineCycles, iCycle );
			nm.ReceptiveFieldMapSparseStreamed ( "SyndactCtl.RFMap.Raw", nm.numBaselineCycles, iCycle );
			nm.SaveRFExpResToFile ( "SyndactCtl.RFMap", nm.numBaselineCycles, iCycle ); cout << "SYN_CTL & RF Map & Files Written.  Cycle: " << iCycle << " Done." << endl;
		} // if (nm.dumpRFMapRawToFile && ((iCycle % nm.doRFMapStepSize) == 0)) {	
		nm.wBeta *= nm.wBetaDecay; nm.wBetaE0 *= nm.wBetaDecay; nm.wBetaEE *= nm.wBetaDecay; nm.wBetaEI *= nm.wBetaDecay; nm.wBetaIE *= nm.wBetaDecay;
	} // 	for ( long iCycle = 1; iCycle <= nm.numSyndactCycles; ++iCycle ) {
	
			//
			//	Part V: Syndactyly Release control.
			//
	nm.wBeta = nm.wBetaCopy; nm.wBetaE0 = nm.wBetaE0Copy; nm.wBetaEE = nm.wBetaEECopy; nm.wBetaEI = nm.wBetaEICopy; nm.wBetaIE = nm.wBetaIECopy;
	for ( long iCycle = 1; iCycle <= nm.numBaselineCycles; ++iCycle ) {															
		nm.BaselineRefinementSparse2 ( "SyndactReleaseCtl.RawPatch", nm.numBaselineCycles, iCycle );	
		if ( ( nm.dumpRFMapRawToFile && ( (iCycle % nm.doRFMapStepSize) == 0 )) || ( iCycle == 1 ) ) {
			nm.SaveNetworkToFileSparse ( "SyndactReleaseCtl", nm.numBaselineCycles, iCycle );
			nm.ReceptiveFieldMapSparseStreamed ( "SyndactReleaseCtl.RFMap.Raw", nm.numBaselineCycles, iCycle );
			nm.SaveRFExpResToFile ("SyndactReleaseCtl.RFMap", nm.numBaselineCycles, iCycle ); cout << "SYN_REL_CTL & RF Map & Files Written.  Cycle: " << iCycle << " Done." << endl;
		} // if (nm.dumpRFMapRawToFile && ((iCycle % nm.doRFMapStepSize) == 0)) {	
		nm.wBeta *= nm.wBetaDecay; nm.wBetaE0 *= nm.wBetaDecay; nm.wBetaEE *= nm.wBetaDecay; nm.wBetaEI *= nm.wBetaDecay; nm.wBetaIE *= nm.wBetaDecay;
	} // 	for ( long iCycle = 1; iCycle <= nm.numBaselineCycles; ++iCycle ) {
			
}; // ExpManager::ExpManager ( NM & nm ) {

	//////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////
	//
	//  END: ExpManager class definitions.
	//
	//////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////