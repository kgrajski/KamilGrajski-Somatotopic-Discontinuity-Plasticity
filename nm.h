//	nm.h

//	Kamil A. Grajski
//	2-November-2014
//

//
//	Header file for the neural network model.
//

#ifndef NM_H
#define NM_H

#include <math.h>

#include <iomanip>
#include <locale>
#include <sstream>
#include <string>
#include <math.h>

#include "ioutil.h"
#include "simmgr.h"
#include "ran01.h"
#include "inpatch.h"
#include "wsparse.h"


	//////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////
	//
	//  START: NM class declarations.
	//
	//////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////

class NM {
public:
		NM () : v0(0), r0(0), v1E(0), r1E(0), v1I(0), r1I(0), r1ERFMapRaw(0), r1IRFMapRaw(0), r1ERFMapExpRes(0), r1IRFMapExpRes(0), w1E0(0), w1IE(0), w1EE(0), w1EI(0),
				v0CBuff(0), r0CBuff(0), v1ECBuff(0), r1ECBuff(0), v1ICBuff(0), r1ICBuff(0), w1E0CBuff(0), w1IECBuff(0), w1EECBuff(0), w1EICBuff(0) {};
		~NM () {};		

			//
			//	Methods, Operators, etc.
			//
							
		void AddNoise ( vector<double> &, const long &, const long &, const double &, const double & );

		void AddVecScalarMult ( vector<double> &, long &, vector <double> &, long &, double &, long & );

		void AddVecMatMult ( vector<double> &, const long &, vector<double> &, const long &, vector<double> &, const long &, const long & );

		void SetVecMatMult(vector<double> &, const long &, vector<double> &, const long &, vector<double> &, const long &, const long &);

		void BaselineRefinement ();

		void CleanUpDig3RefineInputPatchList ( vector<InputPatch> & );

		void CleanUpDig3SyndactylyInputPatchList ( vector<InputPatch> & );

		void CompareNE ( vector<double> &, vector<double> &, const double &, const long &, const long & );

		void CopyLastIterValsToZeroIter ();

		void CopyLastIterWeightsToZeroIter ();
		
		void GenDig3RefineInputPatchList ( vector<InputPatch> & );

		void GenDig3SelectiveStimInputPatchList ( vector<InputPatch> & );

		void GenDig3SyndactylyInputPatchList ( vector<InputPatch> & );

		void GenDig3SyndactylyInputPatchListSparse3 ( vector<InputPatch> & );

		void GenDig3SyndactylyControlInputPatchList ( vector<InputPatch> & );

		void GenUniformRefineInputPatchList ( vector<InputPatch> & );

		long GetCol ( const long &, const long & );

		long GetLegal ( const double &, const double &, const double & );

		long GetLin  ( const long &, const long &, const long & );

		long GetRow ( const long &, const long & );

		void InitializeNetworkAndWorkspace ();

		void InitPreSynWeights ( vector<double> &, const long &, const double &, const double &, const long &, const long &  );

		void IterateSystem4 ( const long & );

		void IterateSystem5 ( const long & );

		void LoadExpParms ( const string & );

		void LoadNetworkFromFile ( const string &, const long &, const long & );

		void NormalizePreSynWeights ( vector<double> &, const double &, const long &, const long & );

		void ParkLastIter ();

		void ReceptiveFieldMap ();

		void ReceptiveFieldProbeAnalysis ();

		void RandomPermutation ( vector<long> &, const long &, const long & );

		void RK4 ( vector<double> &, long &, vector <double> &, long &, double &, long &  );

		void SetNoise ( vector<double> &, const long &, const long &, const double &, const double & );

		void SaveNetworkToFile ( const string &, const long &, const long & );

		void SaveRFExpResToFile ( const string &, const long &, const long & );

		void SaveRFExpResToFileFocalStim ( const string &, const long &, const long &, const long &, const long & );

		void Sigmoid ( vector<double> &, vector<double> &, const double &, const long &, const long & );

		double Sigmoid1 ( const double &, const double & );

		void Syndactyly ();

		void SyndactylyControl ();

		void StimRespMag ( vector<double> &, vector<double>::iterator & );

		void StimRespMag2 ( vector<double> &, vector<double>::iterator & );

		void SubtractVecMatMult ( vector<double> &, const long &, vector<double> &, const long &, vector<double> &, const long &, const long & );

		void SumByCol ( vector<double> &, vector<double> &, const long &, const long &, const long & );

		void SumByRow ( vector<double> &, vector<double> &, const long &, const long &, const long &, const long & );

		void SystemSpinUp ();

		void UnParkLastIterToFirstIter ();

		void VecAdd (vector<double> &, const long &, vector<double> &, const long &, vector<double> &, const long &, const long &);

		void VecCopy ( vector<double> &, const long &, vector<double> &, const long &, const long & );

		double VecMaxValue ( vector<double> &, const long &, const long & );
		
		void VecScalarMult ( vector<double> &, long &, vector <double> &, long &, double &, long & );

		double VecSum ( vector<double> & );

		void VecVecMult ( vector<double> &, long &, vector<double> &, long &, vector<double> &, long &, const double &, const long & );

		void VecVecMultWeights ( vector<double> &, long &, vector<double> &, long &, vector<double> &, long &, long &, const double &, const long & );

					//
					//	Additional methods for handling sparse weight matrix representations.
					//	This is a very conservative approach to programming.  Keep all the old routines, because we know they work.
					//	Add the new functionality and double-check that it produces exactly the same results as before.
					//	Once verified, the original functions (see above) won't be maintained.
					//
		void AdaptSparseWeightMatrix ( vector<SparseWeightMatrix> &, const long &, const long &, const double &, const double &, vector<double>::iterator &, vector<double>::iterator & );
		
		void AdaptSparseWeightMatrixRK4 ( vector<SparseWeightMatrix> &, const long &, const long &, const double &, const double &, vector<double>::iterator &, vector<double>::iterator & );
		
		void AddInputPatchSparse ( long &, InputPatch &, vector<double>::iterator & );
		
		void AddVecMatMultSparse ( vector<double> &, const long &, vector<double> &, const long &, vector<SparseWeightMatrix> &, const long &, const long & );
		
		void BaselineRefinementSparse ();

		void BaselineRefinementSparse2 ( const string &, const long &, const long & );

		void ConvertToSparse ( vector<SparseWeightMatrix> &, vector<double> & );

		void GenDig3RefineInputPatchListSparse ( vector<InputPatch> & );

		void GenDig3RefineInputPatchListSparse2 (vector<InputPatch> &);

		void GenDig3RefineInputPatchListSparse3 ( vector<InputPatch> & );

		void GenDig3SelectiveStimInputPatchListSparse ( vector<InputPatch> & );

		void GenDig3SyndactylyInputPatchListSparse ( vector<InputPatch> & );

		void InitializeNetworkAndWorkspaceSparse ();
		
		void InitializeSparseWeightMatrix ( vector<SparseWeightMatrix> &, const long &, const long & );
		
		void IterateSystem5Sparse ( const long & );

		void IterateSystem5SparseNewNorm ( const long & );

		void IterateSystem5ASparseNewNorm ( const long & );

		void IterateSystem5SparseNewNormRK4 ( const long & );

		void IterateSystem5SparseUpdateECellWeightsOnly ( const long & );

		void LoadNetworkFromFileSparse ( const string &, const long &, const long & );
		
		void NormalizeSparseWeightMatrix ( vector<SparseWeightMatrix> &, const long &, const long &, const double & );

		void ReadSparseWeightMatrixFromFile ( const string &, vector<SparseWeightMatrix>::iterator &, vector<SparseWeightMatrix>::iterator & );
		
		void ReceptiveFieldMapSparse ();

		void ReceptiveFieldMapSparseStreamed ( const string &, const long &, const long & );

		void SaveNetworkToFileSparse ( const string &, const long &, const long & );

		void SaveNetworkToFileSparseFocalStim ( const string &, const long &, const long &, const long &, const long & );

		void SelectiveStimSparse ();
		
		void SubtractVecMatMultSparse ( vector<double> &, const long &, vector<double> &, const long &, vector<SparseWeightMatrix> &, const long &, const long & );
		
		void SyndactylySparse ();

		void SyndactylySparse2 ( const string &, const long &, const long & );

		void SyndactylyControlSparse ();

		void SystemSpinUpSparse ();

		void UniformRefinementSparse ();

		void WriteSparseWeightMatrixToFile ( const string &, vector<SparseWeightMatrix>::iterator &, vector<SparseWeightMatrix>::iterator & );

			//
			//		Constants
			//

		static const long	kWBetaDecayCycleTrigger = 10;
		static const long	kAllWeightsNetConfig = 5;
		static const long	kNoEEWeightsNetConfig = 4;

			//
			//		Variables
			//

		long		plasticityFlag;

		long		networkConfigType;
		long		experimentType;

		long		dumpRefineStepsRawToFile;
		long		dumpRefineStepsRawToFileStepSize;

		long		dumpRFMapRawToFile;
		long		doRFMapStepSize;

		string		baselineFileName;

		long		N;		
		long		g0;
		long		g1EE;
		long		g1EI;
		long		g1IE;

		double		noiseLevel;
		double		noiseLevelDotDeltaT;
		double		beta;
		double		tau;
		double		tauDivisor;
		
		double		wghtNormalizeFlag;
		double		eWResource;
		double		iWResource;
		double		eeWResource;		
		double		e0WResource;
		double		eiWResource;
		double		ieWResource;
		double		wghtMinValue;
		double		wghtMaxValue;
		double		wBeta;
		double		wBetaDecay;
		double		wBetaE0;
		double		wBetaEE;
		double		wBetaEI;
		double		wBetaIE;
		double		vETau;
		double		vITau;

		long		numSpinUpTrials;

		double		magRFProbe;
		double		durRFProbe;
		double		trialDurRFProbe;
		double		offsetToRFProbe;
		double		kRFPeakToEdgeDetect;

		long		numBaselineCycles;
		long		refinementPatchSize;
		double		refinementPatchNormalizedMagnitude;
		double		refinementOffsetToPatchStart;
		double		refinementStimDuration;

		long		numSyndactCycles;
		long		syndactPatchSize;
		double		syndactPatchNormalizedMagnitude;
		double		syndactOffsetToPatchStart;
		double		syndactStimDuration;
		long		syndactStimZoneID;

		long		numSelStimCycles;
		long		selectiveStimPatchSize;
		double		selectiveStimDuration;
		long		selectiveStimFactor;
		long		selectiveStimZoneID;

		long		numAmputCycles;
		long		refinementAmputPatchSize;
		double		refinementAmputStimDuration;
		long		kAmputDigit;
		long		kAmputZoneMax;

		long		numCLesionCycles;
		long		refinementCLesionPatchSize;
		double		refinementCLesionStimDuration;
		long		kCLesionDigit;
		long		kCLesionZoneMaxList;

					//	These are computed once all other parameters are known.

		double		deltaT;
		double		deltaTD2;
		double		deltaTD6;		//	Use with RK4
		double		minusOneOverTau;
		long		oneSecondNumIter;
		double		v0Tau;
		double		v0Alpha;
		double		v1EAlpha;
		double		v1IAlpha;
		double		wTau;
		double		wTauAlpha;
		double		minusOneOverWTau;

		long		numIterPerTrial;
		long		numIterPerTrialMinusOne;
		long		trialDurRFProbeInIters;
		long		offsetToRFProbeInIters;
		long		refinementOffsetToPatchInIters;
		long		syndactOffsetToPatchInIters;
		long		durRFProbeInIters;
		long		timeOffRFProbeInIters;

					//	Some useful intermediate variables

		long		N2;
		long		trialLengthInIters;
		long		numValsPerIter;
		long		numValsPerTrial;
		long		numValsPerRFMapExp;
		long		numValsPerRFMapRes;
		long		numWeightsPerIter;
		long		numWeightsPerTrial;
		double		wBetaCopy;
		double		wBetaE0Copy;
		double		wBetaEECopy;
		double		wBetaEICopy;
		double		wBetaIECopy;

					//	These are the time-series variables.

		vector<double>	vWork;
		vector<double>	v0;
		vector<double>	r0;
		vector<double>	v1E;
		vector<double>	r1E;
		vector<double>	v1I;
		vector<double>	r1I;
		vector<double>	r1ERFMapRaw;
		vector<double>	r1IRFMapRaw;
		vector<double>	r1ERFMapExpRes;
		vector<double>	r1IRFMapExpRes;

					//	These are the weight matrices.

		vector<double>	w1E0;
		vector<double>	w1IE;
		vector<double>	w1EE;
		vector<double>	w1EI;

					//	These are the matrices that support circular buffering.
		
		vector<double>	v0CBuff;
		vector<double>	r0CBuff;
		vector<double>	v1ECBuff;
		vector<double>	r1ECBuff;
		vector<double>	v1ICBuff;
		vector<double>	r1ICBuff;
		vector<double>	w1E0CBuff;
		vector<double>	w1IECBuff;
		vector<double>	w1EECBuff;
		vector<double>	w1EICBuff;

					//	Input Patches

		vector<InputPatch>	inputPatchList;
				
					//	These are the weight matrices - sparse representation.

		vector<SparseWeightMatrix>	w1E0_Sparse;
		vector<SparseWeightMatrix>	w1IE_Sparse;
		vector<SparseWeightMatrix>	w1EE_Sparse;
		vector<SparseWeightMatrix>	w1EI_Sparse;

					//	Preallocate for use in RK4.

		vector<double>	k1;
		vector<double>	k2;
		vector<double>	k3;
		vector<double>	k4;

private:

}; // class NM


	//////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////
	//
	//  End: NM class declarations.
	//
	//////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////

#endif