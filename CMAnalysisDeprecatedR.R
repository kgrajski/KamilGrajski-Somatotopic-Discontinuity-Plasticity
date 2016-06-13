#
#	CM analysis scripts.  

#
#	These came from: CorticalMagnificationAnalysis.CPP.SelStim.R
#
#
#	These are here, because while they may have proved useful they were
#	stepping stones towards the mainline analysis methods.
#

		#
		#	Version 1.0: By sector.
		#
iControlSector = 2; iControlSectorCellList = GenD3SectorCellList ( N, iControlSector );
iExpSector = 2; iExpSectorCellList = GenD3SectorCellList ( N, iExpSector );

init.CMag.Fig.e = rep ( 0, N2 );
init.CMag.Fig.e[MapMagToCellList( init.cortical.Amp.e, iControlSectorCellList)] = 1;
init.CMag.Fig.e[MapMagToCellList( init.cortical.Amp.e, iExpSectorCellList)] = 1;

init.CMag.Fig.i = rep ( 0, N2 );
init.CMag.Fig.i[MapMagToCellList( init.cortical.Amp.i, iControlSectorCellList)] = 1;
init.CMag.Fig.i[MapMagToCellList( init.cortical.Amp.i, iExpSectorCellList)] = 1;

final.CMag.Fig.e = rep ( 0, N2 );
final.CMag.Fig.e[MapMagToCellList( final.cortical.Amp.e, iControlSectorCellList)] = 1;
final.CMag.Fig.e[MapMagToCellList( final.cortical.Amp.e, iExpSectorCellList)] = 1;

final.CMag.Fig.i = rep ( 0, N2 );
final.CMag.Fig.i[MapMagToCellList( final.cortical.Amp.i, iControlSectorCellList)] = 1;
final.CMag.Fig.i[MapMagToCellList( final.cortical.Amp.i, iExpSectorCellList)] = 1;

x11(); par(mfrow=c(2,2));
ShowVecAsMap ( init.CMag.Fig.e, "E" ); abline ( h = c ( (N/3)+0.5, 2*(N/3)+0.5 ) );
ShowVecAsMap ( init.CMag.Fig.i, "I" ); abline ( h = c ( (N/3)+0.5, 2*(N/3)+0.5 ) );
ShowVecAsMap ( final.CMag.Fig.e, "E" ); abline ( h = c ( (N/3)+0.5, 2*(N/3)+0.5 ) );
ShowVecAsMap ( final.CMag.Fig.i, "I" ); abline ( h = c ( (N/3)+0.5, 2*(N/3)+0.5 ) );

		#
		#	Version 1.1: By sector with frequency count cutoff.
		#
ctlVal = 25;
expVal = 200;
iControlSector = 6; iControlSectorCellList = GenD3SectorCellList ( N, iControlSector, 1 );
iExpSector = 4; iExpSectorCellList = GenD3SectorCellList ( N, iExpSector, 1 );

init.CMag.Fig.e = rep ( 0, N2 );
init.CMag.Fig.e[MapMagToCellList1GE( init.cortical.Amp.e, iControlSectorCellList, init.stimCount, ctlVal )] = 1;
init.CMag.Fig.e[MapMagToCellList1GE( init.cortical.Amp.e, iExpSectorCellList, final.stimCount, expVal )] = 1;

init.CMag.Fig.i = rep ( 0, N2 );
init.CMag.Fig.i[MapMagToCellList1GE( init.cortical.Amp.i, iControlSectorCellList, init.stimCount, ctlVal )] = 1;
init.CMag.Fig.i[MapMagToCellList1GE( init.cortical.Amp.i, iExpSectorCellList, final.stimCount, expVal )] = 1;

final.CMag.Fig.e = rep ( 0, N2 );
final.CMag.Fig.e[MapMagToCellList1GE( final.cortical.Amp.e, iControlSectorCellList, init.stimCount, ctlVal )] = 1;
final.CMag.Fig.e[MapMagToCellList1GE( final.cortical.Amp.e, iExpSectorCellList, final.stimCount, expVal )] = 1;

final.CMag.Fig.i = rep ( 0, N2 );
final.CMag.Fig.i[MapMagToCellList1GE( final.cortical.Amp.i, iControlSectorCellList, init.stimCount, ctlVal )] = 1;
final.CMag.Fig.i[MapMagToCellList1GE( final.cortical.Amp.i, iExpSectorCellList, final.stimCount, expVal )] = 1;

x11(); par(mfrow=c(2,2));
ShowVecAsMap ( init.CMag.Fig.e, "E" ); abline ( h = c ( (N/3)+0.5, 2*(N/3)+0.5 ) );
ShowVecAsMap ( init.CMag.Fig.i, "I" );
ShowVecAsMap ( final.CMag.Fig.e, "E" ); abline ( h = c ( (N/3)+0.5, 2*(N/3)+0.5 ) );
ShowVecAsMap ( final.CMag.Fig.i, "I" );

		#
		#	Version 1.2: By hand-picked input layer node locations.
		#
iControlSectorCellList = (N) * (N/3) + as.integer ( (N/3) );
iExpSectorCellList = (N) * (N/3) + as.integer ( (N/3) );
iControlSectorCellList = seq ( 21 * N + 1, 23 * N );
iExpSectorCellList = seq ( 10 * N + 1, 11 * N );

iControlSectorCellList = tmp;
iExpSectorCellList = tmp;

x11(); par(mfrow=c(2,2));

for ( i in 1:length(iControlSectorCellList ) ) {

	init.CMag.Fig.e = rep ( 0, N2 );
	init.CMag.Fig.e[MapMagToCellList( init.cortical.Amp.e, iControlSectorCellList[i])] = 1;
	init.CMag.Fig.e[MapMagToCellList( init.cortical.Amp.e, iExpSectorCellList[i])] = 1;

	init.CMag.Fig.i = rep ( 0, N2 );
	init.CMag.Fig.i[MapMagToCellList( init.cortical.Amp.i, iControlSectorCellList[i])] = 1;
	init.CMag.Fig.i[MapMagToCellList( init.cortical.Amp.i, iExpSectorCellList[i])] = 1;

	final.CMag.Fig.e = rep ( 0, N2 );
	final.CMag.Fig.e[MapMagToCellList( final.cortical.Amp.e, iControlSectorCellList[i])] = 1;
	final.CMag.Fig.e[MapMagToCellList( final.cortical.Amp.e, iExpSectorCellList[i])] = 1;

	final.CMag.Fig.i = rep ( 0, N2 );
	final.CMag.Fig.i[MapMagToCellList( final.cortical.Amp.i, iControlSectorCellList[i])] = 1;
	final.CMag.Fig.i[MapMagToCellList( final.cortical.Amp.i, iExpSectorCellList[i])] = 1;

	ShowVecAsMap ( init.CMag.Fig.e, "E" ); abline ( h = c ( (N/3)+0.5, 2*(N/3)+0.5 ) );
	ShowVecAsMap ( init.CMag.Fig.i, "I" ); abline ( h = c ( (N/3)+0.5, 2*(N/3)+0.5 ) );
	ShowVecAsMap ( final.CMag.Fig.e, "E" ); abline ( h = c ( (N/3)+0.5, 2*(N/3)+0.5 ) );
	ShowVecAsMap ( final.CMag.Fig.i, "I" ); abline ( h = c ( (N/3)+0.5, 2*(N/3)+0.5 ) );

	Sys.sleep ( 2 );

} # for ( i in 1:length(iControlSector) ) {

