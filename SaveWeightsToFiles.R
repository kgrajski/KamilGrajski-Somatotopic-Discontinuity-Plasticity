

	#
	#	Checking function for the C++
	#	Dump w1.e.0, w1.e.e, w1.i.e, w1.e.i.
	#

tmp = w1.e.0[,,1]; tmp = as.vector(tmp[]);
fName = "w1e0.dat1";
towrite = file ( fName, "wb" );
writeBin ( tmp, towrite, size=8, endian="little" );
close ( towrite );

tmp = w1.e.e[,,1]; tmp = as.vector(tmp[]);
fName = "w1ee.dat1";
towrite = file ( fName, "wb" );
writeBin ( tmp, towrite, size=8 );
close ( towrite );

tmp = w1.e.i[,,1]; tmp = as.vector(tmp[]);
fName = "w1ei.dat1";
towrite = file ( fName, "wb" );
writeBin ( tmp, towrite, size=8 );
close ( towrite );

tmp = w1.i.e[,,1]; tmp = as.vector(tmp[]);
fName = "w1ie.dat1";
towrite = file ( fName, "wb" );
writeBin ( tmp, towrite, size=8 );
close ( towrite );

