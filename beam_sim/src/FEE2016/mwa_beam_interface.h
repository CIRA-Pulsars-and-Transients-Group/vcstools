
// az - is assumed to be astro azimuth N=0, E=90, S=180, W=270 :
extern "C" double CalcMWABeam( double az, double za, double freq_hz, char beam, int gridpoint, int zenith_norm=1 );

double CalcMWABeam_FindClosest( double az, double za, double freq_hz, char beam='X', int zenith_norm=0, int find_closest=1, int speed_test =1 );
