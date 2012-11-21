/** Reads line data from Robert Kurucz files called gf<elcode>.all, and saves
 * them in HDF format. */
int make_linelist_hdf(const char* filename);

/** Data structure containing all data from the Kurucz 160-character linelists.
 * */
struct LineStruct
{
    double wl;    /** Wavelength (nm) */
    double loggf; /** Log10 (gf) (RHS evaluated in CGS units) */
    double elcode; /** Kurucz element code (\f$ = Z * 100 + I \f$) */
    double energy1; /** Upper level energy (cm^{-1}) */
    double energy2; /** Lower level energy (cm^{-1}) */
    double J1; /** Angular momentum of upper level */
    double J2; /** Angular momentum of lower level */
};
