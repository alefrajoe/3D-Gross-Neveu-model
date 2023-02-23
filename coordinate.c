#include "coordinate.h"

//////////////////////////////////////////////////
//             coordinate stuff
//////////////////////////////////////////////////
/**
 * Print to stdout the DIM coordinates of the array : array dimensions[LSPACE, LSPACE, .. , LTIME].
 */
void PrintCoordinate(const int * coor_array)
{
    printf("Array is ");

    for(int d=0; d<DIM; d++)
    printf("[%d]", coor_array[d]);

    printf("\n");
}

/**
 * The function takes a number in [0, VOLUME) and converts it into a DIM-component coordinate.
 * The three coordinates are saved into the coor_array.
 * Note that the TIME-direction (which is usually the different one) is the last entry of the array.
 * Indeed in Euclidean notations the time direction is usually the last entry.
 */
void ConvertLessicographicIntoCoordinate(const int number, int * coor_array)
{
    // initialize the starting base and copy the number
    int base = LTIME, temp = number;

    // first step is the different one
    coor_array[DIM - 1] = (temp % LTIME); 
    temp = (temp / LTIME);

    // convert the last (DIM - 1) components into coordinates 
    for(int d=1; d<DIM; d++)
    {
        coor_array[DIM - 1 - d] = (temp % LSPACE);
        temp = (temp / LSPACE); 
    }
}

/**
 * The function takes an array of length DIM containing a lattice coordinate and converts it into a Lessicographic coordinate.
 * The new Lessicographic coordinate is saved into *number.
 * Remember that the last lattice coordinate represents the lattice coordinate along the temporal direction.
 */
void ConvertCoordinateIntoLessicographic(const int const * array, int * number)
{
    // initialize the base
    int base = 1, temp = 0;

    // cycle to reconstruct the base
    for(int d=0; d<DIM; d++)
    {
        // convert from lattice coordinate to Lessicographic number
        temp += array[DIM - 1 - d] * base;
        
        // change the base properly
        if(d == 0) base *= LTIME;
        else base *= LSPACE;
    }

    // save the result into time
    (*number) = temp;
}

/**
 * A site in Lessicographic basis is passed to the function.
 * The function computes the lattice sites (in Lessicographic coordinate), which are neighbors of the site passed, in the forward directions
 * so that in the array we save the coordinates of the sites passed.
 * The order in the array depends on the lattice directions.
 */
void InitializePPArray(int site, int *destination)
{
    // temporary array required to save the site in lattice coordinates
    int temp[DIM];

    // convert the site into lattice coordinates
    ConvertLessicographicIntoCoordinate(site, temp);

    // for all the directions 
    for(int d=0; d<DIM; d++)
    {
        // add 1 mod LSITE or LSPACE, depending on the direction
        if(d < (DIM - 1)) temp[d] = (temp[d] + 1) % LSPACE;
        else temp[d] = (temp[d] + 1) % LTIME;

        // convert the number and save it into the array
        ConvertCoordinateIntoLessicographic(temp, &(destination[d]));

        // remove 1 mod LSITE or LSPACE, depending on the direction
        if(d < (DIM - 1)) temp[d] = (temp[d] + LSPACE - 1) % LSPACE;
        else temp[d] = (temp[d] + LTIME - 1) % LTIME;
    }
}

/**
 * A site in Lessicographic basis is passed to the function.
 * The function computes the lattice sites (in Lessicographic coordinate), which are neighbors of the site passed, in the backward direction
 * so that in the array we save the coordinates of the sites passed.
 * The order in the array depends on the lattice directions.
 */
void InitializeMMArray(int site, int * destination)
{
    // temporary array required to save the site in lattice coordinates
    int temp[DIM];

    // convert the site into lattice coordinates
    ConvertLessicographicIntoCoordinate(site, temp);

    // for all the directions 
    for(int d=0; d<DIM; d++)
    {
        // subtract 1 mod LSITE or LSPACE, depending on the direction
        if(d < (DIM - 1)) temp[d] = (temp[d] + LSPACE - 1) % LSPACE;
        else temp[d] = (temp[d] + LTIME - 1) % LTIME;

        // convert the number and save it into the array
        ConvertCoordinateIntoLessicographic(temp, &(destination[d]));

        // remove 1 mod LSITE or LSPACE, depending on the direction
        if(d < (DIM - 1)) temp[d] = (temp[d] + 1) % LSPACE;
        else temp[d] = (temp[d] + 1) % LTIME;
    }
}

/**
 * Start the eta_pp-phase according to the parity of the passed site.
 * The site passed is in Lessicographic coordinates.
 * Eta is defined as: eta_\mu (x) = (-)^{n_1 + .. + n_{\mu-1}}.
 * The function also initialize the epsilon phases, being the pariti associated with each site \epsilon(x) = (-1)^{x_1+x_2+x_3}
 */
void StartEtapp(int site, double * destination)
{
    // temporary lattice coordinates
    int temp[DIM];

    // convert the site into lattice coordinates
    ConvertLessicographicIntoCoordinate(site, temp);

    // fill the eta_pp matrix
    for(int d=0; d<DIM; d++)
    {
        // calculate \sum_{\nu=1}^{\mu-1} n_\nu
        int sum = 0.0;
        for(int t=0; t<d; t++) sum += temp[t];
        // save into the destination array
        destination[d] = pow((-1.0), sum);

        // if APBC for the fermion field must be considered
        #ifdef APBC_FERMION_TIME
        // along time direction, if it is the last site along the time direction
        if(d == DIM-1)
        if(temp[(DIM-1)] == (LTIME - 1))
        destination[d] = -destination[d];
        #endif
        #ifdef APBC_ALL
        // along all directions, change sign
        if(temp[d] == (LTIME - 1))
        destination[d] = -destination[d];
        #endif
    }
}

/**
 * Start the eta_pp-phase according to the parity of the passed site.
 * The site passed is in Lessicographic coordinates.
 * Eta is defined as: eta_\mu (x) = (-)^{n_1 + .. + n_{\mu-1}}
 */
void StartEtamm(int site, double * destination)
{
    // temporary lattice coordinates
    int temp[DIM];

    // convert the site into lattice coordinates
    ConvertLessicographicIntoCoordinate(site, temp);

    // fill the eta_pp matrix
    for(int d=0; d<DIM; d++)
    {
        // calculate \sum_{\nu=1}^{\mu-1} n_\nu
        int sum = 0.0;
        for(int t=0; t<d; t++) sum += temp[t];
        // save into the destination array
        destination[d] = pow((-1.0), sum);

        // if APBC for the fermion field must be considered
        #ifdef APBC_FERMION_TIME
        // along time direction, if it is the first site along the time direction
        if(d == (DIM-1))
        if(temp[(DIM-1)] == 0)
        destination[d] = -destination[d];
        #endif
        #ifdef APBC_ALL
        // along all directions
        if(temp[d] == 0)
        destination[d] = -destination[d];
        #endif
    }
}

/**
 * Start the eta_pp-phase according to the parity of the passed site.
 * The site passed is in Lessicographic coordinates.
 * Eta is defined as: eta_\mu (x) = (-)^{n_1 + .. + n_{\mu-1}}
 */
void StartEtapbc(int site, double * destination)
{
    // temporary lattice coordinates
    int temp[DIM];

    // convert the site into lattice coordinates
    ConvertLessicographicIntoCoordinate(site, temp);

    // fill the eta_pp matrix
    for(int d=0; d<DIM; d++)
    {
        // calculate \sum_{\nu=1}^{\mu-1} n_\nu
        int sum = 0.0;
        for(int t=0; t<d; t++) sum += temp[t];
        // save into the destination array
        destination[d] = pow((-1.0), sum);
    }
}

/**
 * Compute the epsilon(x) phase \epsilon(x)=(-1.0)^{x_1+..+x_3} of all lattice sites.
 * The lattice site is passed into Lessicographic coordinates.
 */
void StartEpsilonPi(int site, int fermion_site, double * destination)
{
    // temporary lattice coordinates
    int temp[DIM];

    // convert the site into lattice coordinates
    ConvertLessicographicIntoCoordinate(site, temp);

    // save the epsilon phase at the proper fermion_site position
    destination[fermion_site] = pow((-1.0), ParityCoordinate(temp));
}

///////////////////////////////////////////////////////////////
//                 coordinate fermions
///////////////////////////////////////////////////////////////
/**
 * Retrieve the parity of a coordinate DIM-length []...[] = (-1)^{n_1 + .. + n_D}.
 */
int ParityCoordinate(const int const * array);

/**
 * Convert a number from fermion Lessicographic into standard Lessicographic.
 */
void ConvertFermionLessicographicIntoLessicographic(const int fermion, int * stereo)
{
    // look for all elements in VOLUME
    // this will be the outcome of the function
    for(int x=0; x<VOLUME; x++)
    {
        // consider x as a Lessicographic coordinate and 
        // convert x in fermion Lessicographic (saved in temp)
        // if temp == fermion, this is the result
        int temp, array[DIM];
        ConvertLessicographicIntoCoordinate(x, array);
        ConvertCoordinateIntoFermionLessicographic(array, &temp);
        if(temp == fermion) (*stereo) = x;
    }
}

/**
 * Convert a site in fermionic Lessicographic coordinates into lattice coordinates.
 * The function exploits some functions already defined in standard coordinates.
 * The lattice coordinates saved as : -SPACE : 1 , .. , LSPACE
 *                                    -TIME  : 1 , .. , LTIME   
 */
void ConvertFermionLessicographicIntoCoordinate(const int number, int * coor_array)
{
    // first convert the fermion number into standard Lessicographic
    int standard;
    ConvertFermionLessicographicIntoLessicographic(number, &standard);
    // now use standard functions
    ConvertLessicographicIntoCoordinate(standard, coor_array);
}

/**
 * Convert a Lessicographic coordinate into a fermion Lessicographic one.
 */
void ConvertLessicographicIntoFermionLessicographic(const int stereo, int * fermion)
{
    // to convert a Lessicographic into a fermion Lessicographic
    // we need an intermediate step where the lattice coordinate
    // are indeed required
    int array[DIM];
    ConvertLessicographicIntoCoordinate(stereo, array);
    ConvertCoordinateIntoFermionLessicographic(array, fermion);
}

/**
 * Convert a number into a fermionic coordinate.
 * Note that this is a bijective mapping in VOLUME -> VOLUME coordinates.
 */
void ConvertCoordinateIntoFermionLessicographic(const int const * array, int * fermion)
{
    // temporary integer
    int stereo = 0;

    // first convert the number into Lessicographic
    ConvertCoordinateIntoLessicographic(array, &stereo);
    // consider the parity of stereo (which is better computed from array)
    int parity = ParityCoordinate(array);
    
    if(parity == 0) (*fermion) = stereo / 2;
    else (* fermion) = (stereo / 2) + HALFVOLUME;
}

/**
 * A fermionic site number is passed to the function, in order to creare the fermionpp[DIM] array.
 * The function computes the lattice sites (in fermionic Lessicographic coordinate) that
 * are neighbors of the site passed (in the forward direction).
 * In the array we save the coordinates of the fermionic sites obtained.
 * The order in the array depends on the lattice directions from 1 to DIM
 */
void InitializeFermionPPArray(int site, int *destination)
{
    // temporary array required to save the site in lattice coordinates
    int temp[DIM];

    // convert the site into lattice coordinates
    ConvertFermionLessicographicIntoCoordinate(site, temp);

    // for all the directions 
    for(int d=0; d<DIM; d++)
    {
        // add 1 mod LSITE or LSPACE, depending on the direction
        if(d < (DIM - 1)) temp[d] = (temp[d] + 1) % LSPACE;
        else temp[d] = (temp[d] + 1) % LTIME;

        // convert the number and save it into the array
        ConvertCoordinateIntoFermionLessicographic(temp, &(destination[d]));

        // remove 1 mod LSITE or LSPACE, depending on the direction
        if(d < (DIM - 1)) temp[d] = (temp[d] + LSPACE - 1) % LSPACE;
        else temp[d] = (temp[d] + LTIME - 1) % LTIME;
    }
}

/**
 * A site in Lessicographic basis is passed to the function.
 * The function computes the lattice sites (in Lessicographic coordinate), which are neighbors of the site passed, in the backward direction
 * so that in the array we save the coordinates of the sites passed.
 * The order in the array depends on the lattice directions.
 */
void InitializeFermionMMArray(int site, int * destination)
{
    // temporary array required to save the site in lattice coordinates
    int temp[DIM];

    // convert the site into lattice coordinates
    ConvertFermionLessicographicIntoCoordinate(site, temp);

    // for all the directions 
    for(int d=0; d<DIM; d++)
    {
        // subtract 1 mod LSITE or LSPACE, depending on the direction
        if(d < (DIM - 1)) temp[d] = (temp[d] + LSPACE - 1) % LSPACE;
        else temp[d] = (temp[d] + LTIME - 1) % LTIME;

        // convert the number and save it into the array
        ConvertCoordinateIntoFermionLessicographic(temp, &(destination[d]));

        // remove 1 mod LSITE or LSPACE, depending on the direction
        if(d < (DIM - 1)) temp[d] = (temp[d] + 1) % LSPACE;
        else temp[d] = (temp[d] + 1) % LTIME;
    }
}

/**
 * A fermionic site number is passed to the function, in order to creare the fermionpp[DIM] array.
 * The function computes the lattice sites (in fermionic Lessicographic coordinate) that
 * are neighbors of the site passed (in the forward direction).
 * In the array we save the coordinates of the fermionic sites obtained.
 * The order in the array depends on the lattice directions from 1 to DIM
 */
void InitializeFastFermionPPArray(int site, int *destination, const int const * fermionsite)
{
    // temporary array required to save the site in lattice coordinates
    int temp[DIM];

    // convert the site into lattice coordinates
    ConvertLessicographicIntoCoordinate(site, temp);

    // for all the directions 
    for(int d=0; d<DIM; d++)
    {
        // add 1 mod LSITE or LSPACE, depending on the direction
        if(d < (DIM - 1)) temp[d] = (temp[d] + 1) % LSPACE;
        else temp[d] = (temp[d] + 1) % LTIME;

        // convert the number and save it into the array
        int prova = 0;
        ConvertCoordinateIntoLessicographic(temp, &prova);
        destination[d] = fermionsite[prova];

        // remove 1 mod LSITE or LSPACE, depending on the direction
        if(d < (DIM - 1)) temp[d] = (temp[d] + LSPACE - 1) % LSPACE;
        else temp[d] = (temp[d] + LTIME - 1) % LTIME;
    }
}