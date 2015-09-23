


# define coulombconstant 332.06 // in units of kcal/mol, angstroms, electron charge
# define onedegreeinkcals 0.001987
# define roomTinkcals 0.5922
# define avogadrosnumberover10raised27 0.0006022 // raised to 10^27 to change molarity from per L^3 to per angstrom^3
# define waterpermittivity 80
# define molardensity 0.00060221415
# define heightequalitytolerance 0.000000001
# define maxrforvmmcbuffer 0.1
# define maxmaxneighbors 5000

typedef struct coord_tag{
	double_triple r;
	double_triple n;
	int nodetype;
	int chainid;
	int monomerid;
	int orientationtype;
    int leafid;
    double height;
} coord;

typedef struct sidechain_params_tag{
	
	//	either
	
	double k1;
	double rperp0;
	double rpar0;
	double r0;
	double J10;     //	add J200 and J300 to this
	double J11;
	double J12;
	double J13;
	double J14;
	
	//	perpendicular
	
	//double J200;	//	add J300 to this
	double J201;
	double J202;
	double J203;
	double J204;
	double J205;
	double J206;
	double J220;	//	add J302 to this
	double J221;
	double J222;
	double J240;	//	add J304 to this
	//double J300;
	//double J302;
	//double J304;	//	redundant
	double J311;
	double J313;
	double J320;
	double J322;
	double J331;
	double J340;
	
	//	parallel
	
	double k2;
	double r20;
	double r21;
	double r22;
	double k3;
	double r30;
	double r31;
	double r32;
} sidechain_params;

typedef struct bonded_params_tag{
	double K10;
	double K11;
	double K12;
	double K13;
	double K14;
	double K20;
	double K21;
	double K22;
	double K23;
	double K24;
	double kl;
	double r0l;
	double sl;
	double kr;
	double r0r;
	double sr;
	sidechain_params *sidechain;
	int *sidechainorientationtype;
    double factor;
} bonded_params;

typedef struct GBsolv_tag{		//	fixed mu=1, nu=-2
	double eps0;
	double sigma0;
	double chi;
	double chiprime;
	double chirep;
	double chiattr;
	double xi;
	double amprep;
	double ampattr;
	double w;
	double interpolatemiddle;
	double interpolatewidth;
} GBsolv;

typedef struct GBQvac_tag{		//	fixed mu=1, nu=-2
	double eps0;
	double sigma0;
	double chi;
	double chiprime;
	double xi;
	double Q;
	
	int code;
	GBsolv p11solv;
    double factor;
} GBQvac;

typedef struct electrostaticparam_tag{
	double rhard;
	double charge;
	double dipole;
	double solvationenergy;
    double shift;
    double factor;
} electrostaticparam;

typedef struct onebodyparam_tag{
	double solvationenergy;
	double z0;
	double uinterface;
	double zinterface;
	double sigmainterface;
} onebodyparam;

typedef struct solvation_parameters_tag{
	double shortrangepermittivity;
	double waterradius;
	double firstshellfraction;
	double saturationrange;
	double debyelength;
	int interface;
	double interfacethickness;
	double meshbondlength;
	double maxmeshbond;
	int_double meshsize;
	int meshnumber;
	double surfacetension;
    double Tolmanlength;
    double meancurvaturelength;
	double meshsurfacecutoff;
	double longestmeshbond;
	double meshcrumpleneighborfactor;
    int meshmovecode;
    double *totallength;
	double totalvolume;
	double targetvolume;
	double kvolumerestraint;
    double kuniformbonds;
    
    //  Two parameters for use in monolayer.c (fixed interface):
    
	double interfaceheight1;
	double interfaceheight2;
} solvation_parameters;

typedef struct nonbonded_params_tag{
	double cutoff2;
    double phenylcutoff2;
 	double rhard0;				//	backbone N
	double rhard1;				//	phenyl
	GBQvac p11vac;				//	phenyl-phenyl
	solvation_parameters solvationparams;
	electrostaticparam p2;		//	amino
	electrostaticparam p3;		//	carboxyl
	onebodyparam one0;
	onebodyparam one1;
	onebodyparam one2;
	onebodyparam one3;
	double vacuumthickness;
} nonbonded_params;

typedef struct monomernodes_tag{
	int backbone;
	int sidechain;
} monomernodes;

typedef struct movieparams_tag{
	double phenylspacing;
	double aminoNdepth;
	double aminoCdepth;
	double aminoCwidth;
	double carboxylbackCdepth;
	double carboxylforwardCdepth;
	double carboxylOdepth;
	double carboxylOwidth;
	double Crad;
	double Orad;
	double Nrad;
} movieparams;

typedef struct cgparams_tag{
	double backrad;
	double phenrad;
	double aminrad;
	double carbrad;
} cgparams;

typedef struct energycomponents_tag{
	double backbone;
	double sidechain;

	double nnsamepoly;
	double ccunlikesamepoly;
	double cclikesamepoly;
	double cnsamepoly;

	double nn;
	double ccunlike;
	double cclike;
	double cn;
	
	double nncross;
	double ccunlikecross;
	double cclikecross;
	double cncross;

	double nndifferent;
	double ccunlikedifferent;
	double cclikedifferent;
	double cndifferent;
    
    double *solvation;
	double *interface;
} energycomponents;

typedef struct monomertypes_tag{
	int charged;
	int nonpolar;
	int *type;
	int monomers;
} monomertypes;

typedef struct monolayer_tag{
    double_triple **backboneoffset;      //  (first one for last polar type (0 if self is polar), second one for self)
    double **backbonenz;
    double **backbonephi;
    double_triple **sidechainrelativepos;
    double **sidechainnpar;
    double **sidechainnphi;
    
	double monomerspacing;
	double interchainspacing;
	double terminusspacing;
	double offset;
	double boxheight;
	int phenylcode;
    
	double offsetfraction;
	double monomerscaleoffsetfraction;
	double length;	
} monolayer;

typedef struct monolayer_alt_tag{
    double_triple ***backboneoffset;      //  (first one for leaf, second one for last polar type (0 if self is polar), third one for self)
    double ***backbonenz;
    double ***backbonephi;
    double_triple ***sidechainrelativepos;
    double ***sidechainnpar;
    double ***sidechainnphi;
    
	double monomerspacing;
	double interchainspacing;
	double terminusspacing;
	double offset;
	double boxheight;
	int phenylcode;
    
	double offsetfraction;
	double monomerscaleoffsetfraction;
	double length;
} monolayer_alt;

typedef struct polymer_tag{
    double_triple **backboneoffset;      //  (first one for last polar type (0 if self is polar), second one for self)
    double **backbonenz;
    double **backbonephi;
    double_triple **sidechainrelativepos;
    double **sidechainnpar;
    double **sidechainnphi;
	
	double monomerspacing;
} polymer;

typedef struct bilayer_tag{
    double_triple ***backboneoffset;      //  (first one for leaf, second one for last polar type (0 if self is polar), third one for self)
    double ***backbonenz;
    double ***backbonephi;
    double_triple ***sidechainrelativepos;
    double ***sidechainnpar;
    double ***sidechainnphi;
    
	double monomerspacing;
	double interchainspacing;
	double terminusspacing;
	double offset;                       // these same for both leaves
	double boxheight;
	int phenylcode;
        
	double offsetfraction;
	double monomerscaleoffsetfraction;
	double length;
} bilayer;

typedef struct bilayer_alt_tag{
    double_triple ****backboneoffset;      //  (first one for leaf, second one for alternating index, third one for last polar type (0 if self is polar), fourth one for self)
    double ****backbonenz;
    double ****backbonephi;
    double_triple ****sidechainrelativepos;
    double ****sidechainnpar;
    double ****sidechainnphi;
    
	double monomerspacing;
	double interchainspacing;
	double terminusspacing;
	double offset;                       // these same for both leaves
	double boxheight;
	int phenylcode;
	
	double offsetfraction;
	double monomerscaleoffsetfraction;
	double length;
} bilayer_alt;

typedef struct linkedlist_tag{
	int ***head;
	int *list;
	int_triple *cell;
	int_triple cellsperside;
	int *number_neighbors;
	int **neighbor;
	double **pair_energy;
} linkedlist;

typedef struct linkedlistfull_tag{
	linkedlist core;
	int *reverselist;
	double_triple cellwidth;
	double mincellwidth;
	int maxneighbors;
} linkedlistfull;

typedef struct linkedlistasymmetric_tag{
	linkedlist core;
	int *reverselist;
    int *reverselist_reverse;
	double_triple cellwidth;
	double mincellwidth;
	int maxneighbors;
	int ***headreverse;
	int *listreverse;
	int_triple *cellreverse;
	int *number_neighbors_reverse;
	int **neighbor_reverse;
	int maxreverseneighbors;
} linkedlistasymmetric;

typedef struct linkedlistset_tag{
	linkedlistfull shortrange;
	linkedlistfull phenylphenyl;
	linkedlistfull chargedcharged;
} linkedlistset;

typedef struct linkdoubleparam_tag{
	double shortrange;
	double phenylphenyl;
	double chargedcharged;
} linkdoubleparam;

typedef struct pressureparam_tag{
    int pressuretype;
    double surfacepressure;
    double ypressure;
    double normalforceperunitarea;
} pressureparam;

typedef struct backbonehistoparams_tag{
	double dotwidth;
	double rwidth;
	double rbottom;
	int dotbins;
	int rbins;
} backbonehistoparams;

typedef struct sidechainhistoparams_tag{
    double rparallelmin;
    double rparallelwidth;
    double rperpwidth;
    double rdotnmin;
    double rdotnwidth;
    double ndotnwidth;
    int bins;
} sidechainhistoparams;

typedef struct reptationparams_tag{
	double rside0;
	double rforward0;
	double range;
	double range2;
	double maxconsecutivedirectordotproduct;
} reptationparams;

typedef struct sidereptationparams_tag{
	double r0;
	double range;
	double range2;
	double mindotproduct;
} sidereptationparams;

typedef struct displacementparams_tag{
    double av;
    long long int count;
	double lastsep;
} displacementparams;

typedef struct allatomparams_tag{
    double Cbetaspacing;
    double NterminusHdepth;
    double NterminusHwidth;
    double Cnyldepth;
    double Cnylwidth;
    double Onyldepth;
    double Onylwidth;
    double CterminusHdepth;
    double CterminusHwidth;
    double Cdepth;
    double Cwidth;
    double CHdepth;
    double CHwidth;
    double CHoutofplane;
	double phenylspacing;
    double phenylHspacing;
	double phenylCgammaspacing;
	double aminoNdepth;
	double aminoCdepth;
	double aminoHdepth;
	double aminoHwidth;
	double carboxylbackCdepth;
	double carboxylforwardCdepth;
	double carboxylOdepth;
	double carboxylOwidth;
    double ethylHdepth;
    double ethylHoutofplane;
    double CterminuslongHdistance;
    double CterminuslongHdepth;
    double CterminuslongHwidth;
	double Crad;
	double Orad;
	double Nrad;
	double Hrad;
} allatomparams;

typedef struct triangledata_tag{
	double area;
	double columnvolume;
	double_triple normal;
} triangledata;

typedef struct bonddata_tag{
    double length;
    double meancurvature;
} bonddata;

typedef struct vertexdata_tag{
    double meancurvature;
    double meancurvaturesquared;
    double area;
} vertexdata;

void parse_chemistry(char **chemistry_char, int *pchaintypes, int **Nchainsoftype, int **chainlengthoftype, int ***chainchemistry, int *pNmonomers, int *pNchains, monomertypes *pmonomercount, int sidechaintypesmax, int *pNnodetypes, int terminuscode);
void associate_chains_and_nodes(int Nmonomers, int Nchains, int chaintypes, int *chainsoftype, int *chainlengthoftype, int **chainchemistry, int *pNnodes, coord **coordarray, monomernodes ***monomerid, int **chainlength, bonded_params my_bonded_params);
void initialize_polymer_config(polymer *pconfig, double backbonezzigzag, double backbonenz, double backbonephi, int Nnodetypes, bonded_params my_bonded_params, double phenylheight);
void initialize_2dlattice_fromconfig(int Nmonomers, int Nchains, double density, double_triple *pbox_dimension, int *chainlength, monomernodes **monomerid, coord *coordarray, bonded_params my_bonded_params, nonbonded_params my_nonbonded_params, polymer config);
void input_and_copy_single_polymer(char *filename, int *pNnodes, int Nchains, int **chainlength, monomernodes ***monomerid, coord **coordarray, double_triple *pbox_dimension, int *pinterface, double *pvacuumthickness, int *pNnodetypes, monomertypes *pmonomercount, long long int *pt, int *pframe, double molarity, double cutoff2, int *prunningtime);
void input_bilayer(char *configname, bilayer *pconfig, int chainlength, int Nnodetypes, double residueshift);
void initialize_bilayer_fragment(int Nchains, double_triple *pbox_dimension, int *chainlength, monomernodes **monomerid, coord *coordarray, bonded_params my_bonded_params, nonbonded_params my_nonbonded_params, bilayer config, int leafsymmetry, int xchains, double molarity);
void create_ordered_bilayer(bilayer *pconfig, int chainlength, int Nnodetypes, double residueshift, double chargedbackbonezoffset, double backbonezzigzag, double backbonenz, double backbonephi, bonded_params my_bonded_params, double leafspacing, double leafxoffsetfrac, double leafyoffsetfrac, int leafsymmetry);
void initialize_periodic_bilayer_xchains(int Nchains, double_triple *pbox_dimension, int *chainlength, monomernodes **monomerid, coord *coordarray, bonded_params my_bonded_params, nonbonded_params my_nonbonded_params, bilayer config, int leafsymmetry, int xchains, int randomize_xoffsets);
void input_monolayer(char *configname, monolayer *pconfig, int chainlength, int Nnodetypes, double residueshift);
double_triple rsol(sidechain_params sidechain);
double nparhardsol(int orientationtype, double_triple *prelativepos, sidechain_params sidechain);
double nphihardsol(int orientationtype, double npar, double_triple relativepos, sidechain_params sidechain);
void initialize_periodic_monolayer_xchains(int Nchains, double_triple *pbox_dimension, int *chainlength, monomernodes **monomerid, coord *coordarray, bonded_params my_bonded_params, nonbonded_params my_nonbonded_params, monolayer config, int xchains, int randomize_xoffsets);
void initialize_periodic_monolayer_xchains_bothinterfaces(int Nchains, double_triple *pbox_dimension, int *chainlength, monomernodes **monomerid, coord *coordarray, bonded_params my_bonded_params, nonbonded_params my_nonbonded_params, monolayer config, int xchains, int randomize_xoffsets);
void initialize_2dxylattice_fromconfig(int Nmonomers, int Nchains, double density, double_triple *pbox_dimension, int *chainlength, monomernodes **monomerid, coord *coordarray, bonded_params my_bonded_params, nonbonded_params my_nonbonded_params, polymer config);
void output_coords(char *filename, int Nnodes, int Nchains, int *chainlength, monomernodes **monomerid, coord *coordarray, double_triple box_dimension, double vacuumthickness, int Nnodetypes, monomertypes monomercount, long long int t, int frame, int runningtime, double_triple *meshpositions, triangledata *meshtriangledata, bonddata *meshbonddata, vertexdata *meshdata, solvation_parameters solvationparameters);
void input_coords(char *filename, int *pNnodes, int *pNchains, int **chainlength, monomernodes ***monomerid, coord **coordarray, double_triple *pbox_dimension, double *pvacuumthickness, int *pNnodetypes, monomertypes *pmonomercount, long long int *pt, int *pframe, int *prunningtime, double_triple **pmeshpositions, triangledata **pmeshtriangledata, bonddata **pmeshbonddata, vertexdata **pmeshdata, solvation_parameters *psolvp);
void configure_cells_struct(linkedlistfull *plink, double_triple box_dimension);
void allocate_linklist(linkedlistfull *plink, int Nnodes);
void allocate_linklist_pairenergies(linkedlistfull *plink, int Nnodes);
void cellconstructdoublylinked(coord *coordarray, int N, linkedlistfull *plink, int lownodetype, int highnodetype);
void constructneighborlist(linkedlistfull *plink, int N, coord *coordarray, double_triple box_dimension);
void constructneighborlist_pairenergies_phenyl(linkedlistfull *plink, int N, coord *coordarray, double_triple box_dimension, nonbonded_params my_nonbonded_params);
void constructneighborlist_pairenergies_charged(linkedlistfull *plink, int N, coord *coordarray, double_triple box_dimension, nonbonded_params my_nonbonded_params);
void update_neighbor_list(int Nchains, int *chainlength, coord *coordarray, nonbonded_params my_nonbonded_params, linkedlistset linkset, monomernodes **monomerid, int *number_neighbors, int **neighbor, int max_neighbors, int *assigned);
void initialize_energycomponents_nofiles(energycomponents *pmy_energycomponents, int Nnodetypes);
void deletefile(char *filename);
void output_trajectory(int Nchains, int *chainlength, coord *coordarray, monomernodes **monomerid, double_triple box_dimension, char *filename, int *pframe, int cycle);
void mc_translate_single_hard(int mover, coord *coordarray, double max_translate, double_triple box_dimension, linkedlistfull *plink, bonded_params my_bonded_params, nonbonded_params my_nonbonded_params, double temperature, monomernodes **monomerid, int *chainlength, linkedlistasymmetric *pmeshlink, double_triple *meshpositions, int **mapmeshtoneighbors);
void mc_translate_single_hard_pull(int mover, coord *coordarray, double max_translate, double_triple box_dimension, linkedlistfull *plink, bonded_params my_bonded_params, nonbonded_params my_nonbonded_params, double temperature, monomernodes **monomerid, int *chainlength, int pulledchain, double force, linkedlistasymmetric *pmeshlink, double_triple *meshpositions, int **mapmeshtoneighbors);
void mc_translate_single_phenyl(int mover, coord *coordarray, double max_translate, double_triple box_dimension, linkedlistfull *plink1, linkedlistfull *plink2, bonded_params my_bonded_params, nonbonded_params my_nonbonded_params, double temperature, monomernodes **monomerid, int *chainlength, linkedlistasymmetric *pmeshlink, double_triple *meshpositions, int **mapmeshtoneighbors, int dovmmc, int *number_polymer_neighbors, int **polymer_neighbor, int *assigned, int Nchains, int max_neighbors);
void mc_translate_single_charged(int mover, coord *coordarray, double max_translate, double_triple box_dimension, linkedlistfull *plink1, linkedlistfull *plink2, bonded_params my_bonded_params, nonbonded_params my_nonbonded_params, double temperature, monomernodes **monomerid, int *chainlength, linkedlistasymmetric *pmeshlink, double_triple *meshpositions, int **mapmeshtoneighbors, int dovmmc, int *number_polymer_neighbors, int **polymer_neighbor, int *assigned, int Nchains, int max_neighbors);
void mc_rotate_single_hard(int mover, coord *coordarray, double max_rotate, double_triple box_dimension, linkedlist link, bonded_params my_bonded_params, nonbonded_params my_nonbonded_params, double temperature, monomernodes **monomerid, int *chainlength);
void mc_rotate_single_phenyl(int mover, coord *coordarray, double max_rotate, double_triple box_dimension, linkedlist link, bonded_params my_bonded_params, nonbonded_params my_nonbonded_params, double temperature, monomernodes **monomerid, int *chainlength);
void mc_rotate_single_charged(int mover, coord *coordarray, double max_rotate, double_triple box_dimension, linkedlist link, bonded_params my_bonded_params, nonbonded_params my_nonbonded_params, double temperature, monomernodes **monomerid, int *chainlength);
int mc_aspect_ratios_cellstruct(int Nnodes, int Nchains, coord *coordarray, coord *newcoordarray, double_triple *pbox_dimension, double maxlogaspect, bonded_params my_bonded_params, nonbonded_params *pmy_nonbonded_params, int *chainlength, monomernodes **monomerid, linkedlistset linkset, pressureparam my_pressureparam, double temperature, linkedlistasymmetric meshlink, double_triple *meshpositions, triangledata *meshtriangledata, bonddata *meshbonddata, vertexdata *meshdata, int **maptrianglestomesh, int **mapbondstomesh, int **mapmeshtoneighbors, int **mapmeshtobonds, int **mapmeshtotriangles, int **mapbondstotriangles);
int change_cells_doublylinked(linkedlistfull *plink, double_triple box_dimension);
void vmmc_translate_wholepolymer(int Nchains, int *chainlength, coord *newcoordarray, coord *reversecoordarray, coord *coordarray, double max_translate, double_triple box_dimension, linkedlistset *plinkset, bonded_params my_bonded_params, nonbonded_params my_nonbonded_params, double temperature, monomernodes **monomerid, int *number_neighbors, int **neighbor, int max_frustrated_links, int max_neighbors, int *interior_member, int *exterior_member, int *in_cluster, int_double *frustrated_link, int *assigned, int move3d, linkedlistasymmetric *pmeshlink, double_triple *meshpositions, int **mapmeshtoneighbors);
void mc_translate_wholepolymer(int Nchains, int *chainlength, coord *newcoordarray, coord *coordarray, double max_translate, double_triple box_dimension, linkedlistset *plinkset, bonded_params my_bonded_params, nonbonded_params my_nonbonded_params, double temperature, monomernodes **monomerid, linkedlistasymmetric *pmeshlink, double_triple *meshpositions, int **mapmeshtoneighbors);
void calc_cellstruct_leafid(coord *coordarray, bonded_params my_bonded_params, nonbonded_params my_nonbonded_params, int Nchains, int *chainlength, monomernodes **monomerid, double_triple box_dimension, energycomponents *pmy_energycomponents, linkedlistset linkset, double_triple *pavboxdimension);
void calc_nonbonded_energy_components_cellstruct_leafid(int n, int neighbor1, int neighbor2, int neighbor3, linkedlist link, coord *coordarray, double_triple box_dimension, nonbonded_params my_nonbonded_params, energycomponents *pmy_energycomponents);
void calc_nonbonded_energy_pair_components_leafid(coord a, coord b, double_triple box_dimension, nonbonded_params params, energycomponents *pmy_energycomponents);
void output_timeseries_old(char *totalfilename, energycomponents *pmy_energycomponents, int cycle, monomertypes monomercount, double factor, int Nnodetypes, double_triple *pavboxdimension, int runningtime, double *avmesharea, double *avmesharea2, double *avbondlength, double *avbondlength2, int Ntriangles, double *surfacetensionenergy, double *uniformbondenergy, double *pdvolume, double *pdvolumeenergy);
void output_timeseries(char *totalfilename, energycomponents *pmy_energycomponents, int cycle, monomertypes monomercount, double factor, int Nnodetypes, double_triple *pavboxdimension, int runningtime, double *avmesharea, double *avmesharea2, double *avmeancurvature, double *avmeancurvaturesquared, double *avbondlength, double *avbondlength2, int Ntriangles, double *surfacetensionenergy, double *pboundaryenergy, double *meancurvatureenergy, double *meancurvaturesquaredenergy, double *uniformbondenergy, double *pdvolume, double *pdvolumeenergy);
void output_timeseries_leaves(char *totalfilename, energycomponents *pmy_energycomponents, int cycle, monomertypes monomercount, double factor, int Nnodetypes, double_triple *pavboxdimension, int runningtime);
void calc_displacement(coord left, coord right, displacementparams *pdisplacement, double boxwidth);
void output_displacement_timeseries(char *filename, int cycle, displacementparams *pdisplacement);

double calc_sidechain_bonded_energy(sidechain_params sidechain, coord backbonecoord, coord sidechaincoord, double_triple box_dimension);
double calc_backbone_bonded_energy(bonded_params my_bonded_params, coord centercoord, coord leftcoord, coord rightcoord, double_triple box_dimension);
double calc_backbone_bonded_energy_onlyleft(bonded_params my_bonded_params, coord centercoord, coord leftcoord, double_triple box_dimension);
double calc_backbone_bonded_energy_norim1i(bonded_params my_bonded_params, coord centercoord, coord leftcoord, coord rightcoord, double_triple box_dimension);
double calc_backbone_bonded_difference_fourneighbors_changer(int twoleft, int left, int right, int tworight, bonded_params my_bonded_params, coord *coordarray, coord newcoord, coord oldcoord, double_triple box_dimension);
double calc_backbone_bonded_difference_twoneighbors_changen(int left, int right, bonded_params my_bonded_params, coord *coordarray, coord newcoord, coord oldcoord, double_triple box_dimension);

double calc_energy_difference_changen_hard(int n, coord coord_old, coord coord_new, coord *coordarray, double_triple box_dimension, bonded_params my_bonded_params, nonbonded_params my_nonbonded_params, monomernodes **monomerid, int *chainlength, linkedlist link);
double calc_energy_difference_changen_phenyl(int n, coord coord_old, coord coord_new, coord *coordarray, double_triple box_dimension, bonded_params my_bonded_params, nonbonded_params my_nonbonded_params, monomernodes **monomerid, int *chainlength, linkedlist link);
double calc_energy_difference_changen_charged(int n, coord coord_old, coord coord_new, coord *coordarray, double_triple box_dimension, bonded_params my_bonded_params, nonbonded_params my_nonbonded_params, monomernodes **monomerid, int *chainlength, linkedlist link);
double calc_energy_difference_changer_hard(int n, coord coord_old, coord coord_new, coord *coordarray, double_triple box_dimension, bonded_params my_bonded_params, nonbonded_params my_nonbonded_params, monomernodes **monomerid, int *chainlength, linkedlist link);
double calc_energy_difference_changer_phenyl(int n, coord coord_old, coord coord_new, coord *coordarray, double_triple box_dimension, bonded_params my_bonded_params, nonbonded_params my_nonbonded_params, monomernodes **monomerid, int *chainlength, linkedlist link1, linkedlist link2);
double calc_energy_difference_changer_charged(int n, coord coord_old, coord coord_new, coord *coordarray, double_triple box_dimension, bonded_params my_bonded_params, nonbonded_params my_nonbonded_params, monomernodes **monomerid, int *chainlength, linkedlist link1, linkedlist link2);

double calc_nonbonded_energy_pair(coord a, coord b, double_triple box_dimension, nonbonded_params params);
double calc_energy_pair(coord a, coord b, nonbonded_params my_nonbonded_params, double_triple box_dimension);
double hard_energy(double_triple ar, double_triple br, double arhard, double brhard, double_triple box_dimension);
double hard_energy_sumradii(double_triple ar, double_triple br, double sum, double_triple box_dimension);
double calc_hard_energy_cellstruct(int neighbor1, int neighbor2, int neighbor3, int numberneighbors, int *neighborlist, coord coord_new, coord *coordarray, double_triple box_dimension, nonbonded_params my_nonbonded_params);
double calc_hard_energy_cellstruct_otherpolymer(int numberneighbors, int *neighborlist, coord coord_new, coord *coordarray, double_triple box_dimension, nonbonded_params my_nonbonded_params);
double electrostatic_energy(coord a, coord b, solvation_parameters solvationparams, electrostaticparam pa, electrostaticparam pb, double_triple box_dimension, double cutoff2);

double calc_onebody_difference(coord coord_old, coord coord_new, onebodyparam p, solvation_parameters solvp);
void calc_onebody_energy_components(coord mycoord, onebodyparam p, solvation_parameters solvp, energycomponents *pmy_energycomponents);
onebodyparam chooseonebody(int type, nonbonded_params p);
double onebody_energy(coord mycoord, onebodyparam p, solvation_parameters solvp);

void outputcoords_centered(FILE *outp, char *name, double_triple pos, double_triple box_dimension);
double_triple nearest_image(double_triple a, double_triple b, double_triple box);
double_double nearest_image_3dto2d(double_triple a, double_triple b, double_triple box);

void updatecell(double_triple pos, linkedlistfull *plink, int n);
double total_energy_cellstruct(coord *coordarray, bonded_params my_bonded_params, nonbonded_params my_nonbonded_params, int Nchains, int *chainlength, monomernodes **monomerid, double_triple box_dimension, linkedlistset linkset);

double phenylphenyl_energy(coord a, coord b, GBQvac pvac, double_triple box_dimension, double cutoff2, solvation_parameters solvp);
double calc_phenyl_energy_cellstruct(int neighbor1, int neighbor2, int neighbor3, int numberneighbors, int *neighborlist, coord coord_new, coord *coordarray, double_triple box_dimension, nonbonded_params my_nonbonded_params);
double calc_charged_energy_cellstruct(int neighbor1, int neighbor2, int neighbor3, int numberneighbors, int *neighborlist, coord coord_new, coord *coordarray, double_triple box_dimension, nonbonded_params my_nonbonded_params);
double calc_phenyl_energy_cellstruct_otherpolymer(int numberneighbors, int *neighborlist, coord coord_new, coord *coordarray, double_triple box_dimension, nonbonded_params my_nonbonded_params);
double calc_charged_energy_cellstruct_otherpolymer(int numberneighbors, int *neighborlist, coord coord_new, coord *coordarray, double_triple box_dimension, nonbonded_params my_nonbonded_params);

void updateneighborlist(linkedlistfull *plink, int self, coord *coordarray, double_triple box_dimension);
void updateneighborlist_pairenergies_phenyl(linkedlistfull *plink, int self, coord *coordarray, double_triple box_dimension, nonbonded_params my_nonbonded_params);
void updateneighborlist_pairenergies_charged(linkedlistfull *plink, int self, coord *coordarray, double_triple box_dimension, nonbonded_params my_nonbonded_params);

int count_allatom_atoms(int Nchains, int *chainlength, monomernodes **monomerid, coord *coordarray);
void output_xyz_allatom(int Nchains, int *chainlength, monomernodes **monomerid, int Natoms, coord *coordarray, double_triple box_dimension, char *xyzname, char *sourcename, allatomparams params, int *pframe);

void output_timeseries_sheets(char *totalfilename, energycomponents *pmy_energycomponents, int cycle, monomertypes monomercount, double factor, int Nnodetypes, double_triple *pavboxdimension, int runningtime);

void input_coords_and_replicate_stacks(char *filename, int *pNnodes, int *pNchains, int **chainlength, monomernodes ***monomerid, coord **coordarray, double_triple *pbox_dimension, int *pinterface, double *pvacuumthickness, int *pNnodetypes, monomertypes *pmonomercount, long long int *pt, int *pframe, int *prunningtime, int replicates, double_triple displacement, double minsep);

double calc_hard_energy_cellstruct_otherbilayer(int numberneighbors, int *neighborlist, coord coord_new, coord *coordarray, double_triple box_dimension, nonbonded_params my_nonbonded_params);
double calc_phenyl_energy_cellstruct_otherbilayer(int numberneighbors, int *neighborlist, coord coord_new, coord *coordarray, double_triple box_dimension, nonbonded_params my_nonbonded_params);
double calc_charged_energy_cellstruct_otherbilayer(int numberneighbors, int *neighborlist, coord coord_new, coord *coordarray, double_triple box_dimension, nonbonded_params my_nonbonded_params);
void mc_translate_wholebilayer(int Nnodes, int Nsheets, int *chainlength, coord *newcoordarray, coord *coordarray, double max_translate, double_triple box_dimension, linkedlistset *plinkset, bonded_params my_bonded_params, nonbonded_params my_nonbonded_params, double temperature, monomernodes **monomerid, linkedlistasymmetric *pmeshlink, double_triple *meshpositions, int **mapmeshtoneighbors);

double calc_hard_energy_cellstruct_givenpolymer(int numberneighbors, int *neighborlist, coord coord_new, coord *coordarray, double_triple box_dimension, nonbonded_params my_nonbonded_params, int polymer);
double calc_phenyl_energy_cellstruct_givenpolymer(int numberneighbors, int *neighborlist, coord coord_new, coord *coordarray, double_triple box_dimension, nonbonded_params my_nonbonded_params, int polymer);
double calc_charged_energy_cellstruct_givenpolymer(int numberneighbors, int *neighborlist, coord coord_new, coord *coordarray, double_triple box_dimension, nonbonded_params my_nonbonded_params, int polymer);

int mc_shiftbilayergap(int Nnodes, int Nsheets, int *chainlength, coord *newcoordarray, coord *coordarray, double max_translate, double_triple *pbox_dimension, linkedlistset *plinkset, bonded_params my_bonded_params, nonbonded_params my_nonbonded_params, double temperature, monomernodes **monomerid, double normalforceperunitarea, int sitespersheet);

int count_cg_atoms(int Nchains, int *chainlength, monomernodes **monomerid, coord *coordarray);
void output_xyz_cg(int Nchains, int *chainlength, monomernodes **monomerid, int Natoms, coord *coordarray, double_triple box_dimension, char *xyzname, char *sourcename, cgparams params, int *pframe);

void expose_input_polymer_to_interface(char *filename, int *pNnodes, int *pNchains, int **chainlength, monomernodes ***monomerid, coord **coordarray, double_triple *pbox_dimension, int *pinterface, double *pvacuumthickness, int *pNnodetypes, monomertypes *pmonomercount, long long int *pt, int *pframe, int *prunningtime, double newheight, double newwidth, double newdepth, double newvacuumthickness);
void calc_com_trajectory(int Nnodes, coord *coordarray, double_triple box_dimension, double_triple *pavcom, int *pcomcount);
void output_com(char *filename, double_triple *pavcom, int *pcomcount, int cycle, int interface, double *avmeshheight, int meshnumber);

void add_hardenergy_if_polymer_not_in_list(coord newcoord, linkedlistfull mylinkedlistfull, coord *coordarray, int *list, double *energy, nonbonded_params my_nonbonded_params, double_triple box_dimension);
double_triple com(int firstnode, int lastnode, coord *coordarray, double_triple box_dimension);
int_triple cellofpos(double_triple pos, double_triple cellwidth);
void vmmc_rotate_wholepolymer(int Nchains, int *chainlength, coord *newcoordarray, coord *reversecoordarray, coord *coordarray, double max_rotate, double_triple box_dimension, linkedlistset *plinkset, bonded_params my_bonded_params, nonbonded_params my_nonbonded_params, double temperature, monomernodes **monomerid, int *number_neighbors, int **neighbor, int max_frustrated_links, int max_neighbors, int *interior_member, int *exterior_member, int *in_cluster, int_double *frustrated_link, int *assigned, int move3d, double *pmaxr, linkedlistasymmetric *pmeshlink, double_triple *meshpositions, int **mapmeshtoneighbors);
void add_phenylphenylenergy_if_polymer_not_in_list(coord newcoord, linkedlistfull mylinkedlistfull, coord *coordarray, int *list, double *energy, nonbonded_params my_nonbonded_params, double_triple box_dimension);
void add_chargedchargedenergy_if_polymer_not_in_list(coord newcoord, linkedlistfull mylinkedlistfull, coord *coordarray, int *list, double *energy, nonbonded_params my_nonbonded_params, double_triple box_dimension);
int calccomcheckr(double_triple *pcom, int firstnode, int lastnode, coord *coordarray, double_triple box_dimension, double *pmaxr);

double calc_hard_energy_cellstruct_givenpolymer_recreatecell(linkedlistfull mylinkedlistfull, coord coord_new, coord *coordarray, double_triple box_dimension, nonbonded_params my_nonbonded_params, int polymer);
double calc_phenyl_energy_cellstruct_givenpolymer_recreatecell(linkedlistfull mylinkedlistfull, coord coord_new, coord *coordarray, double_triple box_dimension, nonbonded_params my_nonbonded_params, int polymer);
double calc_charged_energy_cellstruct_givenpolymer_recreatecell(linkedlistfull mylinkedlistfull, coord coord_new, coord *coordarray, double_triple box_dimension, nonbonded_params my_nonbonded_params, int polymer);

double calc_hard_energy_cellstruct_specificbilayer(int numberneighbors, int *neighborlist, coord coord_new, coord *coordarray, double_triple box_dimension, nonbonded_params my_nonbonded_params, int othersheet);
double calc_phenyl_energy_cellstruct_specificbilayer(int numberneighbors, int *neighborlist, coord coord_new, coord *coordarray, double_triple box_dimension, nonbonded_params my_nonbonded_params, int othersheet);
double calc_charged_energy_cellstruct_specificbilayer(int numberneighbors, int *neighborlist, coord coord_new, coord *coordarray, double_triple box_dimension, nonbonded_params my_nonbonded_params, int othersheet);

void initialize_interface_mesh(solvation_parameters *psolvationparameters, double_triple box_dimension, double vacuumthickness, double_triple **pmeshpositions, triangledata **pmeshtriangledata, bonddata **pmeshbonddata, vertexdata **pmeshdata);
double height_relative_to_interface(int n, coord mycoord, solvation_parameters solvationparams, linkedlistasymmetric meshlink, double_triple box_dimension, double_triple *meshpositions, int **mapmeshtoneighbors);

double calc_onebody_energy(coord coord_new, onebodyparam p, solvation_parameters solvp);
double height_relative_to_interface_recreatecell(coord mycoord, solvation_parameters solvp, int n, linkedlistasymmetric meshlink, double_triple box_dimension, double_triple *meshpositions, int **mapmeshtoneighbors);
void output_xyz_mesh(double_triple *meshpositions, char *xyzname, char *sourcename, int *pframe, double siteradius, double_triple box_dimension, solvation_parameters solvp);
int mc_mesh(double_triple *meshpositions, solvation_parameters *psolvp, double max_translate, double_triple *pbox_dimension, double temperature, coord *coordarray, linkedlistasymmetric *pmeshlink, int Nnodes, nonbonded_params my_nonbonded_params, linkedlistfull *pphenlink, int *pmeshflag, triangledata *meshtriangledata, bonddata *meshbonddata, vertexdata *meshdata, linkedlistfull *pmeshmeshlink, int **mapmeshtotriangles, int **mapmeshtobonds, int **mapmeshtoneighbors, int **maptrianglestotriangles, int **maptrianglestomesh, int **mapmeshtoouterbonds, int **mapmeshtooutertriangles);
double triangular_area(double_triple p, double_triple q, double_triple r, double_triple box_dimension);
double projected_hexagon_area(double_triple c0, double_triple c1, double_triple c2, double_triple c3, double_triple c4, double_triple c5, double_triple box);
	
void allocate_linklist_asymmetric(linkedlistasymmetric *plink, int Nnodes, int reverseNodes);
void cellconstructasymmetric(double_triple *positions, int N, int reverseN, linkedlistasymmetric *plink, coord *coordarray);
void constructneighborlistasymmetric(linkedlistasymmetric *plink, int N, int Nreverse, coord *coordarray, double_triple box_dimension, double_triple *position);
void configure_cells_struct_asymmetric(linkedlistasymmetric *plink, double_triple box_dimension);
int change_cells_asymmetric(linkedlistasymmetric *plink, double_triple box_dimension);
void updateneighborlist_asymmetric(linkedlistasymmetric *plink, int self, double_triple selfpos, double_triple *positions, double_triple box_dimension);
void updatecell_asymmetric(double_triple pos, linkedlistasymmetric *plink, int n);
void updatecell_asymmetric_reverse(double_triple pos, linkedlistasymmetric *plink, int n);
void updateneighborlist_asymmetric_reverse(linkedlistasymmetric *plink, int self, double_triple selfpos, coord *coordarray, double_triple box_dimension);
double calc_onebody_difference_heights(double oldheight, double newheight, onebodyparam p, solvation_parameters solvp);
void calcavmesharea_bondlengths(solvation_parameters solvp, double *avarea, double *avarea2, double *avmeancurvature, double *avmeancurvaturesquared, double *avbondlength, double *avbondlength2, double *surfacetensionenergy, double *meancurvatureenergy, double *meancurvaturesquaredenergy, double *pboundaryenergy, double *uniformbondenergy, double *pdvolume, double *pdvolumeenergy, triangledata *meshtriangledata, bonddata *meshbonddata, vertexdata *meshdata, double_triple box_dimension, double surfacepressure);
double phenylphenyl_energy_heightchanged(coord a, coord b, GBQvac pvac, double_triple box_dimension, double cutoff2, solvation_parameters solvp, int heightchangeda, double heightnewa, int heightchangedb, double heightnewb);
triangledata calc_triangledata(double_triple p, double_triple q, double_triple r, double_triple box_dimension);

void height_verbose_debugging(int minindex, double_triple center, double_triple *meshpositions, int *neighbors, double_triple box_dimension, double_triple mycoordr, double_triple relativepos);

void cornersdontmatchdebug(int *edgeties, int *tieindex, int *nearoutsidecorner, int *tiecode);

void cellconstructdoublylinked_position(double_triple *position, int N, linkedlistfull *plink);
void constructneighborlist_positions(linkedlistfull *plink, int N, double_triple *positions, double_triple box_dimension);

int checkoverlap(double_triple *meshpositions, double_triple newpos, int meshnumber, int_double meshsize, int mover, int *neighbors, double_triple box_dimension, linkedlistfull meshmeshlink, triangledata *meshtriangledata, triangledata *newtriangledata, int **mapmeshtotriangles, int **maptrianglestotriangles, int **maptrianglestomesh);
int checktriangleoverlap(double_triple p1, double_triple q1, double_triple r1, double_triple p2, double_triple q2, double_triple r2, double_triple norm1, double_triple norm2, double_triple box_dimension);
int checkedgecrossinplane(double_double a0, double_double a1, double_double b0, double_double b1);
int checkpointintriangle(double_double point, double_double *vertices);
void updateneighborlist_positions(linkedlistfull *plink, int self, double_triple *positions, double_triple box_dimension);

int checktriangleoverlap_sharecorner(double_triple p1, double_triple q1, double_triple r1, double_triple p2, double_triple q2, double_triple r2, double_triple norm1, double_triple norm2, double_triple box_dimension);

void initialize_mesh_maps(int_double meshsize, int ***pmapmeshtotriangles, int ***pmapmeshtobonds, int ***pmaptrianglestomesh, int ***pmapbondstomesh, int ***pmapmeshtoneighbors, int ***pmaptrianglestotriangles, int ***pmapmeshtoouterbonds, int ***pmapmeshtooutertriangles, int ***pmapbondstotriangles);

double height_relative_to_interface_debug(int n, coord mycoord, solvation_parameters solvationparams, linkedlistasymmetric meshlink, double_triple box_dimension, double_triple *meshpositions, int **mapmeshtoneighbors);
int checkpointintrianglelist(double_double point, double_double v0, double_double v1, double_double v2);

void free_linklist(linkedlistfull link, int Nnodes);
void free_linklist_pairenergies(linkedlistfull link, int Nnodes);
void free_linklist_asymmetric(linkedlistasymmetric link, int Nnodes, int reversenodes);
void free_mesh_maps(int_double meshsize, int **mapmeshtotriangles, int **mapmeshtobonds, int **maptrianglestomesh, int **mapbondstomesh, int **mapmeshtoneighbors, int **maptrianglestotriangles, int **mapmeshtoouterbonds, int **mapmeshtooutertriangles, int **mapbondstotriangles);

void calc_cellstruct_leafid_fixedinterface(coord *coordarray, bonded_params my_bonded_params, nonbonded_params my_nonbonded_params, int Nchains, int *chainlength, monomernodes **monomerid, double_triple box_dimension, energycomponents *pmy_energycomponents, linkedlistset linkset, double *avheight, double_triple *pavboxdimension);
void output_timeseries_fixedinterface(char *totalfilename, energycomponents *pmy_energycomponents, int cycle, monomertypes monomercount, double factor, int Nnodetypes, double *avheight, double_triple *pavboxdimension, int runningtime);

void calc_cellstruct_leafid_cluster(coord *coordarray, bonded_params my_bonded_params, nonbonded_params my_nonbonded_params, int Nchains, int *chainlength, monomernodes **monomerid, double_triple box_dimension, energycomponents *pmy_energycomponents, linkedlistset linkset, double_triple *pavboxdimension, int *clusterdistribution);
void calc_nonbonded_energy_components_cellstruct_leafid_cluster(int n, int neighbor1, int neighbor2, int neighbor3, linkedlist link, coord *coordarray, double_triple box_dimension, nonbonded_params my_nonbonded_params, energycomponents *pmy_energycomponents, int *cluster_index_of_chain, int *size_of_cluster, int *head_of_cluster, int *list_of_chain, int *pNclusters);
void calc_nonbonded_energy_pair_components_leafid_cluster(coord a, coord b, double_triple box_dimension, nonbonded_params params, energycomponents *pmy_energycomponents, int *cluster_index_of_chain, int *size_of_cluster, int *head_of_cluster, int *list_of_chain, int *pNclusters);
void output_cluster(char *clusterfile, int *clusterdistribution, int cycle, double factor, int Nchains);

void calc_interface_avheight(int meshnumber, double_triple *meshpositions, double *avmeshheight);

void recalculate_meshbonddata_meshdata(int meshnumber, double_triple *meshpositions, bonddata *meshbonddata, vertexdata *meshdata, triangledata *meshtriangledata, double_triple box_dimension, int **mapbondstomesh, int **mapbondstotriangles, int **mapmeshtotriangles, int **mapmeshtobonds, double *totallength);

void output_seed(char *base);
void input_seed(char *inputfile);

void input_coords_hex(char *filename, int *pNnodes, int *pNchains, int **chainlength, monomernodes ***monomerid, coord **coordarray, double_triple *pbox_dimension, double *pvacuumthickness, int *pNnodetypes, monomertypes *pmonomercount, long long int *pt, int *pframe, int *prunningtime, double_triple **pmeshpositions, triangledata **pmeshtriangledata, bonddata **pmeshbonddata, vertexdata **pmeshdata, solvation_parameters *psolvp);
void output_coords_hex(char *filename, int Nnodes, int Nchains, int *chainlength, monomernodes **monomerid, coord *coordarray, double_triple box_dimension, double vacuumthickness, int Nnodetypes, monomertypes monomercount, long long int t, int frame, int runningtime, double_triple *meshpositions, triangledata *meshtriangledata, bonddata *meshbonddata, vertexdata *meshdata, solvation_parameters solvationparameters);
    
void output_coords_long(char *filename, int Nnodes, int Nchains, int *chainlength, monomernodes **monomerid, coord *coordarray, double_triple box_dimension, double vacuumthickness, int Nnodetypes, monomertypes monomercount, long long int t, int frame, int runningtime, double_triple *meshpositions, triangledata *meshtriangledata, bonddata *meshbonddata, vertexdata *meshdata, solvation_parameters solvationparameters);

double calc_phenyl_energy_cellstruct_polymerlist_singlecount(int numberneighbors, int *neighborlist, coord coord_new, coord *coordarray, double_triple box_dimension, nonbonded_params my_nonbonded_params, int *polymerlist);

