int input_trajectory(int Nchains, int chainlength, coord *coordarray, double_triple *pbox_dimension, FILE *inp, int *pcycle, int lim);
void poormansxrd_fromcoord_bin_leafid(int N, coord *coordarray, double_triple box_dimension, double ******Iqhisto, double xrdqspace, int xrdqn, int Nnodetypes, int ***xrdcount, int **edgeoncount, int *powdercount, int *radiallyaveragededgeoncount, double *****Iqedgeonhisto, double ****Iqpowderhisto, double ****Iqradiallyaveragededgeonhisto, int Nleaves);
void backbone_histo(int Nchains, int chainlength, coord *coordarray, double_triple box_dimension, backbonehistoparams params, int *histo_r, int *histo_rleftdotn, int *histo_rrightdotn, int *histo_tripleproduct, int **histo_r_rleftdotn, int **histo_r_rrightdotn, int **histo_tripleproduct_rleftdotn, int **histo_tripleproduct_rrightdotn, int *pcount_r, int *pcount_rleftdotn, int *pcount_rrightdotn, int *pcount_tripleproduct, int *pcount_r_rleftdotn, int *pcount_r_rrightdotn, int *pcount_tripleproduct_rleftdotn, int *pcount_tripleproduct_rrightdotn);
void gr_fromcoord(int N, coord *coordarray, double_triple box_dimension, double binwidth, double ****histo, int *pcount, int Nchains, int maxgrbin);
void output_gr_3dnorm(char *filename, double ****histo, int Nsites, double binwidth, int Nbins, int count, double *fractions, int type);
void output_gr_2dnorm(char *filename, double ****histo, int Nsites, double binwidth, int Nbins, int count, double *fractions, int type);
void output_dotproduct(char *filename, int bins, int *histo, int count, double binwidth);
void dotproducthistosametype_nomixed(int N, coord *coordarray, double_triple box_dimension, int Nchains, double threshold, double binwidth, int *histo, int *pcount, int type);
void dotproducthistosametype(int N, coord *coordarray, double_triple box_dimension, int Nchains, double threshold, double binwidth, int *histo, int *mixedhisto, int *pcount, int type);
void dotproducthistodifferenttypetwoshells(int N, coord *coordarray, double_triple box_dimension, int Nchains, double threshold, double outerthreshold, double binwidth, int *histo, int *pcount, int *outerhisto, int *poutercount, int typea, int typeb);
void output_backbone_histo(char *base, backbonehistoparams params, int *histo_r, int *histo_rleftdotn, int *histo_rrightdotn, int *histo_tripleproduct, int **histo_r_rleftdotn, int **histo_r_rrightdotn, int **histo_tripleproduct_rleftdotn, int **histo_tripleproduct_rrightdotn, int *pcount_r, int *pcount_rleftdotn, int *pcount_rrightdotn, int *pcount_tripleproduct, int *pcount_r_rleftdotn, int *pcount_r_rrightdotn, int *pcount_tripleproduct_rleftdotn, int *pcount_tripleproduct_rrightdotn);
void output_3d_xrd_leaves(char *base, double ******Iqhisto, double xrdqspace, int xrdqn, int Nnodetypes, int ***xrdcount, int Nnodes, int charlength, double *scatterfactor, int symmetries, double_triple box_dimension);
void output_2d_xrd_from3d_z0_leaves(char *base, double ******Iqhisto, double xrdqspace, int xrdqn, int Nnodetypes, int ***xrdcount, int Nnodes, int charlength, double *scatterfactor, int symmetries, double_triple box_dimension);
void output_2d_xrd_leaves(char *base, double *****Iqhisto, double xrdqspace, int xrdqn, int Nnodetypes, int **xrdcount, int Nnodes, int charlength, double *scatterfactor, int symmetries, double_triple box_dimension);
void output_1d_xrd_from2d_z0_leaves(char *base, double *****Iqhisto, double xrdqspace, int xrdqn, int Nnodetypes, int **xrdcount, int Nnodes, int charlength, double *scatterfactor, int symmetries, double_triple box_dimension);
void output_1d_xrd_leaves(char *base, double ****Iqhisto, double xrdqspace, int xrdqn, int Nnodetypes, int *xrdcount, int Nnodes, int charlength, double *scatterfactor, int symmetries, double_triple box_dimension);
int input_trajectory_v0(int Nchains, int chainlength, coord *coordarray, double_triple *pbox_dimension, FILE *inp, int *pcycle, int lim, int leafflag);
void assign_monomerid_and_chainlength(int Nchains, int uniformchainlength, monomernodes ***monomerid, int **chainlength);
int count_allatom_atoms_atomid(int Nchains, int *chainlength, monomernodes **monomerid, coord *coordarray, int *atomid, int *leaves);
void poormansxrd_fromcoord_bin_leafid_allatom(int Nchains, coord *coordarray, double_triple box_dimension, double ******Iqhisto, double xrdqspace, int xrdqn, int Nnodetypes, int ***xrdcount, int **edgeoncount, int *powdercount, int *radiallyaveragededgeoncount, double *****Iqedgeonhisto, double ****Iqpowderhisto, double ****Iqradiallyaveragededgeonhisto, int allatoms, int *atomid, int *leaves, double_triple *position, int *chainlength, monomernodes **monomerid, allatomparams params, int Nleaves);
int input_trajectory_blank(int Nchains, int chainlength, FILE *inp, int lim);
void poormansxrd_fromcoord_bin_leafid_allatom_no3d_resolutionlimited(int Nchains, coord *coordarray, double_triple box_dimension, double *****Iqhistoxy, double ****Iqhistoz, double_triple xrdqspace, int_triple xrdqn, int Nnodetypes, int *pxrdcount, int allatoms, int *atomid, int *leaves, double_triple *position, int *chainlength, monomernodes **monomerid, allatomparams params, int Nleaves, double_triple *pavbox);
void output_2d_xrd_from3d_z0_leaves_resolutionlimited(char *base, double ******Iqhisto, double_triple xrdqspace, int_triple xrdqn, int Nnodetypes, int ***xrdcount, int Nresidues, int charlength, double *scatterfactor, int symmetries, double_triple avbox, int boxcount);
void output_2d_xrd_leaves_resolutionlimited(char *base, double *****Iqhisto2d, double_triple xrdqspace, int_triple xrdqn, int Nnodetypes, int xrdcount, int Nresidues, int charlength, double *scatterfactor, int symmetries, double_triple avbox);
void output_1d_xrd_leaves_resolutionlimited(char *base, double ****Iqhisto1d, double_triple xrdqspace, int_triple xrdqn, int Nnodetypes, int xrdcount, int Nresidues, int charlength, double *scatterfactor, double *neutronscatterfactor, int symmetries, double_triple avbox);
void output_2d_xrd_leaves_resolutionlimited_andradav(char *base, double *****Iqhisto2d, double_triple xrdqspace, int_triple xrdqn, int Nnodetypes, int xrdcount, int Nresidues, int charlength, double *scatterfactor, double *neutronscatterfactor, int symmetries, double_triple avbox, double binspace);
void dotproducthistoalltypes(int N, coord *coordarray, double_triple box_dimension, int Nchains, double *thresholds, double binwidth, int **histo, int **mixedaddhisto, int **mixedsubtracthisto, int *count, int type, int issigned);
void output_unsigned_dotproduct_severaltypes(char *filename, int bins, int **histos, int *counts, double binwidth, int types, int countfactor);
void sidechain_histo(int Nchains, int chainlength, coord *coordarray, double_triple box_dimension, sidechainhistoparams params, int **histo_rparallel, int **histo_rperp, int **histo_rdotn, int **histo_ndotn, int *count);
void output_sidechain_histo(char *base, sidechainhistoparams params, int **histo_rparallel, int **histo_rperp, int **histo_rdotn, int **histo_ndotn, int *count);

