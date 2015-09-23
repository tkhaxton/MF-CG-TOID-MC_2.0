
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "mystdlib.h"
#include "mymath.h"
#include "peptoid_functions.h"

void parse_chemistry(char **chemistry_char, int *pchaintypes, int **Nchainsoftype, int **chainlengthoftype, int ***chainchemistry, int *pNmonomers, int *pNchains, monomertypes *pmonomercount, int sidechaintypesmax, int *pNnodetypes, int terminuscode){
	(*pNnodetypes)=0;
	(*pNmonomers)=0;
	(*pNchains)=0;
	(*pchaintypes)=atoi(chemistry_char[0]);
	(*Nchainsoftype)=xcalloc((*pchaintypes), sizeof(int));
	(*chainlengthoftype)=xcalloc((*pchaintypes), sizeof(int));
	int i, j, k, patterntype, counter=1, nonpolartype, firstpolartype, secondpolartype, Nblocks, blockcount, maxblocks=10, *polartype, *blocklength, count, extrachainlength;
	polartype=xcalloc(maxblocks, sizeof(int));
	blocklength=xcalloc(maxblocks, sizeof(int));
	(*chainchemistry)=xcalloc((*pchaintypes), sizeof(int *));
	(*pmonomercount).charged=(*pmonomercount).nonpolar=(*pmonomercount).monomers=0;
	for(i=0;i<sidechaintypesmax;i++) (*pmonomercount).type[i]=0;
	for(i=0;i<(*pchaintypes);i++){
		patterntype=atoi(chemistry_char[counter]);
		counter++;
		(*Nchainsoftype)[i]=atoi(chemistry_char[counter]);
		(*pNchains)+=(*Nchainsoftype)[i];
		counter++;
		if(patterntype==0){										//	alternating
			(*chainlengthoftype)[i]=atoi(chemistry_char[counter]);
			counter++;
			nonpolartype=atoi(chemistry_char[counter]);
			if(nonpolartype>=sidechaintypesmax) my_exit("type too big!");
			if(nonpolartype+1>(*pNnodetypes)) (*pNnodetypes)=nonpolartype+1;
			counter++;
			firstpolartype=atoi(chemistry_char[counter]);
			if(firstpolartype>=sidechaintypesmax) my_exit("type too big!");
			if(firstpolartype+1>(*pNnodetypes)) (*pNnodetypes)=firstpolartype+1;
			counter++;
			secondpolartype=atoi(chemistry_char[counter]);
			if(secondpolartype+1>(*pNnodetypes)) (*pNnodetypes)=secondpolartype+1;
			if(secondpolartype>=sidechaintypesmax) my_exit("type too big!");
			counter++;
            if(terminuscode==1) extrachainlength=1;
            else extrachainlength=0;
            (*chainlengthoftype)[i]+=extrachainlength;
			(*chainchemistry)[i]=xcalloc((*chainlengthoftype)[i], sizeof(int));
			for(j=0;j<(*chainlengthoftype)[i]-extrachainlength;j++){
				if(j%4==0){
					(*chainchemistry)[i][j]=firstpolartype;				//	start with polar
					(*pmonomercount).monomers+=(*Nchainsoftype)[i];
					(*pmonomercount).charged+=(*Nchainsoftype)[i];
					(*pmonomercount).nonpolar+=(*Nchainsoftype)[i];						//	backbone site
					(*pmonomercount).type[firstpolartype]+=(*Nchainsoftype)[i];
					(*pmonomercount).type[0]+=(*Nchainsoftype)[i];							//	backbone
				}
				if(j%2==1){
					(*chainchemistry)[i][j]=nonpolartype;
					(*pmonomercount).monomers+=(*Nchainsoftype)[i];
					(*pmonomercount).nonpolar+=2*(*Nchainsoftype)[i];						//	one for sidechain, one for backbone
					(*pmonomercount).type[nonpolartype]+=(*Nchainsoftype)[i];
					(*pmonomercount).type[0]+=(*Nchainsoftype)[i];							//	backbone
				}
				if(j%4==2){
					(*chainchemistry)[i][j]=secondpolartype;
					(*pmonomercount).monomers+=(*Nchainsoftype)[i];
					(*pmonomercount).charged+=(*Nchainsoftype)[i];
					(*pmonomercount).nonpolar+=(*Nchainsoftype)[i];						//	backbone site
					(*pmonomercount).type[secondpolartype]+=(*Nchainsoftype)[i];
					(*pmonomercount).type[0]+=(*Nchainsoftype)[i];							//	backbone
				}
			}
			(*pNmonomers)+=((*chainlengthoftype)[i])*((*Nchainsoftype)[i]);
            if(terminuscode==1){
                (*chainchemistry)[i][j]=-1;                              //  chainchemistry[i][j]=-1 denotes backbone site with no sidechain site
            }
		}
		if(patterntype==1){																//	block
			Nblocks=atoi(chemistry_char[counter]);
			if(Nblocks>maxblocks) my_exit("Nblocks too big!");
			counter++;
			nonpolartype=atoi(chemistry_char[counter]);
			if(nonpolartype>=sidechaintypesmax) my_exit("type too big!");
			if(nonpolartype+1>(*pNnodetypes)) (*pNnodetypes)=nonpolartype+1;
			counter++;
			blockcount=0;
            (*chainlengthoftype)[i]=0;
			while(blockcount<Nblocks){
				polartype[blockcount]=atoi(chemistry_char[counter]);
				if(polartype[blockcount]>=sidechaintypesmax) my_exit("type too big!");
				if(polartype[blockcount]+1>(*pNnodetypes)) (*pNnodetypes)=polartype[blockcount]+1;
				counter++;
				blocklength[blockcount]=atoi(chemistry_char[counter]);			//	 in dimers
				(*chainlengthoftype)[i]+=2*blocklength[blockcount];
				counter++;
				blockcount++;
			}
			count=0;
            
            if(terminuscode==1) extrachainlength=1;
            else extrachainlength=0;
            (*chainlengthoftype)[i]+=extrachainlength;
            
			(*chainchemistry)[i]=xcalloc((*chainlengthoftype)[i], sizeof(int));
			for(j=0;j<blockcount;j++){
				for(k=0;k<blocklength[j];k++){
					(*chainchemistry)[i][count]=polartype[j];				//	start with polar
					(*pmonomercount).monomers+=(*Nchainsoftype)[i];
					(*pmonomercount).charged+=(*Nchainsoftype)[i];
					(*pmonomercount).nonpolar+=(*Nchainsoftype)[i];							//	backbone site
					(*pmonomercount).type[polartype[j]]+=(*Nchainsoftype)[i];
					(*pmonomercount).type[0]+=(*Nchainsoftype)[i];								//	backbone
					count++;
					(*chainchemistry)[i][count]=nonpolartype;
					(*pmonomercount).monomers+=(*Nchainsoftype)[i];
					(*pmonomercount).nonpolar+=2*(*Nchainsoftype)[i];							//	one for sidechain, one for backbone
					(*pmonomercount).type[nonpolartype]+=(*Nchainsoftype)[i];
					(*pmonomercount).type[0]+=(*Nchainsoftype)[i];								//	backbone
					count++;
				}
			}
			(*pNmonomers)+=((*chainlengthoftype)[i])*((*Nchainsoftype)[i]);
            if(terminuscode==1){
                (*chainchemistry)[i][count]=-1;                              //  chainchemistry[i][j]=-1 denotes backbone site with no sidechain site
            }
		}
		if(patterntype==2){																//	block, not necessarily alternating hydrophobic-hydrophilic
			Nblocks=atoi(chemistry_char[counter]);
			if(Nblocks>maxblocks) my_exit("Nblocks too big!");
			counter++;
			blockcount=0;
			(*chainlengthoftype)[i]=0;
			while(blockcount<Nblocks){
				polartype[blockcount]=atoi(chemistry_char[counter]);
				if(polartype[blockcount]>=sidechaintypesmax) my_exit("type too big!");
				if(polartype[blockcount]+1>(*pNnodetypes)) (*pNnodetypes)=polartype[blockcount]+1;
				counter++;
				blocklength[blockcount]=atoi(chemistry_char[counter]);			//	 in MONOMERS
				counter++;
				(*chainlengthoftype)[i]+=blocklength[blockcount];
				blockcount++;
			}
			count=0;
            
            if(terminuscode==1) extrachainlength=1;
            else extrachainlength=0;
            (*chainlengthoftype)[i]+=extrachainlength;
            
			(*chainchemistry)[i]=xcalloc((*chainlengthoftype)[i], sizeof(int));
			for(j=0;j<blockcount;j++){
				for(k=0;k<blocklength[j];k++){
					(*chainchemistry)[i][count]=polartype[j];				//	start with polar
					(*pmonomercount).monomers+=(*Nchainsoftype)[i];
					(*pmonomercount).charged+=(*Nchainsoftype)[i];
					(*pmonomercount).nonpolar+=(*Nchainsoftype)[i];							//	backbone site
					(*pmonomercount).type[polartype[j]]+=(*Nchainsoftype)[i];
					(*pmonomercount).type[0]+=(*Nchainsoftype)[i];								//	backbone
					count++;
				}
			}
			(*pNmonomers)+=((*chainlengthoftype)[i])*((*Nchainsoftype)[i]);
            if(terminuscode==1){
                (*chainchemistry)[i][count]=-1;                              //  chainchemistry[i][j]=-1 denotes backbone site with no sidechain site
            }
		}
	}
    free(polartype);
    free(blocklength);
}

void associate_chains_and_nodes(int Nmonomers, int Nchains, int chaintypes, int *chainsoftype, int *chainlengthoftype, int **chainchemistry, int *pNnodes, coord **coordarray, monomernodes ***monomerid, int **chainlength, bonded_params my_bonded_params){
	int i, j,k,  nodespermonomer=2, chaincounter=0, nodecounter=0;
	(*pNnodes)=Nmonomers*nodespermonomer;
	(*coordarray)=xcalloc((*pNnodes), sizeof(coord));
	(*monomerid)=xcalloc(Nchains, sizeof(monomernodes *));
	(*chainlength)=xcalloc(Nchains, sizeof(int));
	for(i=0;i<chaintypes;i++){
		for(j=0;j<chainsoftype[i];j++){
			(*monomerid)[chaincounter]=xcalloc(chainlengthoftype[i], sizeof(monomernodes));
			(*chainlength)[chaincounter]=chainlengthoftype[i];
			for(k=0;k<chainlengthoftype[i];k++){
				(*monomerid)[chaincounter][k].backbone=nodecounter;
				(*coordarray)[nodecounter].chainid=chaincounter;
				(*coordarray)[nodecounter].monomerid=k;
				(*coordarray)[nodecounter].nodetype=0;							//	backbone
                (*coordarray)[nodecounter].leafid=0;                            //  default
				nodecounter++;
				(*monomerid)[chaincounter][k].sidechain=nodecounter;
				(*coordarray)[nodecounter].chainid=chaincounter;
				(*coordarray)[nodecounter].monomerid=k;
				(*coordarray)[nodecounter].nodetype=chainchemistry[i][k];		//	sidechain
                if((*coordarray)[nodecounter].nodetype==-1) (*coordarray)[nodecounter].orientationtype=-1;
				else (*coordarray)[nodecounter].orientationtype=my_bonded_params.sidechainorientationtype[(*coordarray)[nodecounter].nodetype];
                (*coordarray)[nodecounter].leafid=0;                            //  default
				nodecounter++;
			}
			chaincounter++;
		}
	}
}

void initialize_polymer_config(polymer *pconfig, double backbonezzigzag, double backbonenz, double backbonephi, int Nnodetypes, bonded_params my_bonded_params, double phenylheight){
    int i, j, last, self;
    allocate_matrix_nozero(double_triple, (*pconfig).backboneoffset, Nnodetypes, Nnodetypes);
    allocate_matrix(double, (*pconfig).backbonenz, Nnodetypes, Nnodetypes);
    allocate_matrix(double, (*pconfig).backbonephi, Nnodetypes, Nnodetypes);
    allocate_matrix_nozero(double_triple, (*pconfig).sidechainrelativepos, Nnodetypes, Nnodetypes);
    allocate_matrix(double, (*pconfig).sidechainnpar, Nnodetypes, Nnodetypes);
    allocate_matrix(double, (*pconfig).sidechainnphi, Nnodetypes, Nnodetypes);
    (*pconfig).backbonenz[0][2]=backbonenz;
    (*pconfig).backbonenz[0][3]=backbonenz;
    (*pconfig).backbonenz[2][1]=backbonenz;
    (*pconfig).backbonenz[3][1]=backbonenz;                                      //  perpendicular, keeping phi at 0
    (*pconfig).backbonephi[0][2]=backbonephi;
    (*pconfig).backbonephi[0][3]=backbonephi;
    (*pconfig).backbonephi[2][1]=backbonephi;
    (*pconfig).backbonephi[3][1]=backbonephi;
    (*pconfig).sidechainrelativepos[2][1]=(*pconfig).sidechainrelativepos[3][1]=rsol(my_bonded_params.sidechain[1]);
    (*pconfig).sidechainrelativepos[0][2]=rsol(my_bonded_params.sidechain[2]);
    (*pconfig).sidechainrelativepos[0][3]=rsol(my_bonded_params.sidechain[3]);
    (*pconfig).sidechainnpar[2][1]=(*pconfig).sidechainnpar[3][1]=nparhardsol(my_bonded_params.sidechainorientationtype[1], &((*pconfig).sidechainrelativepos[2][1]), my_bonded_params.sidechain[1]);
    (*pconfig).sidechainnpar[0][2]=nparhardsol(my_bonded_params.sidechainorientationtype[2], &((*pconfig).sidechainrelativepos[0][2]), my_bonded_params.sidechain[2]);
    (*pconfig).sidechainnpar[0][3]=nparhardsol(my_bonded_params.sidechainorientationtype[3], &((*pconfig).sidechainrelativepos[0][3]), my_bonded_params.sidechain[3]);
    (*pconfig).sidechainnphi[2][1]=(*pconfig).sidechainnphi[3][1]=nphihardsol(my_bonded_params.sidechainorientationtype[1], (*pconfig).sidechainnpar[2][1], (*pconfig).sidechainrelativepos[2][1], my_bonded_params.sidechain[1]);
    (*pconfig).sidechainnphi[0][2]=nphihardsol(my_bonded_params.sidechainorientationtype[2], (*pconfig).sidechainnpar[0][2], (*pconfig).sidechainrelativepos[0][2], my_bonded_params.sidechain[2]);
    (*pconfig).sidechainnphi[0][3]=nphihardsol(my_bonded_params.sidechainorientationtype[3], (*pconfig).sidechainnpar[0][3], (*pconfig).sidechainrelativepos[0][3], my_bonded_params.sidechain[3]);
    for(i=0;i<2;i++){
        for(j=2;j<4;j++){
            if(i==0){
                last=0;
                self=j;
            }
            else{
                last=j;
                self=1;
            }
            (*pconfig).backboneoffset[last][self].z=phenylheight-(*pconfig).sidechainrelativepos[2][1].z;
            if(last==0) (*pconfig).backboneoffset[last][self].z-=backbonezzigzag;
            (*pconfig).backboneoffset[last][self].x=(*pconfig).backboneoffset[last][self].y=0;
        }
    }
    
}

void initialize_2dlattice_fromconfig(int Nmonomers, int Nchains, double density, double_triple *pbox_dimension, int *chainlength, monomernodes **monomerid, coord *coordarray, bonded_params my_bonded_params, nonbonded_params my_nonbonded_params, polymer config){
	int chainsonaside, i, j, k, lastpolar, lastpolartype, chaincounter=0, backboneid, sidechainid, selftype;
    double nz, phi, nx, ny, npar, nphi;
    double_triple position, nalong, nside, xvector, sidepos, nback;
    xvector.x=1;xvector.y=xvector.z=0;
	(*pbox_dimension).x=pow(Nmonomers/density, 1.0/3.0);
	(*pbox_dimension).y=(*pbox_dimension).z=(*pbox_dimension).x;
	chainsonaside=(int) floor(sqrt(1.*Nchains))+1;
	double chainspacing=(*pbox_dimension).x/chainsonaside;
	double chainspacingz, offsetz;
	if(my_nonbonded_params.solvationparams.interface==1){
		chainspacingz=((*pbox_dimension).z-my_nonbonded_params.vacuumthickness)/chainsonaside;
		offsetz=0.5*my_nonbonded_params.vacuumthickness;
	}
	else{
		chainspacingz=chainspacing;
		offsetz=0;
	}
	for(i=0;i<chainsonaside;i++){
		for(j=0;j<chainsonaside;j++){
			position.x=0.5*(*pbox_dimension).x;
			position.y=(i+0.5)*chainspacing;
			position.z=offsetz+(j+0.5)*chainspacingz;
			for(k=0;k<chainlength[chaincounter];k++){
				backboneid=monomerid[chaincounter][k].backbone;
				sidechainid=monomerid[chaincounter][k].sidechain;
                if(coordarray[sidechainid].nodetype<0){             //  extra backbone, assume alternating amphiphilic and ending with nonpolar
                    lastpolartype=0;
                    selftype=lastpolar;
                }
				else if(coordarray[sidechainid].nodetype>1){
					lastpolar=coordarray[sidechainid].nodetype;
					lastpolartype=0;                                //  self is polar
                    selftype=coordarray[sidechainid].nodetype;
				}
				else{
                    lastpolartype=lastpolar;
                    selftype=coordarray[sidechainid].nodetype;
                }
				coordarray[backboneid].r=position;
				coordarray[backboneid].r.x+=(-0.5*chainlength[chaincounter]+k)*config.monomerspacing;
				coordarray[backboneid].r=add_double_triple(coordarray[backboneid].r, config.backboneoffset[lastpolartype][selftype]);
				fmod_double_triple(&(coordarray[backboneid].r), (*pbox_dimension));
				nz=config.backbonenz[lastpolartype][selftype];
				coordarray[backboneid].n.z=((k%2)*2-1)*nz;          //  amphiphilic
				phi=config.backbonephi[lastpolartype][selftype];
				coordarray[backboneid].n.x=sqrt(1-nz*nz)*cos(phi);
				coordarray[backboneid].n.y=sqrt(1-nz*nz)*sin(phi);
				
				if(coordarray[sidechainid].nodetype>=0){
					nback=coordarray[backboneid].n;
					sidepos=config.sidechainrelativepos[lastpolartype][coordarray[sidechainid].nodetype];
					nalong=subtract_double_triple(xvector, scalar_multiply_double_triple(nback, nback.x));  //  projection onto r_{back-side} but perp to nhat, assuming r_{back-side} has no y component
					normalize(&nalong);
					nside=cross_product(nback, nalong);
					npar=config.sidechainnpar[lastpolartype][coordarray[sidechainid].nodetype];
					nphi=config.sidechainnphi[lastpolartype][coordarray[sidechainid].nodetype];
					coordarray[sidechainid].r=add_double_triple(add_double_triple(add_double_triple(coordarray[backboneid].r, scalar_multiply_double_triple(nback, sidepos.z)), scalar_multiply_double_triple(nalong, sidepos.x)), scalar_multiply_double_triple(nside, sidepos.y));
					fmod_double_triple(&(coordarray[sidechainid].r), (*pbox_dimension));
					nx=sqrt(1-pow(npar, 2))*cos(nphi);
					ny=sqrt(1-pow(npar, 2))*sin(nphi);
					coordarray[sidechainid].n=add_double_triple(add_double_triple(scalar_multiply_double_triple(nback, npar), scalar_multiply_double_triple(nalong, nx)), scalar_multiply_double_triple(nside, ny));
				}
			}
			chaincounter++;
			if(chaincounter==Nchains) return;
		}
	}
}

void input_and_copy_single_polymer(char *filename, int *pNnodes, int Nchains, int **chainlength, monomernodes ***monomerid, coord **coordarray, double_triple *pbox_dimension, int *pinterface, double *pvacuumthickness, int *pNnodetypes, monomertypes *pmonomercount, long long int *pt, int *pframe, double molarity, double cutoff2, int *prunningtime){
	int i, j, inputchains, inputnodes, success, k, l;
	FILE *inp;
	inp=fopen(filename, "r");
	double_triple input_box_dimension, com, sep, try;
    sep.x=sep.y=sep.z=0;
	double maxr=0, packingfraction;
	fscanf(inp, "%lf %lf %lf", &((input_box_dimension).x), &((input_box_dimension).y), &((input_box_dimension).z));
    fscanf(inp, "%i", &(*pinterface));
    if((*pinterface)==1){
        fscanf(inp, " %lf", pvacuumthickness);
    }   
	fscanf(inp, "%i %i %i %i %i %i", &inputnodes, &inputchains, &((*pmonomercount).charged), &((*pmonomercount).nonpolar), &((*pmonomercount).monomers), &(*pNnodetypes));
	
	((*pmonomercount).charged)*=Nchains;
	((*pmonomercount).nonpolar)*=Nchains;
	((*pmonomercount).monomers)*=Nchains;
	
	if(inputchains!=1){
		printf("inputchains=%i! (%f %f %f %i)\n", inputchains, input_box_dimension.x, input_box_dimension.y, input_box_dimension.z, inputnodes);
		exit(1);
	}
	(*pNnodes)=Nchains*inputnodes;
	
	for(i=0;i<(*pNnodetypes);i++){
		fscanf(inp, " %i", &((*pmonomercount).type[i]));
		((*pmonomercount).type[i])*=Nchains;
	}
	(*coordarray)=xcalloc((*pNnodes), sizeof(coord));
	(*chainlength)=xcalloc((Nchains), sizeof(int));
	(*monomerid)=xcalloc((Nchains), sizeof(monomernodes *));
	fscanf(inp, "%i", &((*chainlength)[0]));								//	scan in one polymer, assign its attributes to all
	((*monomerid)[0])=xcalloc(((*chainlength)[0]), sizeof(monomernodes));
	for(j=0;j<(*chainlength)[0];j++){
		fscanf(inp, "%i %i", &((*monomerid)[0][j].backbone), &((*monomerid)[0][j].sidechain));
	}
	for(i=1;i<Nchains;i++){
		(*chainlength)[i]=(*chainlength)[0];
		((*monomerid)[i])=xcalloc(((*chainlength)[i]), sizeof(monomernodes));
		for(j=0;j<(*chainlength)[0];j++){
			(*monomerid)[i][j].backbone=(*monomerid)[0][j].backbone+i*inputnodes;
			(*monomerid)[i][j].sidechain=(*monomerid)[0][j].sidechain+i*inputnodes;
		}
	}
    
    if((*pinterface)==1){
        (*pbox_dimension).x=(*pbox_dimension).y=pow(molarity*Nchains*chainlength[0][0], 1.0/2.0);
        (*pbox_dimension).z=(input_box_dimension).z;
    }
    else{
        (*pbox_dimension).x=pow((1.*Nchains)/(molarity*molardensity), 1.0/3.0);
        (*pbox_dimension).y=(*pbox_dimension).z=(*pbox_dimension).x;
    }    
    
	com.x=com.y=com.z=0;
    int realnodes=0;
	for(i=0;i<inputnodes;i++){
		fscanf(inp, "%lf %lf %lf %lf %lf %lf %i %i %i %i %i", &((*coordarray)[i].r.x), &((*coordarray)[i].r.y), &((*coordarray)[i].r.z), &((*coordarray)[i].n.x), &((*coordarray)[i].n.y), &((*coordarray)[i].n.z), &((*coordarray)[i].nodetype), &((*coordarray)[i].chainid), &((*coordarray)[i].monomerid), &((*coordarray)[i].orientationtype), &((*coordarray)[i].leafid));
        if((*coordarray)[i].nodetype>=0){
            realnodes++;
            if(i==0){
                com=(*coordarray)[i].r;
            }
            else{
                sep=subtract_double_triple((*coordarray)[i].r, (*coordarray)[i-1].r);
                recenter_double_triple(&sep, input_box_dimension);
                (*coordarray)[i].r=add_double_triple((*coordarray)[i-1].r, sep);
                com=add_double_triple(com, (*coordarray)[i].r);
            }
        }
	}
	com=scalar_multiply_double_triple(com, (1./realnodes));
    if((*pinterface)==1) com.z=0;
    double norm2;
	for(i=0;i<inputnodes;i++){
        if((*coordarray)[i].nodetype>=0){
            (*coordarray)[i].r=subtract_double_triple((*coordarray)[i].r, com);
            recenter_double_triple(&((*coordarray)[i].r), input_box_dimension);
            norm2=(*coordarray)[i].r.x*(*coordarray)[i].r.x+(*coordarray)[i].r.y*(*coordarray)[i].r.y;
            if((*pinterface)==0) norm2+=(*coordarray)[i].r.z*(*coordarray)[i].r.z;
            if(norm2>maxr*maxr) maxr=sqrt(norm2);
        }
	}
	maxr=maxr+0.5*sqrt(cutoff2);
    if((*pinterface)==1) packingfraction=Nchains*M_PI*pow(maxr, 2)/((*pbox_dimension).x*(*pbox_dimension).y);
	else packingfraction=4.*Nchains*M_PI/3.*pow(maxr, 3)/((*pbox_dimension).x*(*pbox_dimension).y*(*pbox_dimension).z);
	if(packingfraction>0.5){
        printf("packing fraction=%f!\n", packingfraction);
        exit(1);
    }
	double_triple *centers;
	centers=xcalloc(Nchains, sizeof(double_triple));
	centers[0].x=centers[0].y=centers[0].z=0;
    for(i=1;i<Nchains;i++){
		success=0;
		while(success==0){
            if((*pinterface)==1){
                try.x=(*pbox_dimension).x*rand_double;
                try.y=(*pbox_dimension).y*rand_double;
                try.z=0;
            }
			else{
                try=rand_unit_cube();
                try=scalar_multiply_double_triple(try, (*pbox_dimension).x);
			}
            success=1;
			for(j=0;j<i;j++){
                sep=subtract_double_triple(centers[j], try);
				recenter_double_triple(&sep, (*pbox_dimension));
				if(norm(sep)<2*maxr) success=0;
			}
		}
		centers[i]=try;
	}
    double_triple axis_vector;
	declare_matrix(double, one_matrix, 3, 3);
	declare_matrix(double, new_matrix, 3, 3);
	declare_matrix(double, total_matrix, 3, 3);
    double angle;
    for(i=0;i<Nchains;i++){
        if((*pinterface)==0){
            for(j=0;j<3;j++){
                axis_vector.x=axis_vector.y=axis_vector.z=0;
                if(j==0) axis_vector.x=1;
                else if(j==1) axis_vector.y=1;
                else axis_vector.z=1;
                angle=2*M_PI*rand_double;
                forward_matrix(angle, axis_vector, one_matrix);
                if(j==0){
                    for(k=0;k<3;k++){
                        for(l=0;l<3;l++){
                            total_matrix[k][l]=one_matrix[k][l];
                        }
                    }
                }
                else{
                    matrix_multiply(total_matrix, one_matrix, new_matrix);
                    for(k=0;k<3;k++){
                        for(l=0;l<3;l++){
                            total_matrix[k][l]=new_matrix[k][l];
                        }
                    }
                }
            }
        }
        else{
            angle=2*M_PI*rand_double;
            axis_vector.x=axis_vector.y=0;
            axis_vector.z=1;
            forward_matrix(angle, axis_vector, total_matrix);
        }
		for(j=0;j<inputnodes;j++){
			(*coordarray)[i*inputnodes+j]=(*coordarray)[j];
			(*coordarray)[i*inputnodes+j].r=add_double_triple(centers[i], rotate_by_matrix((*coordarray)[j].r, total_matrix));
			(*coordarray)[i*inputnodes+j].n=rotate_by_matrix((*coordarray)[j].n, total_matrix);
			(*coordarray)[i*inputnodes+j].chainid=i;
		}
	}
    for(i=0;i<(*pNnodes);i++){
        fmod_double_triple(&((*coordarray)[i].r), (*pbox_dimension));
    }
    free(centers);
    free_matrix(one_matrix, 3);
    free_matrix(new_matrix, 3);
    free_matrix(total_matrix, 3);
	fscanf(inp, "%lli %i %i", &(*pt), &(*pframe), &(*prunningtime));
}

void input_bilayer(char *configname, bilayer *pconfig, int chainlength, int Nnodetypes, double residueshift){
 	FILE *fp;
	fp=fopen(configname, "r");
	int size=1000, i, j, k;
	char *buff;
	buff=xcalloc(size, sizeof(char));
	while(fgets(buff, size, fp));
	fclose(fp);
	int offset;
	sscanf(buff, "%*i%n", &offset);
	buff+=offset;
	sscanf(buff, "%*f%n", &offset);
	buff+=offset;
	sscanf(buff, "%lf%n", &((*pconfig).monomerspacing), &offset);
	buff+=offset;
	sscanf(buff, "%lf%n", &((*pconfig).monomerscaleoffsetfraction), &offset);
	buff+=offset;
	sscanf(buff, "%lf%n", &((*pconfig).terminusspacing), &offset);
	buff+=offset;
	sscanf(buff, "%lf%n", &((*pconfig).interchainspacing), &offset);
	buff+=offset;
    
    (*pconfig).length=chainlength*(*pconfig).monomerspacing+(*pconfig).terminusspacing;
    (*pconfig).offset=(*pconfig).length*(*pconfig).offsetfraction+((*pconfig).monomerscaleoffsetfraction+residueshift)*(*pconfig).monomerspacing;
    
    allocate_3d_tensor_nozero(double_triple, (*pconfig).backboneoffset, 2, Nnodetypes, Nnodetypes);
    allocate_3d_tensor(double, (*pconfig).backbonenz, 2, Nnodetypes, Nnodetypes);
    allocate_3d_tensor(double, (*pconfig).backbonephi, 2, Nnodetypes, Nnodetypes);
    allocate_3d_tensor_nozero(double_triple, (*pconfig).sidechainrelativepos, 2, Nnodetypes, Nnodetypes);
    allocate_3d_tensor(double, (*pconfig).sidechainnpar, 2, Nnodetypes, Nnodetypes);
    allocate_3d_tensor(double, (*pconfig).sidechainnphi, 2, Nnodetypes, Nnodetypes);
    
    int leaf, charged, chargetype, last, self;
	for(leaf=0;leaf<2;leaf++){
        for(charged=0;charged<2;charged++){
            for(chargetype=2;chargetype<4;chargetype++){
                if(charged==1){
                    last=0;
                    self=chargetype;
                }
                else{
                    last=chargetype;
                    self=1;
                }
                sscanf(buff, "%lf%n", &((*pconfig).backboneoffset[leaf][last][self].x), &offset);
                buff+=offset;
                sscanf(buff, "%lf%n", &((*pconfig).backboneoffset[leaf][last][self].y), &offset);
                buff+=offset;
                sscanf(buff, "%lf%n", &((*pconfig).backboneoffset[leaf][last][self].z), &offset);
                buff+=offset;
                sscanf(buff, "%lf%n", &((*pconfig).backbonenz[leaf][last][self]), &offset);
                buff+=offset;
                sscanf(buff, "%lf%n", &((*pconfig).backbonephi[leaf][last][self]), &offset);
                buff+=offset;
                sscanf(buff, "%lf%n", &((*pconfig).sidechainrelativepos[leaf][last][self].x), &offset);
                buff+=offset;
                sscanf(buff, "%lf%n", &((*pconfig).sidechainrelativepos[leaf][last][self].y), &offset);
                buff+=offset;
                sscanf(buff, "%lf%n", &((*pconfig).sidechainrelativepos[leaf][last][self].z), &offset);
                buff+=offset;
                sscanf(buff, "%lf%n", &((*pconfig).sidechainnpar[leaf][last][self]), &offset);
                buff+=offset;
                sscanf(buff, "%lf%n", &((*pconfig).sidechainnphi[leaf][last][self]), &offset);
                buff+=offset;
            }
        }
    }
}

void initialize_bilayer_fragment(int Nchains, double_triple *pbox_dimension, int *chainlength, monomernodes **monomerid, coord *coordarray, bonded_params my_bonded_params, nonbonded_params my_nonbonded_params, bilayer config, int leafsymmetry, int xchains, double molarity){
    int j, l, i, k, lastpolar, lastpolartype, chaincounter=0, backboneid, sidechainid, selftype;
    double nz, phi, nx, ny, npar, nphi;
    double_triple position, nalong, nside, xvector, sidepos, nback, periodicbox;
    xvector.x=1;xvector.y=xvector.z=0;
	double spaceperchain=chainlength[0]*config.monomerspacing+config.terminusspacing;
	periodicbox.x=spaceperchain*xchains;
	periodicbox.y=(Nchains/2/xchains)*config.interchainspacing;
	periodicbox.z=config.boxheight;
    
	(*pbox_dimension).x=pow((1.*Nchains)/(molarity*molardensity), 1.0/3.0);
	(*pbox_dimension).y=(*pbox_dimension).z=(*pbox_dimension).x;
	if((periodicbox.x>(*pbox_dimension).x)||(periodicbox.y>(*pbox_dimension).y)) my_exit("molarity too high; sheet doesn't fit\n");
	
	if(Nchains%2!=0) my_exit("Nchains not divisible by 2!");
	int ychains=Nchains/2/xchains;
    for(j=0;j<2;j++){
		for(i=0;i<ychains;i++){
			for(l=0;l<xchains;l++){
				position.x=l*spaceperchain+(mod(i, 2))*config.offset-0.5*periodicbox.x+0.5*(*pbox_dimension).x;
				position.y=(i+0.5)*config.interchainspacing-0.5*periodicbox.y+0.5*(*pbox_dimension).y;
				position.z=0.5*(*pbox_dimension).z;
				for(k=0;k<chainlength[chaincounter];k++){
					backboneid=monomerid[chaincounter][k].backbone;
					sidechainid=monomerid[chaincounter][k].sidechain;
                    coordarray[backboneid].leafid=j;
                    coordarray[sidechainid].leafid=j;
                    if(coordarray[sidechainid].nodetype<0){             //  extra backbone, assume alternating amphiphilic and ending with nonpolar
                        lastpolartype=0;
                        selftype=lastpolar;
                    }
					else if(coordarray[sidechainid].nodetype>1){
						lastpolar=coordarray[sidechainid].nodetype;
						lastpolartype=0;                                //  self is polar
                        selftype=coordarray[sidechainid].nodetype;
					}
					else{
                        lastpolartype=lastpolar;
                        selftype=coordarray[sidechainid].nodetype;
                    }
					coordarray[backboneid].r=position;
					if((j==1)&&(leafsymmetry==-1)) coordarray[backboneid].r.x-=(-0.5*chainlength[chaincounter]+k)*config.monomerspacing;
					else coordarray[backboneid].r.x+=(-0.5*chainlength[chaincounter]+k)*config.monomerspacing;
					coordarray[backboneid].r=add_double_triple(coordarray[backboneid].r, config.backboneoffset[j][lastpolartype][selftype]);
					fmod_double_triple(&(coordarray[backboneid].r), (*pbox_dimension));
					nz=config.backbonenz[j][lastpolartype][selftype];
					coordarray[backboneid].n.z=(1-2*j)*((k%2)*2-1)*nz;          //  leaves facing and amphiphilic
					phi=config.backbonephi[j][lastpolartype][selftype];
					coordarray[backboneid].n.x=sqrt(1-nz*nz)*cos(phi);
					coordarray[backboneid].n.y=sqrt(1-nz*nz)*sin(phi);
					if(coordarray[sidechainid].nodetype>=0){
						nback=coordarray[backboneid].n;
						nalong=subtract_double_triple(xvector, scalar_multiply_double_triple(nback, nback.x));
						normalize(&nalong);
						nside=scalar_multiply_double_triple(cross_product(nback, nalong), (1-2*j));
						npar=config.sidechainnpar[j][lastpolartype][coordarray[sidechainid].nodetype];
						nphi=config.sidechainnphi[j][lastpolartype][coordarray[sidechainid].nodetype];
						sidepos=config.sidechainrelativepos[j][lastpolartype][coordarray[sidechainid].nodetype];
						coordarray[sidechainid].r=add_double_triple(add_double_triple(add_double_triple(coordarray[backboneid].r, scalar_multiply_double_triple(nback, sidepos.z)), scalar_multiply_double_triple(nalong, sidepos.x)), scalar_multiply_double_triple(nside, sidepos.y));
						fmod_double_triple(&(coordarray[sidechainid].r), (*pbox_dimension));
						nx=sqrt(1-pow(npar, 2))*cos(nphi);
						ny=sqrt(1-pow(npar, 2))*sin(nphi);
						coordarray[sidechainid].n=add_double_triple(add_double_triple(scalar_multiply_double_triple(nback, npar), scalar_multiply_double_triple(nalong, nx)), scalar_multiply_double_triple(nside, ny));
					}
				}
				chaincounter++;
				if(chaincounter==Nchains) return;
			}
		}
    }
}

void create_ordered_bilayer(bilayer *pconfig, int chainlength, int Nnodetypes, double residueshift, double chargedbackbonezoffset, double backbonezzigzag, double backbonenz, double backbonephi, bonded_params my_bonded_params, double leafspacing, double leafxoffsetfrac, double leafyoffsetfrac, int leafsymmetry){
	int i, j, k;
    
    (*pconfig).length=chainlength*(*pconfig).monomerspacing+(*pconfig).terminusspacing;
    (*pconfig).offset=(*pconfig).length*(*pconfig).offsetfraction+((*pconfig).monomerscaleoffsetfraction+residueshift)*(*pconfig).monomerspacing;
    
    allocate_3d_tensor_nozero(double_triple, (*pconfig).backboneoffset, 2, Nnodetypes, Nnodetypes);
    allocate_3d_tensor(double, (*pconfig).backbonenz, 2, Nnodetypes, Nnodetypes);
    allocate_3d_tensor(double, (*pconfig).backbonephi, 2, Nnodetypes, Nnodetypes);
    allocate_3d_tensor_nozero(double_triple, (*pconfig).sidechainrelativepos, 2, Nnodetypes, Nnodetypes);
    allocate_3d_tensor(double, (*pconfig).sidechainnpar, 2, Nnodetypes, Nnodetypes);
    allocate_3d_tensor(double, (*pconfig).sidechainnphi, 2, Nnodetypes, Nnodetypes);
    
    int leaf, charged, chargetype, last, self;
	for(leaf=0;leaf<2;leaf++){
        for(charged=0;charged<2;charged++){
            for(chargetype=2;chargetype<4;chargetype++){
                if(charged==1){
                    last=0;
                    self=chargetype;
                }
                else{
                    last=chargetype;
                    self=1;
                }
				(*pconfig).backboneoffset[leaf][last][self].x=0+leaf*leafxoffsetfrac*(*pconfig).monomerspacing;
				(*pconfig).backboneoffset[leaf][last][self].y=0+leaf*leafyoffsetfrac*(*pconfig).interchainspacing;
				(*pconfig).backboneoffset[leaf][last][self].z=(1-leaf*2)*chargedbackbonezoffset+(leaf-0.5)*leafspacing;
				if(self==1) (*pconfig).backboneoffset[leaf][last][self].z+=(1-leaf*2)*backbonezzigzag;
				(*pconfig).backbonenz[leaf][last][self]=backbonenz;
				if((leaf==1)&&(leafsymmetry==-1)){
					(*pconfig).backbonephi[leaf][last][self]=M_PI-backbonephi;
				}
				else{
					(*pconfig).backbonephi[leaf][last][self]=backbonephi;
				}
                (*pconfig).sidechainrelativepos[leaf][last][self]=rsol(my_bonded_params.sidechain[self]);
                if((leaf==1)&&(leafsymmetry==-1)) (*pconfig).sidechainrelativepos[leaf][last][self].x=-(*pconfig).sidechainrelativepos[leaf][last][self].x;
				(*pconfig).sidechainnpar[leaf][last][self]=nparhardsol(my_bonded_params.sidechainorientationtype[self], &((*pconfig).sidechainrelativepos[leaf][last][self]), my_bonded_params.sidechain[self]);
				if((i==1)&&(leafsymmetry==-1)){
					(*pconfig).sidechainnphi[leaf][last][self]=M_PI-nphihardsol(my_bonded_params.sidechainorientationtype[self], (*pconfig).sidechainnpar[leaf][last][self], (*pconfig).sidechainrelativepos[leaf][last][self], my_bonded_params.sidechain[self]);
				}
				else{
					(*pconfig).sidechainnphi[leaf][last][self]=nphihardsol(my_bonded_params.sidechainorientationtype[self], (*pconfig).sidechainnpar[leaf][last][self], (*pconfig).sidechainrelativepos[leaf][last][self], my_bonded_params.sidechain[self]);
				}
			}
        }
    }
}

void initialize_periodic_bilayer_xchains(int Nchains, double_triple *pbox_dimension, int *chainlength, monomernodes **monomerid, coord *coordarray, bonded_params my_bonded_params, nonbonded_params my_nonbonded_params, bilayer config, int leafsymmetry, int xchains, int randomize_xoffsets){
    int j, l, i, k, lastpolar, lastpolartype, chaincounter=0, backboneid, sidechainid, selftype;
    double nz, phi, nx, ny, npar, nphi, chainoffset;
    double_triple position, nalong, nside, xvector, sidepos, nback;
    xvector.x=1;xvector.y=xvector.z=0;
	double spaceperchain=chainlength[0]*config.monomerspacing+config.terminusspacing;
	(*pbox_dimension).x=spaceperchain*xchains;
	(*pbox_dimension).y=(Nchains/2/xchains)*config.interchainspacing;
	(*pbox_dimension).z=config.boxheight;
    if(Nchains%2!=0) my_exit("Nchains not divisible by 2!");
	int ychains=Nchains/2/xchains;
    for(j=0;j<2;j++){
		for(i=0;i<ychains;i++){
            if(randomize_xoffsets==1){
                if(i==0){
                    if(j==0) chainoffset=0;
                    else chainoffset=(rand_int(chainlength[0]/2+1)*2-chainlength[0]/4)*config.monomerspacing;
                }
                else{
                    chainoffset+=((rand_int(2)*2-1)*config.offset + (rand_int(chainlength[0]/2+1)*2-chainlength[0]/4)*config.monomerspacing);
                }
            }
            else{
                chainoffset=(mod(i, 2))*config.offset;
            }
			for(l=0;l<xchains;l++){
				position.x=l*spaceperchain+chainoffset;
				position.y=(i+0.5)*config.interchainspacing;
				position.z=0.5*(*pbox_dimension).z;
				for(k=0;k<chainlength[chaincounter];k++){
					backboneid=monomerid[chaincounter][k].backbone;
					sidechainid=monomerid[chaincounter][k].sidechain;
                    coordarray[backboneid].leafid=j;
                    coordarray[sidechainid].leafid=j;
                    if(coordarray[sidechainid].nodetype<0){             //  extra backbone, assume alternating amphiphilic and ending with nonpolar
                        lastpolartype=0;
                        selftype=lastpolar;
                    }
					else if(coordarray[sidechainid].nodetype>1){
						lastpolar=coordarray[sidechainid].nodetype;
						lastpolartype=0;                                //  self is polar
                        selftype=coordarray[sidechainid].nodetype;
					}
					else{
                        lastpolartype=lastpolar;
                        selftype=coordarray[sidechainid].nodetype;
                    }
					coordarray[backboneid].r=position;
					if((j==1)&&(leafsymmetry==-1)) coordarray[backboneid].r.x-=(-0.5*chainlength[chaincounter]+k)*config.monomerspacing;
					else coordarray[backboneid].r.x+=(-0.5*chainlength[chaincounter]+k)*config.monomerspacing;
					coordarray[backboneid].r=add_double_triple(coordarray[backboneid].r, config.backboneoffset[j][lastpolartype][selftype]);
					fmod_double_triple(&(coordarray[backboneid].r), (*pbox_dimension));
					nz=config.backbonenz[j][lastpolartype][selftype];
					coordarray[backboneid].n.z=(1-2*j)*((k%2)*2-1)*nz;          //  leaves facing and amphiphilic
					phi=config.backbonephi[j][lastpolartype][selftype];
					coordarray[backboneid].n.x=sqrt(1-nz*nz)*cos(phi);
					coordarray[backboneid].n.y=sqrt(1-nz*nz)*sin(phi);
					if(coordarray[sidechainid].nodetype>=0){
						nback=coordarray[backboneid].n;
						nalong=subtract_double_triple(xvector, scalar_multiply_double_triple(nback, nback.x));
						normalize(&nalong);
						nside=scalar_multiply_double_triple(cross_product(nback, nalong), (1-2*j));
						npar=config.sidechainnpar[j][lastpolartype][coordarray[sidechainid].nodetype];
						nphi=config.sidechainnphi[j][lastpolartype][coordarray[sidechainid].nodetype];
						sidepos=config.sidechainrelativepos[j][lastpolartype][coordarray[sidechainid].nodetype];
						coordarray[sidechainid].r=add_double_triple(add_double_triple(add_double_triple(coordarray[backboneid].r, scalar_multiply_double_triple(nback, sidepos.z)), scalar_multiply_double_triple(nalong, sidepos.x)), scalar_multiply_double_triple(nside, sidepos.y));
						fmod_double_triple(&(coordarray[sidechainid].r), (*pbox_dimension));
						nx=sqrt(1-pow(npar, 2))*cos(nphi);
						ny=sqrt(1-pow(npar, 2))*sin(nphi);
						coordarray[sidechainid].n=add_double_triple(add_double_triple(scalar_multiply_double_triple(nback, npar), scalar_multiply_double_triple(nalong, nx)), scalar_multiply_double_triple(nside, ny));
					}
				}
				chaincounter++;
				if(chaincounter==Nchains) return;
			}
		}
    }
}

void input_monolayer(char *configname, monolayer *pconfig, int chainlength, int Nnodetypes, double residueshift){
 	FILE *fp;
	fp=fopen(configname, "r");
	int size=1000, i, j;
	char *buff;
	buff=xcalloc(size, sizeof(char));
	while(fgets(buff, size, fp));
	fclose(fp);
	int offset;
	sscanf(buff, "%*i%n", &offset);
	buff+=offset;
	sscanf(buff, "%*f%n", &offset);
	buff+=offset;
	sscanf(buff, "%lf%n", &((*pconfig).monomerspacing), &offset);
	buff+=offset;
	sscanf(buff, "%lf%n", &((*pconfig).monomerscaleoffsetfraction), &offset);
	buff+=offset;
	sscanf(buff, "%lf%n", &((*pconfig).terminusspacing), &offset);
	buff+=offset;
	sscanf(buff, "%lf%n", &((*pconfig).interchainspacing), &offset);
	buff+=offset;
    
    (*pconfig).length=chainlength*(*pconfig).monomerspacing+(*pconfig).terminusspacing;
    (*pconfig).offset=(*pconfig).length*(*pconfig).offsetfraction+((*pconfig).monomerscaleoffsetfraction+residueshift)*(*pconfig).monomerspacing;
    
	if(Nnodetypes<4) Nnodetypes=4;
    
    allocate_matrix_nozero(double_triple, ((*pconfig).backboneoffset), Nnodetypes, Nnodetypes);
    allocate_matrix(double, (*pconfig).backbonenz, Nnodetypes, Nnodetypes);
    allocate_matrix(double, (*pconfig).backbonephi, Nnodetypes, Nnodetypes);
    allocate_matrix_nozero(double_triple, (*pconfig).sidechainrelativepos, Nnodetypes, Nnodetypes);
    allocate_matrix(double, (*pconfig).sidechainnpar, Nnodetypes, Nnodetypes);
    allocate_matrix(double, (*pconfig).sidechainnphi, Nnodetypes, Nnodetypes);
    
    int charged, chargetype, last, self;
    for(charged=0;charged<2;charged++){
        for(chargetype=2;chargetype<4;chargetype++){
            if(charged==1){
                last=0;
                self=chargetype;
            }
            else{
                last=chargetype;
                self=1;
            }
            sscanf(buff, "%lf%n", &((*pconfig).backboneoffset[last][self].x), &offset);
            buff+=offset;
            sscanf(buff, "%lf%n", &((*pconfig).backboneoffset[last][self].y), &offset);
            buff+=offset;
            sscanf(buff, "%lf%n", &((*pconfig).backboneoffset[last][self].z), &offset);
            buff+=offset;
            sscanf(buff, "%lf%n", &((*pconfig).backbonenz[last][self]), &offset);
            buff+=offset;
            sscanf(buff, "%lf%n", &((*pconfig).backbonephi[last][self]), &offset);
            buff+=offset;
            sscanf(buff, "%lf%n", &((*pconfig).sidechainrelativepos[last][self].x), &offset);
            buff+=offset;
            sscanf(buff, "%lf%n", &((*pconfig).sidechainrelativepos[last][self].y), &offset);
            buff+=offset;
            sscanf(buff, "%lf%n", &((*pconfig).sidechainrelativepos[last][self].z), &offset);
            buff+=offset;
            sscanf(buff, "%lf%n", &((*pconfig).sidechainnpar[last][self]), &offset);
            buff+=offset;
            sscanf(buff, "%lf%n", &((*pconfig).sidechainnphi[last][self]), &offset);
            buff+=offset;
        }
    }
}

double_triple rsol(sidechain_params sidechain){
    double arctan=quartic_global_minimum(sidechain.J14, sidechain.J13, sidechain.J12, sidechain.J11);
    double rparoverrperp=tan(arctan);
    
    //  solve quadratic
    
    double a=pow(rparoverrperp, 2)+1;
    double b=-2*(rparoverrperp*sidechain.rpar0+sidechain.rperp0);
    double c=pow(sidechain.rperp0, 2)+pow(sidechain.rpar0, 2)-pow(sidechain.r0, 2);
    double rperp=(-b+sqrt(b*b-4*a*c))/(2*a);
    
    double rpar=rparoverrperp*rperp;
    double_triple result;
    result.z=rpar;
    result.x=rperp;
    result.y=0;
    return result;
}

double nparhardsol(int orientationtype, double_triple *pr, sidechain_params sidechain){
    if(orientationtype==1) return 0;
    double result=sidechain.r20+(*pr).z*sidechain.r21+pow((*pr).z, 2)*sidechain.r22, rperpdif, rpardif;
    if(result>1){
        printf("changing sidechainrelativepos.z from %f to", (*pr).z);
        (*pr).z=(-sidechain.r21+sqrt(pow(sidechain.r21, 2)-4*sidechain.r22*(sidechain.r20-1)))/(2*sidechain.r22);
        rpardif=(*pr).z-sidechain.rpar0;
        rperpdif=sqrt(pow(sidechain.r0, 2)-pow(rpardif, 2));
        (*pr).x=sidechain.rperp0+rperpdif;
        printf(" %f\n", (*pr).z);
        return 1;
    }
    else if(result<-1){
        printf("changing sidechainrelativepos.z from %f to", (*pr).z);
        (*pr).z=(-sidechain.r21+sqrt(pow(sidechain.r21, 2)-4*sidechain.r22*(sidechain.r20-1)))/(2*sidechain.r22);
        rpardif=(*pr).z-sidechain.rpar0;
        rperpdif=sqrt(pow(sidechain.r0, 2)-pow(rpardif, 2));
        (*pr).x=sidechain.rperp0+rperpdif;
        printf(" %f\n", (*pr).z);
        return -1;
    }
    return result;
}

double nphihardsol(int orientationtype, double npar, double_triple relativepos, sidechain_params sidechain){
    double rpar=relativepos.z;
    double rperp=relativepos.x;
    double rdotn, rdotnperphat;//, rmag=norm(relativepos);
    double nparmag=fabs(npar);
    if(orientationtype==0){
        if(nparmag==1){
            rdotn=0;
            return 0;           //  phi is unconstrained; choose 0 arbitrarily
        }
        rdotn=sidechain.r30+npar*sidechain.r31+npar*npar*sidechain.r32;
    }
    else if(orientationtype==1){
        double p4, p3, p2, p1;
        if(nparmag==1){
            rdotn=0;
            return 0;           //  phi is unconstrained; choose 0 arbitrarily
        }
        else{
            
            //  find minimum of quartic
            
            p4=sidechain.J340;
            p3=sidechain.J331*nparmag;
            p2=(sidechain.J320+sidechain.J322*nparmag*nparmag);
            p1=sidechain.J311*nparmag+sidechain.J313*nparmag*nparmag*nparmag;
            rdotn=quartic_global_minimum(p4, p3, p2, p1);
            if(npar<0) rdotn=-rdotn;
        }
    }
    rdotnperphat=(rdotn-rpar*npar)/(rperp*sqrt(1-npar*npar));
    if(rdotnperphat>1){
        rdotnperphat=1;
        return 0;
        
        //  can't adjust relativepos.x to fix this; would have to make it bigger
        
    }
    else if(rdotnperphat<-1){
        rdotnperphat=-1;
        
        //  can't adjust relativepos.x to fix this; would have to make it bigger
        
        return M_PI;
    }
    return acos(rdotnperphat);
}

void initialize_periodic_monolayer_xchains(int Nchains, double_triple *pbox_dimension, int *chainlength, monomernodes **monomerid, coord *coordarray, bonded_params my_bonded_params, nonbonded_params my_nonbonded_params, monolayer config, int xchains, int randomize_xoffsets){
    int i, j, k, lastpolar=2, lastpolartype, chaincounter=0, backboneid, sidechainid, selftype;
    double nz, phi, nx, ny, npar, nphi, chainoffset;
    double_triple position, nalong, nside, xvector, sidepos, nback;
    xvector.x=1;xvector.y=xvector.z=0;
	double spaceperchain=chainlength[0]*config.monomerspacing+config.terminusspacing;
	(*pbox_dimension).x=spaceperchain*xchains;
	(*pbox_dimension).y=(Nchains/xchains)*config.interchainspacing;
	(*pbox_dimension).z=config.boxheight;
	
	int ychains=Nchains/xchains;
	for(i=0;i<ychains;i++){
        if(randomize_xoffsets==1){
            if(i==0) chainoffset=0;
            else{
                chainoffset+=((rand_int(2)*2-1)*config.offset + (rand_int(chainlength[0]/2+1)*2-chainlength[0]/4)*config.monomerspacing);
            }
        }
        else{
            chainoffset=(mod(i, 2))*config.offset;
        }
		for(j=0;j<xchains;j++){
			position.x=j*spaceperchain+chainoffset;
			position.y=(i+0.5)*config.interchainspacing;
			position.z=(*pbox_dimension).z-0.5*my_nonbonded_params.vacuumthickness;
			for(k=0;k<chainlength[chaincounter];k++){
				backboneid=monomerid[chaincounter][k].backbone;
				sidechainid=monomerid[chaincounter][k].sidechain;
                if(coordarray[sidechainid].nodetype<0){             //  extra backbone, assume alternating amphiphilic and ending with nonpolar
                    lastpolartype=0;
                    selftype=lastpolar;
                }
				else if(coordarray[sidechainid].nodetype>1){
					lastpolar=coordarray[sidechainid].nodetype;
					lastpolartype=0;                                //  self is polar
                    selftype=coordarray[sidechainid].nodetype;
				}
				else{
                    lastpolartype=lastpolar;
                    selftype=coordarray[sidechainid].nodetype;
                }
				coordarray[backboneid].r=position;
				coordarray[backboneid].r.x+=(-0.5*chainlength[chaincounter]+k)*config.monomerspacing;
				coordarray[backboneid].r=add_double_triple(coordarray[backboneid].r, config.backboneoffset[lastpolartype][selftype]);
				
				fmod_double_triple(&(coordarray[backboneid].r), (*pbox_dimension));
				nz=config.backbonenz[lastpolartype][selftype];
				coordarray[backboneid].n.z=((k%2)*2-1)*nz;          //  amphiphilic
				phi=config.backbonephi[lastpolartype][selftype];
				coordarray[backboneid].n.x=sqrt(1-nz*nz)*cos(phi);
				coordarray[backboneid].n.y=sqrt(1-nz*nz)*sin(phi);
				if(coordarray[sidechainid].nodetype>=0){
					nback=coordarray[backboneid].n;
					nalong=subtract_double_triple(xvector, scalar_multiply_double_triple(nback, nback.x));
					normalize(&nalong);
					nside=cross_product(nback, nalong);
					npar=config.sidechainnpar[lastpolartype][coordarray[sidechainid].nodetype];
					nphi=config.sidechainnphi[lastpolartype][coordarray[sidechainid].nodetype];
					sidepos=config.sidechainrelativepos[lastpolartype][coordarray[sidechainid].nodetype];
					coordarray[sidechainid].r=add_double_triple(add_double_triple(add_double_triple(coordarray[backboneid].r, scalar_multiply_double_triple(nback, sidepos.z)), scalar_multiply_double_triple(nalong, sidepos.x)), scalar_multiply_double_triple(nside, sidepos.y));
					fmod_double_triple(&(coordarray[sidechainid].r), (*pbox_dimension));
					nx=sqrt(1-pow(npar, 2))*cos(nphi);
					ny=sqrt(1-pow(npar, 2))*sin(nphi);
					coordarray[sidechainid].n=add_double_triple(add_double_triple(scalar_multiply_double_triple(nback, npar), scalar_multiply_double_triple(nalong, nx)), scalar_multiply_double_triple(nside, ny));
				}
			}
			chaincounter++;
			if(chaincounter==Nchains) return;
		}
	}
}

void initialize_periodic_monolayer_xchains_bothinterfaces(int Nchains, double_triple *pbox_dimension, int *chainlength, monomernodes **monomerid, coord *coordarray, bonded_params my_bonded_params, nonbonded_params my_nonbonded_params, monolayer config, int xchains, int randomize_xoffsets){
    int i, j, k, interface, lastpolar=2, lastpolartype, chaincounter=0, backboneid, sidechainid, firstdirection, selftype;
    double nz, phi, nx, ny, npar, nphi, chainoffset;
    double_triple position, nalong, nside, xvector, sidepos, nback;
    xvector.x=1;xvector.y=xvector.z=0;
	double spaceperchain=chainlength[0]*config.monomerspacing+config.terminusspacing;
	(*pbox_dimension).x=spaceperchain*xchains;
	int ychains=Nchains/xchains/2;
	(*pbox_dimension).y=ychains*config.interchainspacing;
	(*pbox_dimension).z=config.boxheight;
	
    for(interface=0;interface<2;interface++){
        for(i=0;i<ychains;i++){
            if(randomize_xoffsets==1){
                if(i==0) chainoffset=0;
                else{
                    chainoffset+=((rand_int(2)*2-1)*config.offset + (rand_int(chainlength[0]/2+1)*2-chainlength[0]/4)*config.monomerspacing);
                }
            }
            else{
                chainoffset=(mod(i, 2))*config.offset;
            }
            for(j=0;j<xchains;j++){
                //position.x=j*spaceperchain+(mod(i, 2))*config.offset;
                position.x=j*spaceperchain+chainoffset;
                position.y=(i+0.5)*config.interchainspacing;
                if(interface==0) position.z=(*pbox_dimension).z-0.5*my_nonbonded_params.vacuumthickness;
                else position.z=0.5*my_nonbonded_params.vacuumthickness;
                for(k=0;k<chainlength[chaincounter];k++){
                    backboneid=monomerid[chaincounter][k].backbone;
                    sidechainid=monomerid[chaincounter][k].sidechain;
                    coordarray[backboneid].leafid=interface;
                    coordarray[sidechainid].leafid=interface;
                    if(coordarray[sidechainid].nodetype==1) firstdirection=1;           //  phenyl toward air
                    else firstdirection=-1;
                    if(coordarray[sidechainid].nodetype<0){             //  extra backbone, assume alternating amphiphilic and ending with nonpolar
                        lastpolartype=0;
                        selftype=lastpolar;
                    }
                    else if(coordarray[sidechainid].nodetype>1){
                        lastpolar=coordarray[sidechainid].nodetype;
                        lastpolartype=0;                                //  self is polar
                        selftype=coordarray[sidechainid].nodetype;
                    }
                    else{
                        lastpolartype=lastpolar;
                        selftype=coordarray[sidechainid].nodetype;
                    }
                    coordarray[backboneid].r=position;
                    coordarray[backboneid].r.x+=(-0.5*chainlength[chaincounter]+k)*config.monomerspacing;
                    coordarray[backboneid].r=add_double_triple(coordarray[backboneid].r, config.backboneoffset[lastpolartype][selftype]);
                    
                    fmod_double_triple(&(coordarray[backboneid].r), (*pbox_dimension));
                    nz=(2*interface-1)*firstdirection*config.backbonenz[lastpolartype][selftype];
                    coordarray[backboneid].n.z=((k%2)*2-1)*nz;          //  amphiphilic
                    phi=config.backbonephi[lastpolartype][selftype];
                    coordarray[backboneid].n.x=sqrt(1-nz*nz)*cos(phi);
                    coordarray[backboneid].n.y=sqrt(1-nz*nz)*sin(phi);
					if(coordarray[sidechainid].nodetype>=0){
						nback=coordarray[backboneid].n;
						nalong=subtract_double_triple(xvector, scalar_multiply_double_triple(nback, nback.x));
						normalize(&nalong);
						nside=cross_product(nback, nalong);
						npar=config.sidechainnpar[lastpolartype][coordarray[sidechainid].nodetype];
						nphi=config.sidechainnphi[lastpolartype][coordarray[sidechainid].nodetype];
						sidepos=config.sidechainrelativepos[lastpolartype][coordarray[sidechainid].nodetype];
						coordarray[sidechainid].r=add_double_triple(add_double_triple(add_double_triple(coordarray[backboneid].r, scalar_multiply_double_triple(nback, sidepos.z)), scalar_multiply_double_triple(nalong, sidepos.x)), scalar_multiply_double_triple(nside, sidepos.y));
						fmod_double_triple(&(coordarray[sidechainid].r), (*pbox_dimension));
						nx=sqrt(1-pow(npar, 2))*cos(nphi);
						ny=sqrt(1-pow(npar, 2))*sin(nphi);
						coordarray[sidechainid].n=add_double_triple(add_double_triple(scalar_multiply_double_triple(nback, npar), scalar_multiply_double_triple(nalong, nx)), scalar_multiply_double_triple(nside, ny));
					}
				}
				chaincounter++;
                if(chaincounter==Nchains) return;
            }
        }
    }
}

void initialize_2dxylattice_fromconfig(int Nmonomers, int Nchains, double density, double_triple *pbox_dimension, int *chainlength, monomernodes **monomerid, coord *coordarray, bonded_params my_bonded_params, nonbonded_params my_nonbonded_params, polymer config){
	int chainsonaside, i, j, k, lastpolar, lastpolartype, chaincounter=0, backboneid, sidechainid, selftype;
    double nz, phi, nx, ny, npar, nphi;
    double_triple position, nalong, nside, xvector, sidepos, nback;
    xvector.x=1;xvector.y=xvector.z=0;
	(*pbox_dimension).x=pow(Nmonomers/density, 1.0/3.0);
	(*pbox_dimension).y=(*pbox_dimension).z=(*pbox_dimension).x;
	chainsonaside=(int) floor(sqrt(1.*Nchains))+1;
	double chainspacing=(*pbox_dimension).x/chainsonaside;
	position.z=(*pbox_dimension).z-0.5*my_nonbonded_params.vacuumthickness;
    for(i=0;i<chainsonaside;i++){
		for(j=0;j<chainsonaside;j++){
			position.x=(i+0.5)*chainspacing;
			position.y=(j+0.5)*chainspacing;
			for(k=0;k<chainlength[chaincounter];k++){
				backboneid=monomerid[chaincounter][k].backbone;
				sidechainid=monomerid[chaincounter][k].sidechain;
                if(coordarray[sidechainid].nodetype<0){             //  extra backbone, assume alternating amphiphilic and ending with nonpolar
                    lastpolartype=0;
                    selftype=lastpolar;
                }
				else if(coordarray[sidechainid].nodetype>1){
					lastpolar=coordarray[sidechainid].nodetype;
					lastpolartype=0;                                //  self is polar
                    selftype=coordarray[sidechainid].nodetype;
				}
				else{
                    lastpolartype=lastpolar;
                    selftype=coordarray[sidechainid].nodetype;
                }
				coordarray[backboneid].r=position;
				coordarray[backboneid].r.x+=(-0.5*chainlength[chaincounter]+k)*config.monomerspacing;
				coordarray[backboneid].r=add_double_triple(coordarray[backboneid].r, config.backboneoffset[lastpolartype][selftype]);
				fmod_double_triple(&(coordarray[backboneid].r), (*pbox_dimension));
				nz=config.backbonenz[lastpolartype][selftype];
				coordarray[backboneid].n.z=((k%2)*2-1)*nz;          //  amphiphilic
				phi=config.backbonephi[lastpolartype][selftype];
				coordarray[backboneid].n.x=sqrt(1-nz*nz)*cos(phi);
				coordarray[backboneid].n.y=sqrt(1-nz*nz)*sin(phi);
				if(coordarray[sidechainid].nodetype>=0){
					nback=coordarray[backboneid].n;
					nalong=subtract_double_triple(xvector, scalar_multiply_double_triple(nback, nback.x));
					normalize(&nalong);
					nside=cross_product(nback, nalong);
					npar=config.sidechainnpar[lastpolartype][coordarray[sidechainid].nodetype];
					nphi=config.sidechainnphi[lastpolartype][coordarray[sidechainid].nodetype];
					sidepos=config.sidechainrelativepos[lastpolartype][coordarray[sidechainid].nodetype];
					coordarray[sidechainid].r=add_double_triple(add_double_triple(add_double_triple(coordarray[backboneid].r, scalar_multiply_double_triple(nback, sidepos.z)), scalar_multiply_double_triple(nalong, sidepos.x)), scalar_multiply_double_triple(nside, sidepos.y));
					fmod_double_triple(&(coordarray[sidechainid].r), (*pbox_dimension));
					nx=sqrt(1-pow(npar, 2))*cos(nphi);
					ny=sqrt(1-pow(npar, 2))*sin(nphi);
					coordarray[sidechainid].n=add_double_triple(add_double_triple(scalar_multiply_double_triple(nback, npar), scalar_multiply_double_triple(nalong, nx)), scalar_multiply_double_triple(nside, ny));
				}
			}
			chaincounter++;
			if(chaincounter==Nchains) return;
		}
	}
}

void input_coords(char *filename, int *pNnodes, int *pNchains, int **chainlength, monomernodes ***monomerid, coord **coordarray, double_triple *pbox_dimension, double *pvacuumthickness, int *pNnodetypes, monomertypes *pmonomercount, long long int *pt, int *pframe, int *prunningtime, double_triple **pmeshpositions, triangledata **pmeshtriangledata, bonddata **pmeshbonddata, vertexdata **pmeshdata, solvation_parameters *psolvp){
	int i, j;
    double val;
	FILE *inp;
	inp=fopen(filename, "r");
	fscanf(inp, "%lf %lf %lf", &((*pbox_dimension).x), &((*pbox_dimension).y), &((*pbox_dimension).z));
	fscanf(inp, "%i", &((*psolvp).interface));
	if(((*psolvp).interface)==1)fscanf(inp, " %lf", pvacuumthickness);
	fscanf(inp, "%i %i %i %i %i %i", &(*pNnodes), &(*pNchains), &((*pmonomercount).charged), &((*pmonomercount).nonpolar), &((*pmonomercount).monomers), &(*pNnodetypes));
	for(i=0;i<(*pNnodetypes);i++) fscanf(inp, " %i", &((*pmonomercount).type[i]));
	(*coordarray)=xcalloc((*pNnodes), sizeof(coord));
	(*chainlength)=xcalloc((*pNchains), sizeof(int));
	(*monomerid)=xcalloc((*pNchains), sizeof(monomernodes *));
	for(i=0;i<*pNchains;i++){
		fscanf(inp, "%i", &((*chainlength)[i]));
		((*monomerid)[i])=xcalloc(((*chainlength)[i]), sizeof(monomernodes));
		for(j=0;j<(*chainlength)[i];j++){
			fscanf(inp, "%i %i", &((*monomerid)[i][j].backbone), &((*monomerid)[i][j].sidechain));
		}
	}
	for(i=0;i<(*pNnodes);i++){
		fscanf(inp, "%lf %lf %lf %lf %lf %lf %i %i %i %i %i", &((*coordarray)[i].r.x), &((*coordarray)[i].r.y), &((*coordarray)[i].r.z), &((*coordarray)[i].n.x), &((*coordarray)[i].n.y), &((*coordarray)[i].n.z), &((*coordarray)[i].nodetype), &((*coordarray)[i].chainid), &((*coordarray)[i].monomerid), &((*coordarray)[i].orientationtype), &((*coordarray)[i].leafid));
	}
    fscanf(inp, "%lli %i %i", &(*pt), &(*pframe), &(*prunningtime));
	if(((*psolvp).interface)==1){
		fscanf(inp, "%i %i", &(((*psolvp).meshsize).x), &(((*psolvp).meshsize).y));
		((*psolvp).meshnumber)=2*((*psolvp).meshsize).x*((*psolvp).meshsize).y;
		allocate_array(double_triple, (*pmeshpositions), ((*psolvp).meshnumber));
		allocate_array(triangledata, (*pmeshtriangledata), (2*(*psolvp).meshnumber));
		allocate_array(bonddata, (*pmeshbonddata), (3*(*psolvp).meshnumber));
		allocate_array(vertexdata, (*pmeshdata), ((*psolvp).meshnumber));
		fscanf(inp, "%lf", &((*psolvp).meshbondlength));
		fscanf(inp, "%lf", &((*psolvp).maxmeshbond));
		fscanf(inp, "%lf", &((*psolvp).longestmeshbond));
		fscanf(inp, "%lf", &((*psolvp).totalvolume));
		fscanf(inp, "%lf", &((*psolvp).targetvolume));
		fscanf(inp, "%lf %lf", &((*psolvp).totallength[0]), &((*psolvp).totallength[1]));
        fscanf(inp, "%lf", &val);
        if(val>0){
            (*psolvp).totallength[2]=val;
            fscanf(inp, "%lf %lf %lf", &((*psolvp).totallength[3]), &((*psolvp).totallength[4]), &((*psolvp).totallength[5]));
            j=0;
        }
        else{
            i=0;
            j=0;
            (*pmeshpositions)[i*((*psolvp).meshsize).x+j].x=val;
            fscanf(inp, "%lf %lf", &((*pmeshpositions)[i*((*psolvp).meshsize).x+j].y), &((*pmeshpositions)[i*((*psolvp).meshsize).x+j].z));
            j=1;
        }
		for(i=0;i<((*psolvp).meshsize).y;i++){
			for(;j<((*psolvp).meshsize).x;j++){
				fscanf(inp, "%lf %lf %lf", &((*pmeshpositions)[i*((*psolvp).meshsize).x+j].x), &((*pmeshpositions)[i*((*psolvp).meshsize).x+j].y), &((*pmeshpositions)[i*((*psolvp).meshsize).x+j].z));
			}
            j=0;
		}        
		for(i=0;i<((*psolvp).meshsize).y;i++){
			for(j=0;j<((*psolvp).meshsize).x;j++){
				fscanf(inp, "%lf %lf %lf", &((*pmeshpositions)[(((*psolvp).meshsize).y+i)*((*psolvp).meshsize).x+j].x), &((*pmeshpositions)[(((*psolvp).meshsize).y+i)*((*psolvp).meshsize).x+j].y), &((*pmeshpositions)[(((*psolvp).meshsize).y+i)*((*psolvp).meshsize).x+j].z));
			}
		}
        for(i=0;i<2*(*psolvp).meshnumber;i++){
            fscanf(inp, "%lf %lf %lf %lf %lf", &((*pmeshtriangledata)[i].area), &((*pmeshtriangledata)[i].columnvolume), &((*pmeshtriangledata)[i].normal.x), &((*pmeshtriangledata)[i].normal.y), &((*pmeshtriangledata)[i].normal.z));
        }
        for(i=0;i<3*(*psolvp).meshnumber;i++){
            fscanf(inp, "%lf %lf", &((*pmeshbonddata)[i].length), &((*pmeshbonddata)[i].meancurvature));
        }
		for(i=0;i<(*psolvp).meshnumber;i++){
            fscanf(inp, "%lf %lf %lf", &((*pmeshdata)[i].area), &((*pmeshdata)[i].meancurvature), &((*pmeshdata)[i].meancurvaturesquared));
        }
	}
}

void recalculate_meshbonddata_meshdata(int meshnumber, double_triple *meshpositions, bonddata *meshbonddata, vertexdata *meshdata, triangledata *meshtriangledata, double_triple box_dimension, int **mapbondstomesh, int **mapbondstotriangles, int **mapmeshtotriangles, int **mapmeshtobonds, double *totallength){
    int i, j, k;
    double bondlength;
    double_triple bond;
    for(i=0;i<3*meshnumber;i++){
        bond=subtract_double_triple(meshpositions[mapbondstomesh[i][0]], meshpositions[mapbondstomesh[i][1]]);
        recenter_double_triple(&bond, box_dimension);
        bondlength=norm(bond);
        meshbonddata[i].length=bondlength;
        meshbonddata[i].meancurvature=bondlength*mysafeacos(dot_product(meshtriangledata[mapbondstotriangles[i][0]].normal, meshtriangledata[mapbondstotriangles[i][1]].normal));
    }
    for(i=0;i<meshnumber;i++){
        meshdata[i].area=0;
        meshdata[i].meancurvature=0;
        for(j=0;j<6;j++){
            meshdata[i].area+=meshtriangledata[mapmeshtotriangles[i][j]].area;
            meshdata[i].meancurvature+=meshbonddata[mapmeshtobonds[i][j]].meancurvature;
        }
        meshdata[i].meancurvaturesquared=pow(meshdata[i].meancurvature, 2)/meshdata[i].area;
    }
    for(i=0;i<6;i++) totallength[i]=0;
    for(i=0;i<2;i++){
        for(j=0;j<meshnumber/2;j++){
            for(k=0;k<3;k++){
                bondlength=meshbonddata[mapmeshtobonds[i*meshnumber/2+j][k]].length;
                totallength[i*3+k]+=bondlength;
            }
        }
    }
}

void input_coords_and_replicate_stacks(char *filename, int *pNnodes, int *pNchains, int **chainlength, monomernodes ***monomerid, coord **coordarray, double_triple *pbox_dimension, int *pinterface, double *pvacuumthickness, int *pNnodetypes, monomertypes *pmonomercount, long long int *pt, int *pframe, int *prunningtime, int replicates, double_triple displacement, double minsep){
	int i, j, Nchainsinput, Nnodesinput, r, replicateindex, count=0;
    double dif, av=0, max=0;
	FILE *inp;
	inp=fopen(filename, "r");
	fscanf(inp, "%lf %lf %lf", &((*pbox_dimension).x), &((*pbox_dimension).y), &((*pbox_dimension).z));
	fscanf(inp, "%i", &(*pinterface));
	if((*pinterface)==1)fscanf(inp, " %lf", pvacuumthickness);
	fscanf(inp, "%i %i %i %i %i %i", &Nnodesinput, &Nchainsinput, &((*pmonomercount).charged), &((*pmonomercount).nonpolar), &((*pmonomercount).monomers), &(*pNnodetypes));
    (*pNnodes)=replicates*Nnodesinput;
    (*pNchains)=replicates*Nchainsinput;
    ((*pmonomercount).charged)*=replicates;
    ((*pmonomercount).nonpolar)*=replicates;
    ((*pmonomercount).monomers)*=replicates;
	for(i=0;i<(*pNnodetypes);i++){
        fscanf(inp, " %i", &((*pmonomercount).type[i]));
        ((*pmonomercount).type[i])*=replicates;
    }
	(*coordarray)=xcalloc((*pNnodes), sizeof(coord));
	(*chainlength)=xcalloc((*pNchains), sizeof(int));
	(*monomerid)=xcalloc((*pNchains), sizeof(monomernodes *));
	for(i=0;i<Nchainsinput;i++){
		fscanf(inp, "%i", &((*chainlength)[i]));
		((*monomerid)[i])=xcalloc(((*chainlength)[i]), sizeof(monomernodes));
		for(j=0;j<(*chainlength)[i];j++){
			fscanf(inp, "%i %i", &((*monomerid)[i][j].backbone), &((*monomerid)[i][j].sidechain));
		}
	}
	for(i=0;i<Nnodesinput;i++){
		fscanf(inp, "%lf %lf %lf %lf %lf %lf %i %i %i %i %i", &((*coordarray)[i].r.x), &((*coordarray)[i].r.y), &((*coordarray)[i].r.z), &((*coordarray)[i].n.x), &((*coordarray)[i].n.y), &((*coordarray)[i].n.z), &((*coordarray)[i].nodetype), &((*coordarray)[i].chainid), &((*coordarray)[i].monomerid), &((*coordarray)[i].orientationtype), &((*coordarray)[i].leafid));
        if((*coordarray)[i].nodetype>=0){                           //  not empty
            dif=(*coordarray)[i].r.z-(*coordarray)[0].r.z;
            recenter(dif, ((*pbox_dimension).z));
            av+=dif;
            count++;
        }
	}
    av=av/(1.*count)+(*coordarray)[0].r.z;
    recenter(av, ((*pbox_dimension).z));
    for(i=0;i<Nnodesinput;i++){
        if((*coordarray)[i].nodetype>=0){                           //  not empty
            dif=(*coordarray)[i].r.z-av;
            recenter(dif, ((*pbox_dimension).z));
            if(dif*dif>max) max=dif*dif;
        }
    }
    max=sqrt(max);
	printf("max displacement 2*%f+%f=%f\n", max, minsep, 2*max+minsep);
    if(2*max+minsep>displacement.z){
        printf("intersheet z displacement %f not bigger than 2*%f+%f!\n", displacement.z, max, minsep);
        exit(1);
    }
    
    ((*pbox_dimension).z)=replicates*displacement.z;
    
    for(i=0;i<Nnodesinput;i++){
        fmod_double_triple(&((*coordarray)[i].r), (*pbox_dimension));
    }
    
    replicateindex=Nnodesinput;
    for(r=1;r<replicates;r++){
        for(i=0;i<Nchainsinput;i++){
            (*chainlength)[r*Nchainsinput+i]=(*chainlength)[i];
            ((*monomerid)[r*Nchainsinput+i])=xcalloc(((*chainlength)[i]), sizeof(monomernodes));
            for(j=0;j<(*chainlength)[i];j++){
                (*coordarray)[replicateindex]=(*coordarray)[(*monomerid)[i][j].backbone];
                if(i!=(*coordarray)[replicateindex].chainid){
                    printf("input chainids %i and %i don't agree!\n", i, (*coordarray)[replicateindex].chainid);
                    exit(1);
                }
                (*coordarray)[replicateindex].chainid=i+r*Nchainsinput;
				(*coordarray)[replicateindex].leafid+=2*r;
                (*monomerid)[i+r*Nchainsinput][j].backbone=replicateindex;
                (*coordarray)[replicateindex].r=add_double_triple((*coordarray)[replicateindex].r, scalar_multiply_double_triple(displacement, r));
                fmod_double_triple(&((*coordarray)[replicateindex].r), (*pbox_dimension));
                replicateindex++;
                (*coordarray)[replicateindex]=(*coordarray)[(*monomerid)[i][j].sidechain];
                if(i!=(*coordarray)[replicateindex].chainid){
                    printf("input chainids %i and %i don't agree!\n", i, (*coordarray)[replicateindex].chainid);
                    exit(1);
                }
                (*coordarray)[replicateindex].chainid=i+r*Nchainsinput;
 				(*coordarray)[replicateindex].leafid+=2*r;
                (*monomerid)[i+r*Nchainsinput][j].sidechain=replicateindex;
                (*coordarray)[replicateindex].r=add_double_triple((*coordarray)[replicateindex].r, scalar_multiply_double_triple(displacement, r));
                fmod_double_triple(&((*coordarray)[replicateindex].r), (*pbox_dimension));
                replicateindex++;
            }
        }
    }
    
    fscanf(inp, "%lli %i %i", &(*pt), &(*pframe), &(*prunningtime));
}

void output_coords(char *filename, int Nnodes, int Nchains, int *chainlength, monomernodes **monomerid, coord *coordarray, double_triple box_dimension, double vacuumthickness, int Nnodetypes, monomertypes monomercount, long long int t, int frame, int runningtime, double_triple *meshpositions, triangledata *meshtriangledata, bonddata *meshbonddata, vertexdata *meshdata, solvation_parameters solvationparameters){
	int i, j;
	FILE *outp;
	outp=fopen(filename, "w");
	fprintf(outp, "%.8f %.8f %.8f\n", box_dimension.x, box_dimension.y, box_dimension.z);
	fprintf(outp, "%i", solvationparameters.interface);
	if(solvationparameters.interface==1)fprintf(outp, " %.8f\n", vacuumthickness);
	else fprintf(outp, "\n");
	fprintf(outp, "%i %i %i %i %i %i", Nnodes, Nchains, monomercount.charged, monomercount.nonpolar, monomercount.monomers, Nnodetypes);
	for(i=0;i<Nnodetypes;i++) fprintf(outp, " %i", monomercount.type[i]);
	fprintf(outp, "\n");
	for(i=0;i<Nchains;i++){
		fprintf(outp, "%i", chainlength[i]);
		for(j=0;j<chainlength[i];j++){
			fprintf(outp, " %i %i", monomerid[i][j].backbone, monomerid[i][j].sidechain);
		}
		fprintf(outp, "\n");
	}
	for(i=0;i<Nnodes;i++){
		fprintf(outp, "%.8f %.8f %.8f %.8f %.8f %.8f %i %i %i %i %i\n", coordarray[i].r.x, coordarray[i].r.y, coordarray[i].r.z, coordarray[i].n.x, coordarray[i].n.y, coordarray[i].n.z, coordarray[i].nodetype, coordarray[i].chainid, coordarray[i].monomerid, coordarray[i].orientationtype, coordarray[i].leafid);
	}
    fprintf(outp, "%lli %i %i\n", t, frame, runningtime);
	if(solvationparameters.interface==1){
		fprintf(outp, "%i %i\n", solvationparameters.meshsize.x, solvationparameters.meshsize.y);
		fprintf(outp, "%lf\n", solvationparameters.meshbondlength);
		fprintf(outp, "%lf\n", solvationparameters.maxmeshbond);
		fprintf(outp, "%lf\n", solvationparameters.longestmeshbond);
		fprintf(outp, "%lf\n", solvationparameters.totalvolume);
		fprintf(outp, "%lf\n", solvationparameters.targetvolume);
		fprintf(outp, "%lf %lf %lf %lf %lf %lf\n", solvationparameters.totallength[0], solvationparameters.totallength[1], solvationparameters.totallength[2], solvationparameters.totallength[3], solvationparameters.totallength[4], solvationparameters.totallength[5]);
		for(i=0;i<solvationparameters.meshsize.y;i++){
			for(j=0;j<solvationparameters.meshsize.x;j++){
				fprintf(outp, "%.8f %.8f %.8f\n", meshpositions[i*solvationparameters.meshsize.x+j].x, meshpositions[i*solvationparameters.meshsize.x+j].y, meshpositions[i*solvationparameters.meshsize.x+j].z);
			}
		}
		for(i=0;i<solvationparameters.meshsize.y;i++){
			for(j=0;j<solvationparameters.meshsize.x;j++){
				fprintf(outp, "%.8f %.8f %.8f\n", meshpositions[(solvationparameters.meshsize.y+i)*solvationparameters.meshsize.x+j].x, meshpositions[(solvationparameters.meshsize.y+i)*solvationparameters.meshsize.x+j].y, meshpositions[(solvationparameters.meshsize.y+i)*solvationparameters.meshsize.x+j].z);
			}
		}
        for(i=0;i<2*solvationparameters.meshnumber;i++){
            fprintf(outp, "%.8f %.8f %.8f %.8f %.8f\n", meshtriangledata[i].area, meshtriangledata[i].columnvolume, meshtriangledata[i].normal.x, meshtriangledata[i].normal.y, meshtriangledata[i].normal.z);
        }
        for(i=0;i<3*solvationparameters.meshnumber;i++){
            fprintf(outp, "%.8f %.8f\n", meshbonddata[i].length, meshbonddata[i].meancurvature);
        }
        for(i=0;i<solvationparameters.meshnumber;i++){
            fprintf(outp, "%.8f %.8f %.8f\n", meshdata[i].area, meshdata[i].meancurvature, meshdata[i].meancurvaturesquared);
        }
	}
	
	fclose(outp);
}

void output_seed(char *base){
    FILE *outp;
    int i;
    declare_array(char, seedfile, 400);
    sprintf(seedfile, "%s/random_seed", base);
    outp=fopen(seedfile, "w");
    fprintf(outp, "%lu\n", getseed());
    fclose(outp);
    free(seedfile);
}

void input_seed(char *inputfile){
    int i;
    unsigned long inputseed;
    declare_array(char, directory, 400);
    declare_array(char, seedfile, 400);
    memcpy(directory, inputfile, strlen(inputfile) - 6);
    directory[strlen(inputfile) - 6] = '\0';
    sprintf(seedfile, "%s/random_seed", directory);
    FILE *inp;
    inp=fopen(seedfile, "r");
    if(inp){
        fscanf(inp, "%lu", &(inputseed));
        fclose(inp);
        mynextseed=inputseed;
        printf("mynextseed=%lu\n", mynextseed);
    }
    free(directory);
    free(seedfile);
}    

void configure_cells_struct(linkedlistfull *plink, double_triple box_dimension){
	(*plink).core.cellsperside.x=floor(box_dimension.x/(*plink).mincellwidth)+1;
	(*plink).core.cellsperside.y=floor(box_dimension.y/(*plink).mincellwidth)+1;
	(*plink).core.cellsperside.z=floor(box_dimension.z/(*plink).mincellwidth)+1;
	if(((*plink).core.cellsperside.x<3)||((*plink).core.cellsperside.y<3)||((*plink).core.cellsperside.y<3)){
		printf("code not written for small cellsperside (%i x %i x %i)\n", (*plink).core.cellsperside.x, (*plink).core.cellsperside.y, (*plink).core.cellsperside.z);
		printf("box_dimension=%f x %f x %f\n", box_dimension.x, box_dimension.y, box_dimension.z);
		exit(1);
	}
	(*plink).cellwidth.x=box_dimension.x/(*plink).core.cellsperside.x;
	(*plink).cellwidth.y=box_dimension.y/(*plink).core.cellsperside.y;
	(*plink).cellwidth.z=box_dimension.z/(*plink).core.cellsperside.z;
}

void configure_cells_struct_asymmetric(linkedlistasymmetric *plink, double_triple box_dimension){
	(*plink).core.cellsperside.x=floor(box_dimension.x/(*plink).mincellwidth)+1;
	(*plink).core.cellsperside.y=floor(box_dimension.y/(*plink).mincellwidth)+1;
	(*plink).core.cellsperside.z=floor(box_dimension.z/(*plink).mincellwidth)+1;
	if(((*plink).core.cellsperside.x<3)||((*plink).core.cellsperside.y<3)||((*plink).core.cellsperside.y<3)){
		printf("code not written for small cellsperside (%i x %i x %i)\n", (*plink).core.cellsperside.x, (*plink).core.cellsperside.y, (*plink).core.cellsperside.z);
		printf("box_dimension=%f x %f x %f\n", box_dimension.x, box_dimension.y, box_dimension.z);
		exit(1);
	}
	(*plink).cellwidth.x=box_dimension.x/(*plink).core.cellsperside.x;
	(*plink).cellwidth.y=box_dimension.y/(*plink).core.cellsperside.y;
	(*plink).cellwidth.z=box_dimension.z/(*plink).core.cellsperside.z;
}

void allocate_linklist(linkedlistfull *plink, int Nnodes){
	int i, j, k;
	allocate_3d_tensor(int, (*plink).core.head, (*plink).core.cellsperside.x, (*plink).core.cellsperside.y, (*plink).core.cellsperside.z);
	allocate_array(int, (*plink).core.list, Nnodes);
	allocate_array(int, (*plink).reverselist, Nnodes);
	allocate_array(int_triple, (*plink).core.cell, Nnodes);
	
	allocate_array(int, (*plink).core.number_neighbors, Nnodes);
	allocate_matrix(int, (*plink).core.neighbor, Nnodes, (*plink).maxneighbors);
}

void allocate_linklist_pairenergies(linkedlistfull *plink, int Nnodes){
	int i, j, k;
	allocate_3d_tensor(int, (*plink).core.head, (*plink).core.cellsperside.x, (*plink).core.cellsperside.y, (*plink).core.cellsperside.z);
	allocate_array(int, (*plink).core.list, Nnodes);
	allocate_array(int, (*plink).reverselist, Nnodes);
	allocate_array(int_triple, (*plink).core.cell, Nnodes);
	
	allocate_array(int, (*plink).core.number_neighbors, Nnodes);
	allocate_matrix(int, (*plink).core.neighbor, Nnodes, (*plink).maxneighbors);
	allocate_matrix(double, (*plink).core.pair_energy, Nnodes, (*plink).maxneighbors);
}

void free_linklist(linkedlistfull link, int Nnodes){
    int i, j;
    free_3d_tensor(link.core.head, link.core.cellsperside.x, link.core.cellsperside.y);
    free(link.core.list);
    free(link.reverselist);
    free(link.core.cell);
    free(link.core.number_neighbors);
    free_matrix(link.core.neighbor, Nnodes);
}

void free_linklist_pairenergies(linkedlistfull link, int Nnodes){
    int i, j;
    free_3d_tensor(link.core.head, link.core.cellsperside.x, link.core.cellsperside.y);
    free(link.core.list);
    free(link.reverselist);
    free(link.core.cell);
    free(link.core.number_neighbors);
    free_matrix(link.core.neighbor, Nnodes);
    free_matrix(link.core.pair_energy, Nnodes);
}

void cellconstructdoublylinked(coord *coordarray, int N, linkedlistfull *plink, int lownodetype, int highnodetype){
	int i, j, k, n;
	for(i=0;i<(*plink).core.cellsperside.x;i++){
		for(j=0;j<(*plink).core.cellsperside.y;j++){
			for(k=0;k<(*plink).core.cellsperside.z;k++){
				(*plink).core.head[i][j][k]=-1;						//	empty
			}
		}
	}
	for(n=0;n<N;n++) (*plink).reverselist[n]=-1;
	for(n=0;n<N;n++){
		if((coordarray[n].nodetype>=lownodetype)&&(coordarray[n].nodetype<=highnodetype)){
			i=(int) floor(coordarray[n].r.x/(*plink).cellwidth.x);
			j=(int) floor(coordarray[n].r.y/(*plink).cellwidth.y);
			k=(int) floor(coordarray[n].r.z/(*plink).cellwidth.z);
			i=mod(i, (*plink).core.cellsperside.x);
			j=mod(j, (*plink).core.cellsperside.y);
			k=mod(k, (*plink).core.cellsperside.z);
			(*plink).core.cell[n].x=i;
			(*plink).core.cell[n].y=j;
			(*plink).core.cell[n].z=k;
			(*plink).core.list[n]=(*plink).core.head[i][j][k];
			if((*plink).core.list[n]!=-1) (*plink).reverselist[(*plink).core.list[n]]=n;
			(*plink).core.head[i][j][k]=n;
		}
	}
}

void cellconstructdoublylinked_position(double_triple *position, int N, linkedlistfull *plink){
	int i, j, k, n;
	for(i=0;i<(*plink).core.cellsperside.x;i++){
		for(j=0;j<(*plink).core.cellsperside.y;j++){
			for(k=0;k<(*plink).core.cellsperside.z;k++){
				(*plink).core.head[i][j][k]=-1;						//	empty
			}
		}
	}
	for(n=0;n<N;n++) (*plink).reverselist[n]=-1;
	for(n=0;n<N;n++){
        i=(int) floor(position[n].x/(*plink).cellwidth.x);
        j=(int) floor(position[n].y/(*plink).cellwidth.y);
        k=(int) floor(position[n].z/(*plink).cellwidth.z);
        i=mod(i, (*plink).core.cellsperside.x);
        j=mod(j, (*plink).core.cellsperside.y);
        k=mod(k, (*plink).core.cellsperside.z);
        (*plink).core.cell[n].x=i;
        (*plink).core.cell[n].y=j;
        (*plink).core.cell[n].z=k;
        (*plink).core.list[n]=(*plink).core.head[i][j][k];
        if((*plink).core.list[n]!=-1) (*plink).reverselist[(*plink).core.list[n]]=n;
        (*plink).core.head[i][j][k]=n;
    }
}

void constructneighborlist(linkedlistfull *plink, int N, coord *coordarray, double_triple box_dimension){
	int i, j, k, self, neighbor, cellshifti, cellshiftj, cellshiftk, neighborcelli, neighborcellj, neighborcellk;
	double_triple selfpos;
	for(i=0;i<N;i++) (*plink).core.number_neighbors[i]=0;
	for(i=0;i<(*plink).core.cellsperside.x;i++){
		for(j=0;j<(*plink).core.cellsperside.y;j++){
			for(k=0;k<(*plink).core.cellsperside.z;k++){
				self=(*plink).core.head[i][j][k];
				while(self!=-1){
					selfpos=coordarray[self].r;
					neighbor=(*plink).core.list[self];										//	searching within own cell, only downstream in linked list
					while(neighbor!=-1){
						if(norm(subtract_double_triple(coordarray[neighbor].r, selfpos))<(*plink).mincellwidth){
							(*plink).core.neighbor[self][(*plink).core.number_neighbors[self]]=neighbor;
							(*plink).core.neighbor[neighbor][(*plink).core.number_neighbors[neighbor]]=self;
							(*plink).core.number_neighbors[self]++;
							(*plink).core.number_neighbors[neighbor]++;
						}
						neighbor=(*plink).core.list[neighbor];
					}
					
					cellshifti=1;														//	searching through 13 of 26 neighboring cells
					if(i==(*plink).core.cellsperside.x-1){
						neighborcelli=0;
						selfpos.x-=box_dimension.x;
					}
					else neighborcelli=i+cellshifti;
					for(cellshiftj=-1;cellshiftj<=1;cellshiftj++){
						if((j==0)&&(cellshiftj==-1)){
							neighborcellj=(*plink).core.cellsperside.y-1;
							selfpos.y+=box_dimension.y;
						}
						else if((j==(*plink).core.cellsperside.y-1)&&(cellshiftj==1)){
							neighborcellj=0;
							selfpos.y-=box_dimension.y;
						}
						else neighborcellj=j+cellshiftj;
						for(cellshiftk=-1;cellshiftk<=1;cellshiftk++){
							if((k==0)&&(cellshiftk==-1)){
								neighborcellk=(*plink).core.cellsperside.z-1;
								selfpos.z+=box_dimension.z;
							}
							else if((k==(*plink).core.cellsperside.z-1)&&(cellshiftk==1)){
								neighborcellk=0;
								selfpos.z-=box_dimension.z;
							}
							else neighborcellk=k+cellshiftk;
							neighbor=(*plink).core.head[neighborcelli][neighborcellj][neighborcellk];
							while(neighbor!=-1){
								if(norm(subtract_double_triple(coordarray[neighbor].r, selfpos))<(*plink).mincellwidth){
									(*plink).core.neighbor[self][(*plink).core.number_neighbors[self]]=neighbor;
									(*plink).core.neighbor[neighbor][(*plink).core.number_neighbors[neighbor]]=self;
									(*plink).core.number_neighbors[self]++;
									(*plink).core.number_neighbors[neighbor]++;
								}
								neighbor=(*plink).core.list[neighbor];
							}
							if((k==0)&&(cellshiftk==-1)) selfpos.z-=box_dimension.z;
							else if((k==(*plink).core.cellsperside.z-1)&&(cellshiftk==1)) selfpos.z+=box_dimension.z;
						}
						if((j==0)&&(cellshiftj==-1)) selfpos.y-=box_dimension.y;
						else if((j==(*plink).core.cellsperside.y-1)&&(cellshiftj==1)) selfpos.y+=box_dimension.y;
					}
					if(i==(*plink).core.cellsperside.x-1) selfpos.x+=box_dimension.x;
					
					neighborcelli=i;
					cellshiftj=1;
					if(j==(*plink).core.cellsperside.y-1){
						neighborcellj=0;
						selfpos.y-=box_dimension.y;
					}
					else neighborcellj=j+cellshiftj;
					for(cellshiftk=-1;cellshiftk<=1;cellshiftk++){
						if((k==0)&&(cellshiftk==-1)){
							neighborcellk=(*plink).core.cellsperside.z-1;
							selfpos.z+=box_dimension.z;
						}
						else if((k==(*plink).core.cellsperside.z-1)&&(cellshiftk==1)){
							neighborcellk=0;
							selfpos.z-=box_dimension.z;
						}
						else neighborcellk=k+cellshiftk;
						neighbor=(*plink).core.head[neighborcelli][neighborcellj][neighborcellk];
						while(neighbor!=-1){
							if(norm(subtract_double_triple(coordarray[neighbor].r, selfpos))<(*plink).mincellwidth){
								(*plink).core.neighbor[self][(*plink).core.number_neighbors[self]]=neighbor;
								(*plink).core.neighbor[neighbor][(*plink).core.number_neighbors[neighbor]]=self;
								(*plink).core.number_neighbors[self]++;
								(*plink).core.number_neighbors[neighbor]++;
							}
							neighbor=(*plink).core.list[neighbor];
						}
						if((k==0)&&(cellshiftk==-1)) selfpos.z-=box_dimension.z;
						else if((k==(*plink).core.cellsperside.z-1)&&(cellshiftk==1)) selfpos.z+=box_dimension.z;
					}
					if(j==(*plink).core.cellsperside.y-1) selfpos.y+=box_dimension.y;
                    
					neighborcelli=i;
					neighborcellj=j;
					cellshiftk=1;
					if(k==(*plink).core.cellsperside.z-1){
						neighborcellk=0;
						selfpos.z-=box_dimension.z;
					}
					else neighborcellk=k+cellshiftk;
					neighbor=(*plink).core.head[neighborcelli][neighborcellj][neighborcellk];
					while(neighbor!=-1){
						if(norm(subtract_double_triple(coordarray[neighbor].r, selfpos))<(*plink).mincellwidth){
							(*plink).core.neighbor[self][(*plink).core.number_neighbors[self]]=neighbor;
							(*plink).core.neighbor[neighbor][(*plink).core.number_neighbors[neighbor]]=self;
							(*plink).core.number_neighbors[self]++;
							(*plink).core.number_neighbors[neighbor]++;
						}
						neighbor=(*plink).core.list[neighbor];
					}
					if(k==(*plink).core.cellsperside.z-1) selfpos.z+=box_dimension.z;
					
					self=(*plink).core.list[self];
				}
			}
		}
	}
}

void constructneighborlist_positions(linkedlistfull *plink, int N, double_triple *positions, double_triple box_dimension){
	int i, j, k, self, neighbor, cellshifti, cellshiftj, cellshiftk, neighborcelli, neighborcellj, neighborcellk;
	double_triple selfpos;
	for(i=0;i<N;i++) (*plink).core.number_neighbors[i]=0;
	for(i=0;i<(*plink).core.cellsperside.x;i++){
		for(j=0;j<(*plink).core.cellsperside.y;j++){
			for(k=0;k<(*plink).core.cellsperside.z;k++){
				self=(*plink).core.head[i][j][k];
				while(self!=-1){
					selfpos=positions[self];
					neighbor=(*plink).core.list[self];										//	searching within own cell, only downstream in linked list
					while(neighbor!=-1){
						if(norm(subtract_double_triple(positions[neighbor], selfpos))<(*plink).mincellwidth){
							(*plink).core.neighbor[self][(*plink).core.number_neighbors[self]]=neighbor;
							(*plink).core.neighbor[neighbor][(*plink).core.number_neighbors[neighbor]]=self;
							(*plink).core.number_neighbors[self]++;
							(*plink).core.number_neighbors[neighbor]++;
                            if((*plink).core.number_neighbors[self]==(*plink).maxneighbors){
                                printf("number_neighbors[%i]=%i!\n", self, (*plink).maxneighbors);
                                exit(1);
                            }
                            if((*plink).core.number_neighbors[neighbor]==(*plink).maxneighbors){
                                printf("number_neighbors[%i]=%i!\n", self, (*plink).maxneighbors);
                                exit(1);
                            }
						}
						neighbor=(*plink).core.list[neighbor];
					}
					
					cellshifti=1;														//	searching through 13 of 26 neighboring cells
					if(i==(*plink).core.cellsperside.x-1){
						neighborcelli=0;
						selfpos.x-=box_dimension.x;
					}
					else neighborcelli=i+cellshifti;
					for(cellshiftj=-1;cellshiftj<=1;cellshiftj++){
						if((j==0)&&(cellshiftj==-1)){
							neighborcellj=(*plink).core.cellsperside.y-1;
							selfpos.y+=box_dimension.y;
						}
						else if((j==(*plink).core.cellsperside.y-1)&&(cellshiftj==1)){
							neighborcellj=0;
							selfpos.y-=box_dimension.y;
						}
						else neighborcellj=j+cellshiftj;
						for(cellshiftk=-1;cellshiftk<=1;cellshiftk++){
							if((k==0)&&(cellshiftk==-1)){
								neighborcellk=(*plink).core.cellsperside.z-1;
								selfpos.z+=box_dimension.z;
							}
							else if((k==(*plink).core.cellsperside.z-1)&&(cellshiftk==1)){
								neighborcellk=0;
								selfpos.z-=box_dimension.z;
							}
							else neighborcellk=k+cellshiftk;
							neighbor=(*plink).core.head[neighborcelli][neighborcellj][neighborcellk];
							while(neighbor!=-1){
								if(norm(subtract_double_triple(positions[neighbor], selfpos))<(*plink).mincellwidth){
									(*plink).core.neighbor[self][(*plink).core.number_neighbors[self]]=neighbor;
									(*plink).core.neighbor[neighbor][(*plink).core.number_neighbors[neighbor]]=self;
									(*plink).core.number_neighbors[self]++;
									(*plink).core.number_neighbors[neighbor]++;
                                    if((*plink).core.number_neighbors[self]==(*plink).maxneighbors){
                                        printf("number_neighbors[%i]=%i!\n", self, (*plink).maxneighbors);
                                        exit(1);
                                    }
                                    if((*plink).core.number_neighbors[neighbor]==(*plink).maxneighbors){
                                        printf("number_neighbors[%i]=%i!\n", self, (*plink).maxneighbors);
                                        exit(1);
                                    }
								}
								neighbor=(*plink).core.list[neighbor];
							}
							if((k==0)&&(cellshiftk==-1)) selfpos.z-=box_dimension.z;
							else if((k==(*plink).core.cellsperside.z-1)&&(cellshiftk==1)) selfpos.z+=box_dimension.z;
						}
						if((j==0)&&(cellshiftj==-1)) selfpos.y-=box_dimension.y;
						else if((j==(*plink).core.cellsperside.y-1)&&(cellshiftj==1)) selfpos.y+=box_dimension.y;
					}
					if(i==(*plink).core.cellsperside.x-1) selfpos.x+=box_dimension.x;
					
					neighborcelli=i;
					cellshiftj=1;
					if(j==(*plink).core.cellsperside.y-1){
						neighborcellj=0;
						selfpos.y-=box_dimension.y;
					}
					else neighborcellj=j+cellshiftj;
					for(cellshiftk=-1;cellshiftk<=1;cellshiftk++){
						if((k==0)&&(cellshiftk==-1)){
							neighborcellk=(*plink).core.cellsperside.z-1;
							selfpos.z+=box_dimension.z;
						}
						else if((k==(*plink).core.cellsperside.z-1)&&(cellshiftk==1)){
							neighborcellk=0;
							selfpos.z-=box_dimension.z;
						}
						else neighborcellk=k+cellshiftk;
						neighbor=(*plink).core.head[neighborcelli][neighborcellj][neighborcellk];
						while(neighbor!=-1){
							if(norm(subtract_double_triple(positions[neighbor], selfpos))<(*plink).mincellwidth){
								(*plink).core.neighbor[self][(*plink).core.number_neighbors[self]]=neighbor;
								(*plink).core.neighbor[neighbor][(*plink).core.number_neighbors[neighbor]]=self;
								(*plink).core.number_neighbors[self]++;
								(*plink).core.number_neighbors[neighbor]++;
                                if((*plink).core.number_neighbors[self]==(*plink).maxneighbors){
                                    printf("number_neighbors[%i]=%i!\n", self, (*plink).maxneighbors);
                                    exit(1);
                                }
                                if((*plink).core.number_neighbors[neighbor]==(*plink).maxneighbors){
                                    printf("number_neighbors[%i]=%i!\n", self, (*plink).maxneighbors);
                                    exit(1);
                                }
							}
							neighbor=(*plink).core.list[neighbor];
						}
						if((k==0)&&(cellshiftk==-1)) selfpos.z-=box_dimension.z;
						else if((k==(*plink).core.cellsperside.z-1)&&(cellshiftk==1)) selfpos.z+=box_dimension.z;
					}
					if(j==(*plink).core.cellsperside.y-1) selfpos.y+=box_dimension.y;
                    
					neighborcelli=i;
					neighborcellj=j;
					cellshiftk=1;
					if(k==(*plink).core.cellsperside.z-1){
						neighborcellk=0;
						selfpos.z-=box_dimension.z;
					}
					else neighborcellk=k+cellshiftk;
					neighbor=(*plink).core.head[neighborcelli][neighborcellj][neighborcellk];
					while(neighbor!=-1){
						if(norm(subtract_double_triple(positions[neighbor], selfpos))<(*plink).mincellwidth){
							(*plink).core.neighbor[self][(*plink).core.number_neighbors[self]]=neighbor;
							(*plink).core.neighbor[neighbor][(*plink).core.number_neighbors[neighbor]]=self;
							(*plink).core.number_neighbors[self]++;
							(*plink).core.number_neighbors[neighbor]++;
                            if((*plink).core.number_neighbors[self]==(*plink).maxneighbors){
                                printf("number_neighbors[%i]=%i!\n", self, (*plink).maxneighbors);
                                exit(1);
                            }
                            if((*plink).core.number_neighbors[neighbor]==(*plink).maxneighbors){
                                printf("number_neighbors[%i]=%i!\n", self, (*plink).maxneighbors);
                                exit(1);
                            }
						}
						neighbor=(*plink).core.list[neighbor];
					}
					if(k==(*plink).core.cellsperside.z-1) selfpos.z+=box_dimension.z;
					
					self=(*plink).core.list[self];
				}
			}
		}
	}
}

void constructneighborlist_pairenergies_phenyl(linkedlistfull *plink, int N, coord *coordarray, double_triple box_dimension, nonbonded_params my_nonbonded_params){
	int i, j, k, self, neighbor, cellshifti, cellshiftj, cellshiftk, neighborcelli, neighborcellj, neighborcellk;
	double energy;
	double_triple selfpos;
	for(i=0;i<N;i++) (*plink).core.number_neighbors[i]=0;
	for(i=0;i<(*plink).core.cellsperside.x;i++){
		for(j=0;j<(*plink).core.cellsperside.y;j++){
			for(k=0;k<(*plink).core.cellsperside.z;k++){
				self=(*plink).core.head[i][j][k];
				while(self!=-1){
					selfpos=coordarray[self].r;
					neighbor=(*plink).core.list[self];										//	searching within own cell, only downstream in linked list
					while(neighbor!=-1){
						if(norm(subtract_double_triple(coordarray[neighbor].r, selfpos))<(*plink).mincellwidth){
							(*plink).core.neighbor[self][(*plink).core.number_neighbors[self]]=neighbor;
							(*plink).core.neighbor[neighbor][(*plink).core.number_neighbors[neighbor]]=self;
							
							energy=phenylphenyl_energy(coordarray[self], coordarray[neighbor], my_nonbonded_params.p11vac, box_dimension, my_nonbonded_params.phenylcutoff2, my_nonbonded_params.solvationparams);
							
							(*plink).core.pair_energy[self][(*plink).core.number_neighbors[self]]=energy;
							(*plink).core.pair_energy[neighbor][(*plink).core.number_neighbors[neighbor]]=energy;
							(*plink).core.number_neighbors[self]++;
							(*plink).core.number_neighbors[neighbor]++;
						}
						neighbor=(*plink).core.list[neighbor];
					}
					
					cellshifti=1;														//	searching through 13 of 26 neighboring cells
					if(i==(*plink).core.cellsperside.x-1){
						neighborcelli=0;
						selfpos.x-=box_dimension.x;
					}
					else neighborcelli=i+cellshifti;
					for(cellshiftj=-1;cellshiftj<=1;cellshiftj++){
						if((j==0)&&(cellshiftj==-1)){
							neighborcellj=(*plink).core.cellsperside.y-1;
							selfpos.y+=box_dimension.y;
						}
						else if((j==(*plink).core.cellsperside.y-1)&&(cellshiftj==1)){
							neighborcellj=0;
							selfpos.y-=box_dimension.y;
						}
						else neighborcellj=j+cellshiftj;
						for(cellshiftk=-1;cellshiftk<=1;cellshiftk++){
							if((k==0)&&(cellshiftk==-1)){
								neighborcellk=(*plink).core.cellsperside.z-1;
								selfpos.z+=box_dimension.z;
							}
							else if((k==(*plink).core.cellsperside.z-1)&&(cellshiftk==1)){
								neighborcellk=0;
								selfpos.z-=box_dimension.z;
							}
							else neighborcellk=k+cellshiftk;
							neighbor=(*plink).core.head[neighborcelli][neighborcellj][neighborcellk];
							while(neighbor!=-1){
								if(norm(subtract_double_triple(coordarray[neighbor].r, selfpos))<(*plink).mincellwidth){
									(*plink).core.neighbor[self][(*plink).core.number_neighbors[self]]=neighbor;
									(*plink).core.neighbor[neighbor][(*plink).core.number_neighbors[neighbor]]=self;
									//energy=phenylphenyl_energy(coordarray[self], coordarray[neighbor], my_nonbonded_params.p11vac, box_dimension, my_nonbonded_params.phenylcutoff2);
                                    
									energy=phenylphenyl_energy(coordarray[self], coordarray[neighbor], my_nonbonded_params.p11vac, box_dimension, my_nonbonded_params.phenylcutoff2, my_nonbonded_params.solvationparams);
									
									(*plink).core.pair_energy[self][(*plink).core.number_neighbors[self]]=energy;
									(*plink).core.pair_energy[neighbor][(*plink).core.number_neighbors[neighbor]]=energy;
									(*plink).core.number_neighbors[self]++;
									(*plink).core.number_neighbors[neighbor]++;
								}
								neighbor=(*plink).core.list[neighbor];
							}
							if((k==0)&&(cellshiftk==-1)) selfpos.z-=box_dimension.z;
							else if((k==(*plink).core.cellsperside.z-1)&&(cellshiftk==1)) selfpos.z+=box_dimension.z;
						}
						if((j==0)&&(cellshiftj==-1)) selfpos.y-=box_dimension.y;
						else if((j==(*plink).core.cellsperside.y-1)&&(cellshiftj==1)) selfpos.y+=box_dimension.y;
					}
					if(i==(*plink).core.cellsperside.x-1) selfpos.x+=box_dimension.x;
					
					neighborcelli=i;
					cellshiftj=1;
					if(j==(*plink).core.cellsperside.y-1){
						neighborcellj=0;
						selfpos.y-=box_dimension.y;
					}
					else neighborcellj=j+cellshiftj;
					for(cellshiftk=-1;cellshiftk<=1;cellshiftk++){
						if((k==0)&&(cellshiftk==-1)){
							neighborcellk=(*plink).core.cellsperside.z-1;
							selfpos.z+=box_dimension.z;
						}
						else if((k==(*plink).core.cellsperside.z-1)&&(cellshiftk==1)){
							neighborcellk=0;
							selfpos.z-=box_dimension.z;
						}
						else neighborcellk=k+cellshiftk;
						neighbor=(*plink).core.head[neighborcelli][neighborcellj][neighborcellk];
						while(neighbor!=-1){
							if(norm(subtract_double_triple(coordarray[neighbor].r, selfpos))<(*plink).mincellwidth){
								(*plink).core.neighbor[self][(*plink).core.number_neighbors[self]]=neighbor;
								(*plink).core.neighbor[neighbor][(*plink).core.number_neighbors[neighbor]]=self;
                                
								energy=phenylphenyl_energy(coordarray[self], coordarray[neighbor], my_nonbonded_params.p11vac, box_dimension, my_nonbonded_params.phenylcutoff2, my_nonbonded_params.solvationparams);
								
								(*plink).core.pair_energy[self][(*plink).core.number_neighbors[self]]=energy;
								(*plink).core.pair_energy[neighbor][(*plink).core.number_neighbors[neighbor]]=energy;
								(*plink).core.number_neighbors[self]++;
								(*plink).core.number_neighbors[neighbor]++;
							}
							neighbor=(*plink).core.list[neighbor];
						}
						if((k==0)&&(cellshiftk==-1)) selfpos.z-=box_dimension.z;
						else if((k==(*plink).core.cellsperside.z-1)&&(cellshiftk==1)) selfpos.z+=box_dimension.z;
					}
					if(j==(*plink).core.cellsperside.y-1) selfpos.y+=box_dimension.y;
					
					neighborcelli=i;
					neighborcellj=j;
					cellshiftk=1;
					if(k==(*plink).core.cellsperside.z-1){
						neighborcellk=0;
						selfpos.z-=box_dimension.z;
					}
					else neighborcellk=k+cellshiftk;
					neighbor=(*plink).core.head[neighborcelli][neighborcellj][neighborcellk];
					while(neighbor!=-1){
						if(norm(subtract_double_triple(coordarray[neighbor].r, selfpos))<(*plink).mincellwidth){
							(*plink).core.neighbor[self][(*plink).core.number_neighbors[self]]=neighbor;
							(*plink).core.neighbor[neighbor][(*plink).core.number_neighbors[neighbor]]=self;
							
							energy=phenylphenyl_energy(coordarray[self], coordarray[neighbor], my_nonbonded_params.p11vac, box_dimension, my_nonbonded_params.phenylcutoff2, my_nonbonded_params.solvationparams);
							
							(*plink).core.pair_energy[self][(*plink).core.number_neighbors[self]]=energy;
							(*plink).core.pair_energy[neighbor][(*plink).core.number_neighbors[neighbor]]=energy;
							(*plink).core.number_neighbors[self]++;
							(*plink).core.number_neighbors[neighbor]++;
						}
						neighbor=(*plink).core.list[neighbor];
					}
					if(k==(*plink).core.cellsperside.z-1) selfpos.z+=box_dimension.z;
					
					self=(*plink).core.list[self];
				}
			}
		}
	}
}

void constructneighborlist_pairenergies_charged(linkedlistfull *plink, int N, coord *coordarray, double_triple box_dimension, nonbonded_params my_nonbonded_params){
	int i, j, k, self, neighbor, cellshifti, cellshiftj, cellshiftk, neighborcelli, neighborcellj, neighborcellk;
	double energy;
	double_triple selfpos;
	for(i=0;i<N;i++) (*plink).core.number_neighbors[i]=0;
	for(i=0;i<(*plink).core.cellsperside.x;i++){
		for(j=0;j<(*plink).core.cellsperside.y;j++){
			for(k=0;k<(*plink).core.cellsperside.z;k++){
				self=(*plink).core.head[i][j][k];
				while(self!=-1){
					selfpos=coordarray[self].r;
					neighbor=(*plink).core.list[self];										//	searching within own cell, only downstream in linked list
					while(neighbor!=-1){
						if(norm(subtract_double_triple(coordarray[neighbor].r, selfpos))<(*plink).mincellwidth){
							(*plink).core.neighbor[self][(*plink).core.number_neighbors[self]]=neighbor;
							(*plink).core.neighbor[neighbor][(*plink).core.number_neighbors[neighbor]]=self;
							if(coordarray[self].nodetype==2){
								if(coordarray[neighbor].nodetype==2){
									energy=electrostatic_energy(coordarray[self], coordarray[neighbor], my_nonbonded_params.solvationparams, my_nonbonded_params.p2, my_nonbonded_params.p2, box_dimension, my_nonbonded_params.cutoff2);
								}
								else if(coordarray[neighbor].nodetype==3){
									energy=electrostatic_energy(coordarray[self], coordarray[neighbor], my_nonbonded_params.solvationparams, my_nonbonded_params.p2, my_nonbonded_params.p3, box_dimension, my_nonbonded_params.cutoff2);
								}
							}
							else if(coordarray[self].nodetype==3){
								if(coordarray[neighbor].nodetype==2){
									energy=electrostatic_energy(coordarray[self], coordarray[neighbor], my_nonbonded_params.solvationparams, my_nonbonded_params.p3, my_nonbonded_params.p2, box_dimension, my_nonbonded_params.cutoff2);
								}
								else if(coordarray[neighbor].nodetype==3){
									energy=electrostatic_energy(coordarray[self], coordarray[neighbor], my_nonbonded_params.solvationparams, my_nonbonded_params.p3, my_nonbonded_params.p3, box_dimension, my_nonbonded_params.cutoff2);
								}
							}
							(*plink).core.pair_energy[self][(*plink).core.number_neighbors[self]]=energy;
							(*plink).core.pair_energy[neighbor][(*plink).core.number_neighbors[neighbor]]=energy;
							(*plink).core.number_neighbors[self]++;
							(*plink).core.number_neighbors[neighbor]++;
						}
						neighbor=(*plink).core.list[neighbor];
					}
					
					cellshifti=1;														//	searching through 13 of 26 neighboring cells
					if(i==(*plink).core.cellsperside.x-1){
						neighborcelli=0;
						selfpos.x-=box_dimension.x;
					}
					else neighborcelli=i+cellshifti;
					for(cellshiftj=-1;cellshiftj<=1;cellshiftj++){
						if((j==0)&&(cellshiftj==-1)){
							neighborcellj=(*plink).core.cellsperside.y-1;
							selfpos.y+=box_dimension.y;
						}
						else if((j==(*plink).core.cellsperside.y-1)&&(cellshiftj==1)){
							neighborcellj=0;
							selfpos.y-=box_dimension.y;
						}
						else neighborcellj=j+cellshiftj;
						for(cellshiftk=-1;cellshiftk<=1;cellshiftk++){
							if((k==0)&&(cellshiftk==-1)){
								neighborcellk=(*plink).core.cellsperside.z-1;
								selfpos.z+=box_dimension.z;
							}
							else if((k==(*plink).core.cellsperside.z-1)&&(cellshiftk==1)){
								neighborcellk=0;
								selfpos.z-=box_dimension.z;
							}
							else neighborcellk=k+cellshiftk;
							neighbor=(*plink).core.head[neighborcelli][neighborcellj][neighborcellk];
							while(neighbor!=-1){
								if(norm(subtract_double_triple(coordarray[neighbor].r, selfpos))<(*plink).mincellwidth){
									(*plink).core.neighbor[self][(*plink).core.number_neighbors[self]]=neighbor;
									(*plink).core.neighbor[neighbor][(*plink).core.number_neighbors[neighbor]]=self;
									if(coordarray[self].nodetype==2){
										if(coordarray[neighbor].nodetype==2){
											energy=electrostatic_energy(coordarray[self], coordarray[neighbor], my_nonbonded_params.solvationparams, my_nonbonded_params.p2, my_nonbonded_params.p2, box_dimension, my_nonbonded_params.cutoff2);
										}
										else if(coordarray[neighbor].nodetype==3){
											energy=electrostatic_energy(coordarray[self], coordarray[neighbor], my_nonbonded_params.solvationparams, my_nonbonded_params.p2, my_nonbonded_params.p3, box_dimension, my_nonbonded_params.cutoff2);
										}
									}
									else if(coordarray[self].nodetype==3){
										if(coordarray[neighbor].nodetype==2){
											energy=electrostatic_energy(coordarray[self], coordarray[neighbor], my_nonbonded_params.solvationparams, my_nonbonded_params.p3, my_nonbonded_params.p2, box_dimension, my_nonbonded_params.cutoff2);
										}
										else if(coordarray[neighbor].nodetype==3){
											energy=electrostatic_energy(coordarray[self], coordarray[neighbor], my_nonbonded_params.solvationparams, my_nonbonded_params.p3, my_nonbonded_params.p3, box_dimension, my_nonbonded_params.cutoff2);
										}
									}
									(*plink).core.pair_energy[self][(*plink).core.number_neighbors[self]]=energy;
									(*plink).core.pair_energy[neighbor][(*plink).core.number_neighbors[neighbor]]=energy;
									(*plink).core.number_neighbors[self]++;
									(*plink).core.number_neighbors[neighbor]++;
								}
								neighbor=(*plink).core.list[neighbor];
							}
							if((k==0)&&(cellshiftk==-1)) selfpos.z-=box_dimension.z;
							else if((k==(*plink).core.cellsperside.z-1)&&(cellshiftk==1)) selfpos.z+=box_dimension.z;
						}
						if((j==0)&&(cellshiftj==-1)) selfpos.y-=box_dimension.y;
						else if((j==(*plink).core.cellsperside.y-1)&&(cellshiftj==1)) selfpos.y+=box_dimension.y;
					}
					if(i==(*plink).core.cellsperside.x-1) selfpos.x+=box_dimension.x;
					
					neighborcelli=i;
					cellshiftj=1;
					if(j==(*plink).core.cellsperside.y-1){
						neighborcellj=0;
						selfpos.y-=box_dimension.y;
					}
					else neighborcellj=j+cellshiftj;
					for(cellshiftk=-1;cellshiftk<=1;cellshiftk++){
						if((k==0)&&(cellshiftk==-1)){
							neighborcellk=(*plink).core.cellsperside.z-1;
							selfpos.z+=box_dimension.z;
						}
						else if((k==(*plink).core.cellsperside.z-1)&&(cellshiftk==1)){
							neighborcellk=0;
							selfpos.z-=box_dimension.z;
						}
						else neighborcellk=k+cellshiftk;
						neighbor=(*plink).core.head[neighborcelli][neighborcellj][neighborcellk];
						while(neighbor!=-1){
							if(norm(subtract_double_triple(coordarray[neighbor].r, selfpos))<(*plink).mincellwidth){
								(*plink).core.neighbor[self][(*plink).core.number_neighbors[self]]=neighbor;
								(*plink).core.neighbor[neighbor][(*plink).core.number_neighbors[neighbor]]=self;
								if(coordarray[self].nodetype==2){
									if(coordarray[neighbor].nodetype==2){
										energy=electrostatic_energy(coordarray[self], coordarray[neighbor], my_nonbonded_params.solvationparams, my_nonbonded_params.p2, my_nonbonded_params.p2, box_dimension, my_nonbonded_params.cutoff2);
									}
									else if(coordarray[neighbor].nodetype==3){
										energy=electrostatic_energy(coordarray[self], coordarray[neighbor], my_nonbonded_params.solvationparams, my_nonbonded_params.p2, my_nonbonded_params.p3, box_dimension, my_nonbonded_params.cutoff2);
									}
								}
								else if(coordarray[self].nodetype==3){
									if(coordarray[neighbor].nodetype==2){
										energy=electrostatic_energy(coordarray[self], coordarray[neighbor], my_nonbonded_params.solvationparams, my_nonbonded_params.p3, my_nonbonded_params.p2, box_dimension, my_nonbonded_params.cutoff2);
									}
									else if(coordarray[neighbor].nodetype==3){
										energy=electrostatic_energy(coordarray[self], coordarray[neighbor], my_nonbonded_params.solvationparams, my_nonbonded_params.p3, my_nonbonded_params.p3, box_dimension, my_nonbonded_params.cutoff2);
									}
								}
								(*plink).core.pair_energy[self][(*plink).core.number_neighbors[self]]=energy;
								(*plink).core.pair_energy[neighbor][(*plink).core.number_neighbors[neighbor]]=energy;
								(*plink).core.number_neighbors[self]++;
								(*plink).core.number_neighbors[neighbor]++;
							}
							neighbor=(*plink).core.list[neighbor];
						}
						if((k==0)&&(cellshiftk==-1)) selfpos.z-=box_dimension.z;
						else if((k==(*plink).core.cellsperside.z-1)&&(cellshiftk==1)) selfpos.z+=box_dimension.z;
					}
					if(j==(*plink).core.cellsperside.y-1) selfpos.y+=box_dimension.y;
					
					neighborcelli=i;
					neighborcellj=j;
					cellshiftk=1;
					if(k==(*plink).core.cellsperside.z-1){
						neighborcellk=0;
						selfpos.z-=box_dimension.z;
					}
					else neighborcellk=k+cellshiftk;
					neighbor=(*plink).core.head[neighborcelli][neighborcellj][neighborcellk];
					while(neighbor!=-1){
						if(norm(subtract_double_triple(coordarray[neighbor].r, selfpos))<(*plink).mincellwidth){
							(*plink).core.neighbor[self][(*plink).core.number_neighbors[self]]=neighbor;
							(*plink).core.neighbor[neighbor][(*plink).core.number_neighbors[neighbor]]=self;
							if(coordarray[self].nodetype==2){
								if(coordarray[neighbor].nodetype==2){
									energy=electrostatic_energy(coordarray[self], coordarray[neighbor], my_nonbonded_params.solvationparams, my_nonbonded_params.p2, my_nonbonded_params.p2, box_dimension, my_nonbonded_params.cutoff2);
								}
								else if(coordarray[neighbor].nodetype==3){
									energy=electrostatic_energy(coordarray[self], coordarray[neighbor], my_nonbonded_params.solvationparams, my_nonbonded_params.p2, my_nonbonded_params.p3, box_dimension, my_nonbonded_params.cutoff2);
								}
							}
							else if(coordarray[self].nodetype==3){
								if(coordarray[neighbor].nodetype==2){
									energy=electrostatic_energy(coordarray[self], coordarray[neighbor], my_nonbonded_params.solvationparams, my_nonbonded_params.p3, my_nonbonded_params.p2, box_dimension, my_nonbonded_params.cutoff2);
								}
								else if(coordarray[neighbor].nodetype==3){
									energy=electrostatic_energy(coordarray[self], coordarray[neighbor], my_nonbonded_params.solvationparams, my_nonbonded_params.p3, my_nonbonded_params.p3, box_dimension, my_nonbonded_params.cutoff2);
								}
							}
							(*plink).core.pair_energy[self][(*plink).core.number_neighbors[self]]=energy;
							(*plink).core.pair_energy[neighbor][(*plink).core.number_neighbors[neighbor]]=energy;
							(*plink).core.number_neighbors[self]++;
							(*plink).core.number_neighbors[neighbor]++;
						}
						neighbor=(*plink).core.list[neighbor];
					}
					if(k==(*plink).core.cellsperside.z-1) selfpos.z+=box_dimension.z;
					self=(*plink).core.list[self];
				}
			}
		}
	}
}

void update_neighbor_list(int Nchains, int *chainlength, coord *coordarray, nonbonded_params my_nonbonded_params, linkedlistset linkset, monomernodes **monomerid, int *number_neighbors, int **neighbor, int max_neighbors, int *assigned){
	int i, j, index, target, targetpolymer, firstnode, lastnode;
	double rsq;
	double_triple r;
	for(index=0;index<Nchains;index++) number_neighbors[index]=0;
	for(index=0;index<Nchains;index++){
		for(i=index+1;i<Nchains;i++) assigned[i]=0;
		assigned[index]=1;
		firstnode=monomerid[index][0].backbone;
		lastnode=monomerid[index][chainlength[index]-1].sidechain;
		for(i=firstnode;i<=lastnode;i++){
			
			//	neighbor if interacting (but don't need to check hard core interactions, since these infinite interactions would not have been allowed)
			
			if(coordarray[i].nodetype==1){
				for(j=0;j<linkset.phenylphenyl.core.number_neighbors[i];j++){
					target=linkset.phenylphenyl.core.neighbor[i][j];
					targetpolymer=coordarray[target].chainid;
					if(targetpolymer>index){							//	to avoid double counting
						if(assigned[targetpolymer]==0){
							r=subtract_double_triple(coordarray[target].r, coordarray[i].r);
							rsq=dot_product(r, r);
							if(rsq<my_nonbonded_params.phenylcutoff2){				//	interacting
								assigned[targetpolymer]=1;
								neighbor[index][number_neighbors[index]]=targetpolymer;
								number_neighbors[index]++;
								if(number_neighbors[index]==max_neighbors) my_exit("max_neighbors not big enough");
								neighbor[targetpolymer][number_neighbors[targetpolymer]]=index;
								number_neighbors[targetpolymer]++;
								if(number_neighbors[targetpolymer]==max_neighbors) my_exit("max_neighbors not big enough");
							}
						}
					}
				}
			}
			else if((coordarray[i].nodetype==2)||(coordarray[i].nodetype==3)){
				for(j=0;j<linkset.chargedcharged.core.number_neighbors[i];j++){
					target=linkset.chargedcharged.core.neighbor[i][j];
					targetpolymer=coordarray[target].chainid;
					if(targetpolymer>index){							//	to avoid double counting
						if(assigned[targetpolymer]==0){
							r=subtract_double_triple(coordarray[target].r, coordarray[i].r);
							rsq=dot_product(r, r);
							if(rsq<my_nonbonded_params.cutoff2){				//	interacting
								assigned[targetpolymer]=1;
								neighbor[index][number_neighbors[index]]=targetpolymer;
								number_neighbors[index]++;
								if(number_neighbors[index]==max_neighbors) my_exit("max_neighbors not big enough");
								neighbor[targetpolymer][number_neighbors[targetpolymer]]=index;
 								number_neighbors[targetpolymer]++;
								if(number_neighbors[targetpolymer]==max_neighbors) my_exit("max_neighbors not big enough");
							}
						}
					}
				}
			}
		}
	}
}

void deletefile(char *filename){
	FILE *outp;
	outp=fopen(filename, "w");
	fclose(outp);
}

void initialize_energycomponents_nofiles(energycomponents *pmy_energycomponents, int Nnodetypes){
	(*pmy_energycomponents).backbone=0;
	(*pmy_energycomponents).sidechain=0;
	
	(*pmy_energycomponents).nn=0;
	(*pmy_energycomponents).ccunlike=0;
	(*pmy_energycomponents).cclike=0;
	(*pmy_energycomponents).cn=0;
    
	(*pmy_energycomponents).nncross=0;
	(*pmy_energycomponents).ccunlikecross=0;
	(*pmy_energycomponents).cclikecross=0;
	(*pmy_energycomponents).cncross=0;
    
	(*pmy_energycomponents).nndifferent=0;
	(*pmy_energycomponents).ccunlikedifferent=0;
	(*pmy_energycomponents).cclikedifferent=0;
	(*pmy_energycomponents).cndifferent=0;
    
	(*pmy_energycomponents).nnsamepoly=0;
	(*pmy_energycomponents).ccunlikesamepoly=0;
	(*pmy_energycomponents).cclikesamepoly=0;
	(*pmy_energycomponents).cnsamepoly=0;
	int i;
	for(i=0;i<Nnodetypes;i++){
		(*pmy_energycomponents).solvation[i]=0;
		(*pmy_energycomponents).interface[i]=0;
	}
}

void output_xyz_cg(int Nchains, int *chainlength, monomernodes **monomerid, int Natoms, coord *coordarray, double_triple box_dimension, char *xyzname, char *sourcename, cgparams params, int *pframe){
	int i, j, jmid, backboneindex, sidechainindex, atomcount=0, lastbackbonecount, lastbackboneindex, firstbackbonecount, newcgsourcefile=0;
	double_triple com, backbonebonded, sidevector, sidechainbonded, lastbackbonebonded, firstbackbonebonded, sep, lastbackbone;
	FILE *outp, *sourceoutp;
	
	if((*pframe)==0){
        newcgsourcefile=1;
		sourceoutp=fopen(sourcename, "w");
		fprintf(sourceoutp, "pbc set {%.6f %.6f %.6f} -first %i -last %i\n", box_dimension.x, box_dimension.y, box_dimension.z, *pframe, *pframe);
		fprintf(sourceoutp, "topo clearbonds\n");
		fprintf(sourceoutp, "set sel [atomselect top \" name C\"]\n$sel set radius %.6f\n", params.backrad);
		fprintf(sourceoutp, "set sel [atomselect top \" name H\"]\n$sel set radius %.6f\n", params.phenrad);
		fprintf(sourceoutp, "set sel [atomselect top \" name N\"]\n$sel set radius %.6f\n", params.aminrad);
		fprintf(sourceoutp, "set sel [atomselect top \" name O\"]\n$sel set radius %.6f\n", params.carbrad);
		fprintf(sourceoutp, "color Name C white\n");
		fprintf(sourceoutp, "color Name H yellow\n");
		fprintf(sourceoutp, "color Name N blue\n");
		fprintf(sourceoutp, "color Name O red\n");
		outp=fopen(xyzname, "w");
	}
	else{
        
		if ((sourceoutp = fopen(sourcename, "r")) == NULL) {
            newcgsourcefile=1;
			sourceoutp=fopen(sourcename, "w");
			fprintf(sourceoutp, "topo clearbonds\n");
			fprintf(sourceoutp, "set sel [atomselect top \" name C\"]\n$sel set radius %.6f\n", params.backrad);
			fprintf(sourceoutp, "set sel [atomselect top \" name H\"]\n$sel set radius %.6f\n", params.phenrad);
			fprintf(sourceoutp, "set sel [atomselect top \" name N\"]\n$sel set radius %.6f\n", params.aminrad);
			fprintf(sourceoutp, "set sel [atomselect top \" name O\"]\n$sel set radius %.6f\n", params.carbrad);
			
		} else {
			fclose(sourceoutp);
			sourceoutp=fopen(sourcename, "a");
		}
		fprintf(sourceoutp, "pbc set {%.6f %.6f %.6f} -first %i -last %i\n", box_dimension.x, box_dimension.y, box_dimension.z, *pframe, *pframe);
		if(newcgsourcefile==0) fclose(sourceoutp);
		outp=fopen(xyzname, "a");
	}
	fprintf(outp, "%i\n%i\t%.6f\t%.6f\t%.6f\n", Natoms, Nchains, box_dimension.x, box_dimension.y, box_dimension.z);
	for(i=0;i<Nchains;i++){
		com.x=com.y=com.z=0;
		for(j=0;j<chainlength[i];j++){
			if(j==0){
				com=coordarray[monomerid[i][j].backbone].r;
				lastbackbone=coordarray[monomerid[i][j].backbone].r;
			}
			else{
				sep=subtract_double_triple(coordarray[monomerid[i][j].backbone].r, lastbackbone);
				recenter_double_triple(&sep, box_dimension);
				lastbackbone=add_double_triple(lastbackbone, sep);
				com=add_double_triple(com, lastbackbone);
			}
		}
		com=scalar_multiply_double_triple(com, (1./chainlength[i]));
		fmod_double_triple(&com, box_dimension);
		for(jmid=0;jmid<chainlength[i];jmid++){
			if(jmid+chainlength[i]/2<chainlength[i]) j=chainlength[i]/2+jmid;
			else j=chainlength[i]/2-(1+jmid-(chainlength[i]-chainlength[i]/2));
			backboneindex=monomerid[i][j].backbone;
			if(j==chainlength[i]/2) firstbackbonebonded=backbonebonded=nearest_image(coordarray[backboneindex].r, com, box_dimension);
			else if(j==chainlength[i]/2-1) backbonebonded=nearest_image(coordarray[backboneindex].r, firstbackbonebonded, box_dimension);
			else backbonebonded=nearest_image(coordarray[backboneindex].r, lastbackbonebonded, box_dimension);
			lastbackbonebonded=backbonebonded;
			lastbackboneindex=backboneindex;
			outputcoords_centered(outp, "C", backbonebonded, box_dimension);				//	backbone N
            if(newcgsourcefile==1){
				if(jmid==0) firstbackbonecount=atomcount;
				else if(jmid==chainlength[i]-chainlength[i]/2){
					fprintf(sourceoutp, "topo addbond %i %i\n", firstbackbonecount, atomcount);
				}
				else{
					fprintf(sourceoutp, "topo addbond %i %i\n", lastbackbonecount, atomcount);
				}
				lastbackbonecount=atomcount;
				atomcount++;
			}
			sidechainindex=monomerid[i][j].sidechain;
            sidevector=subtract_double_triple(coordarray[sidechainindex].r, coordarray[backboneindex].r);
            recenter_double_triple(&sidevector, box_dimension);
            sidechainbonded=add_double_triple(backbonebonded, sidevector);
            if(coordarray[sidechainindex].nodetype>0){
                if(coordarray[sidechainindex].nodetype==1) outputcoords_centered(outp, "H", sidechainbonded, box_dimension);
                else if(coordarray[sidechainindex].nodetype==2) outputcoords_centered(outp, "N", sidechainbonded, box_dimension);
                else if(coordarray[sidechainindex].nodetype==3) outputcoords_centered(outp, "O", sidechainbonded, box_dimension);
                else my_exit("haven't written output_xyz_cg for other nodetypes");
                if(newcgsourcefile==1){
                    fprintf(sourceoutp, "topo addbond %i %i\n", atomcount-1, atomcount);
                    atomcount+=1;
                }
			}
		}
	}
    if(newcgsourcefile==1){
		fclose(sourceoutp);
	}
	fclose(outp);
	(*pframe)++;
}

int checktriangleoverlap(double_triple p1, double_triple q1, double_triple r1, double_triple p2, double_triple q2, double_triple r2, double_triple norm1, double_triple norm2, double_triple box_dimension){
	
	// using the Moller-Trumbore intersection algorithm because it makes sense to already have normals calculated
	
	int i, j;
	declare_matrix_nozero(double_triple, vertices, 2, 3);
	declare_matrix(double, verticesplaneframe, 2, 3);
	declare_array(double, planeoffset, 2);
	declare_array(int, singleindex, 2);
	declare_array(int, onplanecode, 2);
	
	vertices[0][0]=p1;
	vertices[0][1]=nearest_image(q1, p1, box_dimension);
	vertices[0][2]=nearest_image(r1, p1, box_dimension);
	vertices[1][0]=nearest_image(p2, p1, box_dimension);
	vertices[1][1]=nearest_image(q2, p1, box_dimension);
	vertices[1][2]=nearest_image(r2, p1, box_dimension);
	
	planeoffset[0]=-dot_product(norm1, vertices[0][0]);
	planeoffset[1]=-dot_product(norm2, vertices[1][0]);
	
	for(j=0;j<3;j++) verticesplaneframe[0][j]=dot_product(norm2, vertices[0][j])+planeoffset[1];
	for(j=0;j<3;j++) verticesplaneframe[1][j]=dot_product(norm1, vertices[1][j])+planeoffset[0];
	
	for(i=0;i<2;i++){
		if(verticesplaneframe[i][0]>0){
			if(verticesplaneframe[i][1]>0){
				if(verticesplaneframe[i][2]>0){
					free_matrix(vertices, 2);
					free_matrix(verticesplaneframe, 2);
					free(planeoffset);
					free(singleindex);
					free(onplanecode);
					return 0;					//	all on same side; early reject
				}
				else if(verticesplaneframe[i][2]<0){
					singleindex[i]=2;
				}
				else{
					singleindex[i]=2;
					onplanecode[i]=1;			//	1 denotes single vertex on plane
				}
			}
			else if(verticesplaneframe[i][1]<0){
				if(verticesplaneframe[i][2]>0){
					singleindex[i]=1;
				}
				else if(verticesplaneframe[i][2]<0){
					singleindex[i]=0;
				}
				else{
					singleindex[i]=0;
					onplanecode[i]=5;			//	3+i, where i>=0, denotes that i is on plane and other two are on opposite sides
				}
			}
			else{
				if(verticesplaneframe[i][2]>0){
					singleindex[i]=1;
					onplanecode[i]=1;			//	1 denotes single vertex on plane
				}
				else if(verticesplaneframe[i][2]<0){
					singleindex[i]=0;
					onplanecode[i]=4;			//	3+i, where i>=0, denotes that i is on plane and other two are on opposite sides
				}
				else{
					singleindex[i]=0;
					onplanecode[i]=2;			//	2 denotes that two are on plane
				}
			}
		}
		else if(verticesplaneframe[i][0]<0){
			if(verticesplaneframe[i][1]>0){
				if(verticesplaneframe[i][2]>0){
					singleindex[i]=0;
				}
				else if(verticesplaneframe[i][2]<0){
					singleindex[i]=1;
				}
				else{
					singleindex[i]=0;
					onplanecode[i]=5;			//	3+i, where i>=0, denotes that i is on plane and other two are on opposite sides
				}
			}
			else if(verticesplaneframe[i][1]<0){
				if(verticesplaneframe[i][2]>0){
					singleindex[i]=2;
				}
				else if(verticesplaneframe[i][2]<0){
					free_matrix(vertices, 2);
					free_matrix(verticesplaneframe, 2);
					free(planeoffset);
					free(singleindex);
					free(onplanecode);
					return 0;					//	all on same side; early reject
				}
				else{
					singleindex[i]=2;
					onplanecode[i]=1;			//	1 denotes single vertex on plane
				}
			}
			else{
				if(verticesplaneframe[i][2]>0){
					singleindex[i]=0;
					onplanecode[i]=4;			//	3+i, where i>=0, denotes that i is on plane and other two are on opposite sides
				}
				else if(verticesplaneframe[i][2]<0){
					singleindex[i]=1;
					onplanecode[i]=1;			//	1 denotes single vertex on plane
				}
				else{
					singleindex[i]=0;
					onplanecode[i]=2;			//	2 denotes that two are on plane
				}
			}
		}
		else{
			if(verticesplaneframe[i][1]>0){
				if(verticesplaneframe[i][2]>0){
					singleindex[i]=0;
					onplanecode[i]=1;			//	1 denotes single vertex on plane
				}
				else if(verticesplaneframe[i][2]<0){
					singleindex[i]=1;
					onplanecode[i]=3;			//	3+i, where i>=0, denotes that i is on plane and other two are on opposite sides
				}
				else{
					singleindex[i]=1;
					onplanecode[i]=2;			//	2 denotes that two are on plane
				}
			}
			else if(verticesplaneframe[i][1]<0){
				if(verticesplaneframe[i][2]>0){
					singleindex[i]=1;
					onplanecode[i]=3;			//	3+i, where i>=0, denotes that i is on plane and other two are on opposite sides
				}
				else if(verticesplaneframe[i][2]<0){
					singleindex[i]=0;
					onplanecode[i]=1;			//	1 denotes single vertex on plane
				}
				else{
					singleindex[i]=1;
					onplanecode[i]=2;			//	2 denotes that two are on plane
				}
			}
			else{
				if(verticesplaneframe[i][2]>0){
					singleindex[i]=2;
					onplanecode[i]=2;			//	2 denotes that two are on plane
				}
				else if(verticesplaneframe[i][2]<0){
					singleindex[i]=2;
					onplanecode[i]=2;			//	2 denotes that two are on plane
				}
				else{
					onplanecode[i]=-1;			//	-1 denotes that all are on plane
				}
			}
		}
	}
	
	int result;
	declare_array(int, Nbounds, 2);
	declare_matrix(double, bounds, 2, 2);
	double singleprojection, tmp, verticesplaneframesingle;
	double_triple normcross;
    
	if((onplanecode[0]>=0)&&(onplanecode[1]>=0)){   //  triangles not coplanar
		normcross=cross_product(norm1, norm2);
		for(i=0;i<2;i++){
			if(onplanecode[i]==1){					//	1 on plane
				Nbounds[i]=1;
				bounds[i][0]=dot_product(normcross, vertices[i][singleindex[i]]);
			}
			else if(onplanecode[i]==2){				//	2 on plane
				Nbounds[i]=2;
				bounds[i][0]=dot_product(normcross, vertices[i][mod(singleindex[i]+1, 3)]);
				bounds[i][1]=dot_product(normcross, vertices[i][mod(singleindex[i]+2, 3)]);
			}
			else{
				Nbounds[i]=2;
				singleprojection=dot_product(normcross, vertices[i][singleindex[i]]);
				verticesplaneframesingle=verticesplaneframe[i][singleindex[i]];
                
				bounds[i][0]=singleprojection+(dot_product(normcross, vertices[i][mod(singleindex[i]+1, 3)])-singleprojection)*verticesplaneframesingle/(verticesplaneframesingle-verticesplaneframe[i][mod(singleindex[i]+1, 3)]);
				bounds[i][1]=singleprojection+(dot_product(normcross, vertices[i][mod(singleindex[i]+2, 3)])-singleprojection)*verticesplaneframesingle/(verticesplaneframesingle-verticesplaneframe[i][mod(singleindex[i]+2, 3)]);
			}
		}
		for(i=0;i<2;i++){
			if(Nbounds[i]==2){
				if(bounds[i][0]>bounds[i][1]){		//	reorder to make bounds[i][0]<=bounds[i][1]
					tmp=bounds[i][0];
					bounds[i][0]=bounds[i][1];
					bounds[i][1]=tmp;
				}
			}
		}
		if(Nbounds[0]==1){
			if(Nbounds[1]==1){
                result=0;                           //  only count interior intersections, not pointwise
			}
			else{
				if(bounds[0][0]<=bounds[1][0]) result=0;
				else if(bounds[0][0]>=bounds[1][1]) result=0;
				else result=1;
			}
		}
		else if(Nbounds[1]==1){
			if(bounds[1][0]<=bounds[0][0]) result=0;
			else if(bounds[1][0]>=bounds[0][1]) result=0;
			else result=1;
		}
		else{										//	general case: two bounds each
			if(bounds[0][0]<=bounds[1][0]){
				if(bounds[0][1]<=bounds[1][0]) result=0;
				else result=1;
			}
			else{
				if(bounds[0][0]>=bounds[1][1]) result=0;
				else result=1;
			}
		}
		free(Nbounds);
		free_matrix(bounds, 2);
		free_matrix(vertices, 2);
		free_matrix(verticesplaneframe, 2);
		free(planeoffset);
		free(singleindex);
		free(onplanecode);
		return result;
	}
	
	//	planar calculation
	
    for(i=0;i<2;i++){
        for(j=0;j<3;j++) printf("\t\t%f, %f, %f\n", vertices[i][j].x, vertices[i][j].y, vertices[i][j].z);
    }
    
	double_triple x, y;			//	arbitrary vectors perpendicular to norm1
	if(norm1.x==1){
		x.x=0;
		x.y=1;
		x.z=0;
	}
	else{
		x.x=1;
		x.y=0;
		x.z=0;
	}
	y=cross_product(norm1, x);
	y=normed(y);
	x=cross_product(y, norm1);
	x=normed(x);
	declare_matrix_nozero(double_double, projectedvertices, 2, 3);
	for(i=0;i<2;i++){
		for(j=0;j<3;j++){
			projectedvertices[i][j].x=dot_product(x, vertices[i][j]);
			projectedvertices[i][j].y=dot_product(y, vertices[i][j]);
		}
	}
    for(i=0;i<2;i++){
        for(j=0;j<3;j++) printf("\t\t%f, %f\n", projectedvertices[i][j].x, projectedvertices[i][j].y);
    }
	for(i=0;i<3;i++){
		for(j=0;j<3;j++){
			if(checkedgecrossinplane(projectedvertices[0][i], projectedvertices[0][mod(i+1, 3)], projectedvertices[1][j], projectedvertices[1][mod(j+1, 3)])==1){
				free_matrix(projectedvertices, 2);
				free(Nbounds);	
				free_matrix(bounds, 2);
				free_matrix(vertices, 2);
				free_matrix(verticesplaneframe, 2);
				free(planeoffset);
				free(singleindex);
				free(onplanecode);
				return 1;
			}
		}
	}
	
	//	check a point from each triangle: see if it is contained within other triangle
    
    if(checkpointintriangle(projectedvertices[0][0], projectedvertices[1])==1){
        free_matrix(projectedvertices, 2);
        free(Nbounds);
        free_matrix(bounds, 2);
        free_matrix(vertices, 2);
        free_matrix(verticesplaneframe, 2);
        free(planeoffset);
        free(singleindex);
        free(onplanecode);
        return 1;
    }
    if(checkpointintriangle(projectedvertices[1][0], projectedvertices[0])==1){
        free_matrix(projectedvertices, 2);
        free(Nbounds);
        free_matrix(bounds, 2);
        free_matrix(vertices, 2);
        free_matrix(verticesplaneframe, 2);
        free(planeoffset);
        free(singleindex);
        free(onplanecode);
        return 1;
    }
	free_matrix(projectedvertices, 2);
	free(Nbounds);	
	free_matrix(bounds, 2);
	free_matrix(vertices, 2);
	free_matrix(verticesplaneframe, 2);
	free(planeoffset);
	free(singleindex);
	free(onplanecode);
	return 0;
}

int checkpointintriangle(double_double point, double_double *vertices){
    double mag;
    double_double parallel, perp, thirdlab, third, pointframe;
    parallel.x=vertices[1].x-vertices[0].x;
    parallel.y=vertices[1].y-vertices[0].y;
    mag=sqrt(parallel.x*parallel.x+parallel.y*parallel.y);
    parallel.x/=mag;
    parallel.y/=mag;
    perp.x=parallel.y;
    perp.y=-parallel.x;
    thirdlab.x=vertices[2].x-vertices[0].x;
    thirdlab.y=vertices[2].y-vertices[0].y;
    third.x=parallel.x*thirdlab.x+parallel.y*thirdlab.y;
    third.y=perp.x*thirdlab.x+perp.y*thirdlab.y;
    point.x=point.x-vertices[0].x;
    point.y=point.y-vertices[0].y;
    pointframe.x=parallel.x*point.x+parallel.y*point.y;
    pointframe.y=perp.x*point.x+perp.y*point.y;
    if(third.y<0){
        third.y*=-1;
        pointframe.y*=-1;       //  mirror reflect
    }
    if(pointframe.y<=0) return 0;        //  outside of 0-1 edge
    if(pointframe.x<=third.x+pointframe.y*third.x/third.y) return 0;     //  outside of 0-2 edge
    if(pointframe.x>=mag+pointframe.y*(third.x-mag)/third.y) return 0;   //  outside of 1-2 edge
    return 1;
}

int checkpointintrianglelist(double_double point, double_double v0, double_double v1, double_double v2){
    double mag;
    double_double parallel, perp, thirdlab, third, pointframe;
    parallel.x=v1.x-v0.x;
    parallel.y=v1.y-v0.y;
    mag=sqrt(parallel.x*parallel.x+parallel.y*parallel.y);
    parallel.x/=mag;
    parallel.y/=mag;
    perp.x=parallel.y;
    perp.y=-parallel.x;
    thirdlab.x=v2.x-v0.x;
    thirdlab.y=v2.y-v0.y;
    third.x=parallel.x*thirdlab.x+parallel.y*thirdlab.y;
    third.y=perp.x*thirdlab.x+perp.y*thirdlab.y;
    point.x=point.x-v0.x;
    point.y=point.y-v0.y;
    pointframe.x=parallel.x*point.x+parallel.y*point.y;
    pointframe.y=perp.x*point.x+perp.y*point.y;
    if(third.y<0){
        third.y*=-1;
        pointframe.y*=-1;       //  mirror reflect
    }
    if(pointframe.y<=0) return 0;        //  outside of 0-1 edge
    if(pointframe.x<=third.x+pointframe.y*third.x/third.y) return 0;     //  outside of 0-2 edge
    if(pointframe.x>=mag+pointframe.y*(third.x-mag)/third.y) return 0;   //  outside of 1-2 edge
    return 1;
}

int checkedgecrossinplane(double_double a0, double_double a1, double_double b0, double_double b1){
	double amag, intercept;
	double_double parallela, perpa, b0framea, b1framea;
	parallela.x=a1.x-a0.x;
	parallela.y=a1.y-a0.y;
	amag=sqrt(parallela.x*parallela.x+parallela.y*parallela.y);
	parallela.x/=amag;
	parallela.y/=amag;
	perpa.x=parallela.y;
	perpa.y=-parallela.x;
	
	//	put b0, b1 in frame with a0 at origin
	
	b0.x-=a0.x;
	b0.y-=a0.y;
	b1.x-=a0.x;
	b1.y-=a0.y;
	
	//	rotate to frame with a1 in x direction
	
	b0framea.x=b0.x*parallela.x+b0.y*parallela.y;
	b0framea.y=b0.x*perpa.x+b0.y*perpa.y;
	b1framea.x=b1.x*parallela.x+b1.y*parallela.y;
	b1framea.y=b1.x*perpa.x+b1.y*perpa.y;
	
	if(b0framea.y>=0){
		if(b1framea.y>=0) return 0;			//	both on same side
	}
	if(b0framea.y<=0){
		if(b1framea.y<=0) return 0;			//	both on same side
	}
	intercept=b0framea.x-b0framea.y/(b1framea.y-b0framea.y)*(b1framea.x-b0framea.x);
    if((intercept>0)&&(intercept<amag)) return 1;
    else return 0;
}

void free_mesh_maps(int_double meshsize, int **mapmeshtotriangles, int **mapmeshtobonds, int **maptrianglestomesh, int **mapbondstomesh, int **mapmeshtoneighbors, int **maptrianglestotriangles, int **mapmeshtoouterbonds, int **mapmeshtooutertriangles, int **mapbondstotriangles){
    int meshnumber=2*meshsize.x*meshsize.y, i;
	free_matrix(mapmeshtotriangles, meshnumber);
	free_matrix(mapmeshtobonds, meshnumber);
	free_matrix(maptrianglestomesh, 2*meshnumber);
	free_matrix(mapbondstomesh, 3*meshnumber);
	free_matrix(mapmeshtoneighbors, meshnumber);
	free_matrix(maptrianglestotriangles, 2*meshnumber);
	free_matrix(mapmeshtoouterbonds, meshnumber);
	free_matrix(mapmeshtooutertriangles, meshnumber);
	free_matrix(mapbondstotriangles, 3*meshnumber);
}

void initialize_mesh_maps(int_double meshsize, int ***pmapmeshtotriangles, int ***pmapmeshtobonds, int ***pmaptrianglestomesh, int ***pmapbondstomesh, int ***pmapmeshtoneighbors, int ***pmaptrianglestotriangles, int ***pmapmeshtoouterbonds, int ***pmapmeshtooutertriangles, int ***pmapbondstotriangles){
	int meshnumber=2*meshsize.x*meshsize.y, interface, row, column, indexshift, meshpoint, i, j;
	allocate_matrix(int, (*pmapmeshtotriangles), meshnumber, 6);
	allocate_matrix(int, (*pmapmeshtobonds), meshnumber, 6);
	allocate_matrix(int, (*pmaptrianglestomesh), 2*meshnumber, 3);
	allocate_matrix(int, (*pmapbondstomesh), 3*meshnumber, 2);
	allocate_matrix(int, (*pmapmeshtoneighbors), meshnumber, 6);
	allocate_matrix(int, (*pmaptrianglestotriangles), 2*meshnumber, 12);
	allocate_matrix(int, (*pmapmeshtoouterbonds), meshnumber, 6);
	allocate_matrix(int, (*pmapmeshtooutertriangles), meshnumber, 6);
	allocate_matrix(int, (*pmapbondstotriangles), 3*meshnumber, 2);
	for(i=0;i<2*meshnumber;i++){
		for(j=0;j<3;j++){
			(*pmaptrianglestomesh)[i][j]=-1;
		}
		for(j=0;j<12;j++){
			(*pmaptrianglestotriangles)[i][j]=-1;
		}
	}
	for(i=0;i<3*meshnumber;i++){
		for(j=0;j<2;j++){
			(*pmapbondstomesh)[i][j]=-1;
		}
    }
    
	//	in right-handed order
	
	for(interface=0;interface<2;interface++){
		indexshift=interface*meshsize.x*meshsize.y;
		for(row=0;row<meshsize.y;row++){
			for(column=0;column<meshsize.x;column++){
				meshpoint=indexshift+row*meshsize.x+column;
                
				(*pmapmeshtotriangles)[meshpoint][0]=2*meshpoint;
				(*pmaptrianglestomesh)[2*meshpoint][2]=meshpoint;
				(*pmapmeshtobonds)[meshpoint][0]=3*meshpoint;
				(*pmapbondstomesh)[3*meshpoint][0]=meshpoint;
				
				(*pmapmeshtotriangles)[meshpoint][1]=2*meshpoint+1;
				(*pmaptrianglestomesh)[2*meshpoint+1][0]=meshpoint;
				(*pmapmeshtobonds)[meshpoint][1]=3*meshpoint+1;
				(*pmapbondstomesh)[3*meshpoint+1][0]=meshpoint;
				
				if(column==meshsize.x-1) (*pmapmeshtoneighbors)[meshpoint][0]=indexshift+row*meshsize.x;
				else (*pmapmeshtoneighbors)[meshpoint][0]=meshpoint+1;
                
				if(row%2==0) (*pmapmeshtoneighbors)[meshpoint][1]=meshpoint+meshsize.x;
				else{
					if(column==meshsize.x-1){
						if(row==meshsize.y-1) (*pmapmeshtoneighbors)[meshpoint][1]=indexshift;
						else (*pmapmeshtoneighbors)[meshpoint][1]=indexshift+(row+1)*meshsize.x;
					}
					else{
						if(row==meshsize.y-1) (*pmapmeshtoneighbors)[meshpoint][1]=indexshift+column+1;
						else (*pmapmeshtoneighbors)[meshpoint][1]=meshpoint+meshsize.x+1;
					}
				}
				
				if(row%2==0){
					if(column==0) (*pmapmeshtoneighbors)[meshpoint][2]=indexshift+(row+2)*meshsize.x-1;
					else (*pmapmeshtoneighbors)[meshpoint][2]=meshpoint+meshsize.x-1;
				}
				else{
					if(row==meshsize.y-1) (*pmapmeshtoneighbors)[meshpoint][2]=indexshift+column;
					else (*pmapmeshtoneighbors)[meshpoint][2]=meshpoint+meshsize.x;
				}
				
				if(column==0) (*pmapmeshtoneighbors)[meshpoint][3]=meshpoint+meshsize.x-1;
				else (*pmapmeshtoneighbors)[meshpoint][3]=meshpoint-1;
				(*pmapmeshtotriangles)[meshpoint][2]=2*(*pmapmeshtoneighbors)[meshpoint][3];
				(*pmaptrianglestomesh)[2*(*pmapmeshtoneighbors)[meshpoint][3]][0]=meshpoint;
				(*pmapmeshtobonds)[meshpoint][2]=3*meshpoint+2;
				(*pmapbondstomesh)[3*meshpoint+2][0]=meshpoint;
				
				if(row%2==0){
					if(column==0){
						if(row==0) (*pmapmeshtoneighbors)[meshpoint][4]=indexshift+meshsize.x*meshsize.y-1;
						else (*pmapmeshtoneighbors)[meshpoint][4]=meshpoint-1;
					}
					else{
						if(row==0) (*pmapmeshtoneighbors)[meshpoint][4]=indexshift+meshsize.x*(meshsize.y-1)+column-1;
						else (*pmapmeshtoneighbors)[meshpoint][4]=meshpoint-meshsize.x-1;
					}
				}
				else (*pmapmeshtoneighbors)[meshpoint][4]=meshpoint-meshsize.x;
				(*pmapmeshtotriangles)[meshpoint][3]=2*(*pmapmeshtoneighbors)[meshpoint][4]+1;
				(*pmaptrianglestomesh)[2*(*pmapmeshtoneighbors)[meshpoint][4]+1][1]=meshpoint;
                
				(*pmapmeshtobonds)[meshpoint][3]=3*(*pmapmeshtoneighbors)[meshpoint][3];
				(*pmapbondstomesh)[3*(*pmapmeshtoneighbors)[meshpoint][3]][1]=meshpoint;
                
				(*pmapmeshtotriangles)[meshpoint][4]=2*(*pmapmeshtoneighbors)[meshpoint][4];
				(*pmaptrianglestomesh)[2*(*pmapmeshtoneighbors)[meshpoint][4]][1]=meshpoint;
                
                (*pmapmeshtobonds)[meshpoint][4]=3*(*pmapmeshtoneighbors)[meshpoint][4]+1;
				(*pmapbondstomesh)[3*(*pmapmeshtoneighbors)[meshpoint][4]+1][1]=meshpoint;
                
				if(row%2==0){
					if(row==0) (*pmapmeshtoneighbors)[meshpoint][5]=indexshift+meshsize.x*(meshsize.y-1)+column;
					else (*pmapmeshtoneighbors)[meshpoint][5]=meshpoint-meshsize.x;
				}
				else{
					if(column==meshsize.x-1) (*pmapmeshtoneighbors)[meshpoint][5]=indexshift+(row-1)*meshsize.x;
					else (*pmapmeshtoneighbors)[meshpoint][5]=meshpoint-meshsize.x+1;
				}
				(*pmapmeshtotriangles)[meshpoint][5]=2*(*pmapmeshtoneighbors)[meshpoint][5]+1;
				(*pmaptrianglestomesh)[2*(*pmapmeshtoneighbors)[meshpoint][5]+1][2]=meshpoint;
                
                (*pmapmeshtobonds)[meshpoint][5]=3*(*pmapmeshtoneighbors)[meshpoint][5]+2;
				(*pmapbondstomesh)[3*(*pmapmeshtoneighbors)[meshpoint][5]+2][1]=meshpoint;
                
            }
		}
	}
    for(i=0;i<meshnumber;i++){
        for(j=0;j<6;j++){
            (*pmapmeshtoouterbonds)[i][j]=(*pmapmeshtobonds)[(*pmapmeshtoneighbors)[i][j]][mod(j+2, 6)];
            (*pmapmeshtooutertriangles)[i][j]=(*pmapmeshtotriangles)[(*pmapmeshtoneighbors)[i][j]][mod(j+1, 6)];
        }
    }   
    int bond;
    for(i=0;i<meshnumber;i++){
        for(j=0;j<3;j++){
            bond=(*pmapmeshtobonds)[i][j];
            (*pmapbondstotriangles)[bond][0]=(*pmapmeshtotriangles)[i][mod(j-1, 6)];
            (*pmapbondstotriangles)[bond][1]=(*pmapmeshtotriangles)[i][j];
        }
    }
    
	int pointneighbor, lastindex, k;
	for(i=0;i<meshnumber;i++){
		for(j=0;j<6;j++){
			if((*pmapmeshtoneighbors)[i][j]==-1) my_exit("neighbors bad");
            if((*pmapmeshtoneighbors)[i][j]==i){
                printf("mapmeshtoneighbors[%i][%i]=%i!\n", i, j, i);
                exit(1);
            }
            for(k=0;k<j;k++) if((*pmapmeshtoneighbors)[i][j]==(*pmapmeshtoneighbors)[i][k]) my_exit("repeat neighbors");
		}
	}
	for(i=0;i<2*meshnumber;i++){
		for(j=0;j<3;j++){
			if((*pmaptrianglestomesh)[i][j]==-1) my_exit("map bad");
		}
	}
	for(i=0;i<3*meshnumber;i++){
		for(j=0;j<2;j++){
			if((*pmapbondstomesh)[i][j]==-1) my_exit("map bad");
		}
	}
	for(i=0;i<2*meshnumber;i++){
		for(j=0;j<3;j++){
			pointneighbor=(*pmaptrianglestomesh)[i][j];
			lastindex=0;
			while((*pmapmeshtotriangles)[pointneighbor][lastindex]!=i){
				lastindex++;
			}			
			if(j==0){
				for(k=0;k<5;k++){
					(*pmaptrianglestotriangles)[i][k]=(*pmapmeshtotriangles)[pointneighbor][mod(lastindex+k+1, 6)];
				}
			}
			else if(j==1){
				for(k=0;k<4;k++){
					(*pmaptrianglestotriangles)[i][5+k]=(*pmapmeshtotriangles)[pointneighbor][mod(lastindex+k+2, 6)];
				}
			}
			else{
				for(k=0;k<3;k++){
					(*pmaptrianglestotriangles)[i][9+k]=(*pmapmeshtotriangles)[pointneighbor][mod(lastindex+k+2, 6)];
				}
			}
		}
	}
	for(i=0;i<2*meshnumber;i++){
		for(j=0;j<12;j++){
			if((*pmaptrianglestotriangles)[i][j]==-1) my_exit("triangles to triangles bad");
            if((*pmaptrianglestotriangles)[i][j]==i) my_exit("triangles mapping to self");
            for(k=0;k<j;k++) if((*pmaptrianglestotriangles)[i][j]==(*pmaptrianglestotriangles)[i][k]) my_exit("repeat triangles");
		}
	}
}

int checkoverlap(double_triple *meshpositions, double_triple newpos, int meshnumber, int_double meshsize, int mover, int *neighbors, double_triple box_dimension, linkedlistfull meshmeshlink, triangledata *meshtriangledata, triangledata *newtriangledata, int **mapmeshtotriangles, int **maptrianglestotriangles, int **maptrianglestomesh){
	int i, j, k, foundneighbors, leftneighbor, candidate, foundtriangles, centraltriangle, neighbor, meshpoint;
	double_triple rightpos, leftpos;
	declare_array(int, myneighbors, meshmeshlink.maxneighbors);
	declare_array(int, mytriangles, 2*meshmeshlink.maxneighbors);		//	typically 2x more triangles than neighbors
	declare_array(int, found, meshnumber);
	declare_array(int, trianglefound, 2*meshnumber);
	for(i=0;i<6;i++){
		foundneighbors=foundtriangles=0;
		if(i==5) leftneighbor=0;
        else leftneighbor=i+1;
		rightpos=meshpositions[neighbors[i]];
		leftpos=meshpositions[neighbors[leftneighbor]];
        
		//	find close mesh points that aren't part of moving triangle
		
        found[mover]=1;
		found[neighbors[i]]=1;
		found[neighbors[leftneighbor]]=1;
		for(j=0;j<meshmeshlink.core.number_neighbors[mover];j++){
			candidate=meshmeshlink.core.neighbor[mover][j];
			if(found[candidate]==0){
				myneighbors[foundneighbors]=candidate;
				found[candidate]=1;
				foundneighbors++;
			}
		}
		for(j=0;j<meshmeshlink.core.number_neighbors[neighbors[i]];j++){
			candidate=meshmeshlink.core.neighbor[neighbors[i]][j];
			if(found[candidate]==0){
				myneighbors[foundneighbors]=candidate;
				found[candidate]=1;
				foundneighbors++;
			}
		}
		for(j=0;j<meshmeshlink.core.number_neighbors[neighbors[leftneighbor]];j++){
			candidate=meshmeshlink.core.neighbor[neighbors[leftneighbor]][j];
			if(found[candidate]==0){
				myneighbors[foundneighbors]=candidate;
				found[candidate]=1;
				foundneighbors++;
			}
		}
		
		//	find triangles that border these mesh points and do not border central triangle
        
        centraltriangle=mapmeshtotriangles[mover][i];
        for(j=0;j<12;j++) trianglefound[maptrianglestotriangles[centraltriangle][j]]=1;    //  don't calculate overlap with triangles that touch central triangle
		foundtriangles=0;
		for(j=0;j<foundneighbors;j++){
            neighbor=myneighbors[j];
            for(k=0;k<6;k++){
                candidate=mapmeshtotriangles[neighbor][k];
                if(trianglefound[candidate]==0){
                    if(checktriangleoverlap(newpos, rightpos, leftpos, meshpositions[maptrianglestomesh[candidate][0]], meshpositions[maptrianglestomesh[candidate][1]], meshpositions[maptrianglestomesh[candidate][2]], newtriangledata[i].normal, meshtriangledata[candidate].normal, box_dimension)==1){
                        free(myneighbors);
                        free(mytriangles);
                        free(found);
                        free(trianglefound);
                        return 1;
                    }
                    trianglefound[candidate]=1;
                    mytriangles[foundtriangles]=candidate;
                    foundtriangles++;
                }
			}
		}
        for(j=0;j<foundneighbors;j++) found[myneighbors[j]]=0;				//	clearing
        found[mover]=0;
		found[neighbors[i]]=0;
		found[neighbors[leftneighbor]]=0;
        for(j=0;j<foundtriangles;j++) trianglefound[mytriangles[j]]=0;		//	clearing
        for(j=0;j<6;j++) trianglefound[maptrianglestotriangles[mapmeshtotriangles[mover][i]][j]]=0;
	}
	free(myneighbors);
	free(mytriangles);
	free(found);
	free(trianglefound);
	return 0;
}

int mc_mesh(double_triple *meshpositions, solvation_parameters *psolvp, double max_translate, double_triple *pbox_dimension, double temperature, coord *coordarray, linkedlistasymmetric *pmeshlink, int Nnodes, nonbonded_params my_nonbonded_params, linkedlistfull *pphenlink, int *pmeshflag, triangledata *meshtriangledata, bonddata *meshbonddata, vertexdata *meshdata, linkedlistfull *pmeshmeshlink, int **mapmeshtotriangles, int **mapmeshtobonds, int **mapmeshtoneighbors, int **maptrianglestotriangles, int **maptrianglestomesh, int **mapmeshtoouterbonds, int **mapmeshtooutertriangles){
    int mover, leftneighbor, i, j, k, interfaceindex, mynodetype, siteneighbor, flag=0, myphenphenneighbor, outerbond;
    if((*psolvp).meshmovecode==0) mover=rand_int((*psolvp).meshnumber);
    else if(((*psolvp).meshmovecode==1)||((*psolvp).meshmovecode==2)) mover=(*psolvp).meshnumber/2+rand_int((*psolvp).meshnumber/2);
    else my_exit("code not written for this meshmovecode!");
    if(mover<(*psolvp).meshnumber/2) interfaceindex=0;
    else interfaceindex=1;
	double_triple shift;
    if((*psolvp).meshmovecode==2){
        shift.x=shift.y=0;
        shift.z=(rand_double*2-1)*0.6*max_translate;
    }
    else shift=scalar_multiply_double_triple(rand_unit_cube(), (rand_int(2)*2-1)*0.6*max_translate);
    double_triple oldpos=meshpositions[mover];
    double_triple newpos=add_double_triple(oldpos, shift);
    onebodyparam my_onebody_param;
    
	double energy, acceptance_prob, totalarea=0, oldarea=0, volumechange=0, newheight, curvatureterm, curvaturesquaredterm;
    
    //  calculate area change, volume change, uniformbond change
    
    double_triple bond;
	declare_array_nozero(triangledata, newtriangledata, 6);
	declare_array_nozero(bonddata, newbonddata, 6);
    declare_array(double, new_outerbonds_meancurvature, 6);
    vertexdata newcentralvertexdata;
    declare_array_nozero(vertexdata, newoutervertexdata, 6);
    for(i=0;i<6;i++){
        if(i==5) leftneighbor=0;
        else leftneighbor=i+1;
		newtriangledata[i]=calc_triangledata(newpos, meshpositions[mapmeshtoneighbors[mover][i]], meshpositions[mapmeshtoneighbors[mover][leftneighbor]], (*pbox_dimension));
		totalarea+=newtriangledata[i].area;
		volumechange+=newtriangledata[i].columnvolume;
		oldarea+=meshtriangledata[mapmeshtotriangles[mover][i]].area;
		volumechange-=meshtriangledata[mapmeshtotriangles[mover][i]].columnvolume;
		
        bond=subtract_double_triple(meshpositions[mapmeshtoneighbors[mover][i]], newpos);
        recenter_double_triple(&bond, (*pbox_dimension));
        newbonddata[i].length=norm(bond);
	}
    for(i=0;i<6;i++){
        newbonddata[i].meancurvature=newbonddata[i].length*mysafeacos(dot_product(newtriangledata[i].normal, newtriangledata[mod(i-1, 6)].normal));
        outerbond=mapmeshtoouterbonds[mover][i];
        new_outerbonds_meancurvature[i]=meshbonddata[outerbond].length*mysafeacos(dot_product(newtriangledata[i].normal, meshtriangledata[mapmeshtooutertriangles[mover][i]].normal));
    }
    newcentralvertexdata.area=0;
    newcentralvertexdata.meancurvature=0;
    for(i=0;i<6;i++){
        newcentralvertexdata.area+=newtriangledata[i].area;
        newcentralvertexdata.meancurvature+=newbonddata[i].meancurvature;
    }
    newcentralvertexdata.meancurvaturesquared=pow(newcentralvertexdata.meancurvature, 2)/newcentralvertexdata.area;
    for(i=0;i<6;i++){
        newoutervertexdata[i]=meshdata[mapmeshtoneighbors[mover][i]];
        newoutervertexdata[i].area+=(newtriangledata[i].area+newtriangledata[mod(i-1, 6)].area-meshtriangledata[mapmeshtotriangles[mover][i]].area-meshtriangledata[mapmeshtotriangles[mover][mod(i-1, 6)]].area);
        newoutervertexdata[i].meancurvature+=(newbonddata[i].meancurvature-meshbonddata[mapmeshtobonds[mover][i]].meancurvature);
        newoutervertexdata[i].meancurvature+=(new_outerbonds_meancurvature[i]-meshbonddata[mapmeshtoouterbonds[mover][i]].meancurvature);
        newoutervertexdata[i].meancurvature+=(new_outerbonds_meancurvature[mod(i-1, 6)]-meshbonddata[mapmeshtoouterbonds[mover][mod(i-1, 6)]].meancurvature);
        newoutervertexdata[i].meancurvaturesquared=pow(newoutervertexdata[i].meancurvature, 2)/newoutervertexdata[i].area;
    }
    
    if(interfaceindex==0) volumechange*=(-1);           //  column volume measured in +z direction, so need to switch signs for bottom-facing interface
    energy=0;
	energy+=(totalarea-oldarea)*(*psolvp).surfacetension;
	energy+=((*psolvp).kvolumerestraint*(pow((*psolvp).totalvolume+volumechange-(*psolvp).targetvolume, 2)-pow((*psolvp).totalvolume-(*psolvp).targetvolume, 2)));
    
    //  Curvature terms
    
    curvatureterm=newcentralvertexdata.meancurvature-meshdata[mover].meancurvature;
    curvaturesquaredterm=newcentralvertexdata.meancurvaturesquared-meshdata[mover].meancurvaturesquared;
    for(i=0;i<6;i++){
        curvatureterm+=(new_outerbonds_meancurvature[i]-meshbonddata[mapmeshtoouterbonds[mover][i]].meancurvature);
        curvaturesquaredterm+=(newoutervertexdata[i].meancurvaturesquared-meshdata[mapmeshtoneighbors[mover][i]].meancurvaturesquared);
    }
    energy+=(*psolvp).surfacetension*(-0.5*(*psolvp).Tolmanlength*curvatureterm+(0.1875*pow((*psolvp).meancurvaturelength, 2)*curvaturesquaredterm));
    
    //  calculate uniform bond energy
    
    double totallength, oldtotallength, uniformbondenergy=0;
    declare_array(double, newtotallength, 3);
    
    for(i=0;i<3;i++){
        totallength=oldtotallength=0;
        for(j=0;j<2;j++){
            totallength+=newbonddata[i+3*j].length;
            oldtotallength+=meshbonddata[mapmeshtobonds[mover][i+3*j]].length;
        }
        newtotallength[i]=(*psolvp).totallength[interfaceindex*3+i]-oldtotallength+totallength;
        for(j=0;j<2;j++){
            uniformbondenergy+=((pow((*psolvp).meshnumber/2*newbonddata[i+3*j].length/newtotallength[i]-1, 2)-pow((*psolvp).meshnumber/2*meshbonddata[mapmeshtobonds[mover][i+3*j]].length/(*psolvp).totallength[interfaceindex*3+i]-1, 2)));
        }
    }
    
    energy+=(*psolvp).kuniformbonds*uniformbondenergy;
    
    int Nsiteneighborsthatchangeheight=0, Nphenylsiteneighborsthatchangeheight=0;
    int *siteneighborsthatchangeheight;
    int *phenylsiteneighborsthatchangeheight;
    double *newheights;
    int *phenylheightchanged;
    double *newphenylheight;
    double **newphenphenenergies;
    if(Nnodes>0){
        
        //	identify site neighbors of mover and neighbors of mover
        
        declare_array(int, siteneighbors, 3*(*pmeshlink).maxreverseneighbors);		//	3x, not 7x bigger
        int Nsiteneighbors=0;
        for(i=0;i<(*pmeshlink).number_neighbors_reverse[mover];i++){
            siteneighbors[Nsiteneighbors]=(*pmeshlink).neighbor_reverse[mover][i];
            Nsiteneighbors++;
        }
        for(i=0;i<6;i++){
            for(j=0;j<(*pmeshlink).number_neighbors_reverse[mapmeshtoneighbors[mover][i]];j++){
                siteneighbor=(*pmeshlink).neighbor_reverse[mapmeshtoneighbors[mover][i]][j];
                for(k=0;k<Nsiteneighbors;k++){
                    if(siteneighbors[k]==siteneighbor) break;
                }
                if(k==Nsiteneighbors){		//	siteneighbor not already in siteneighbors list
                    siteneighbors[Nsiteneighbors]=siteneighbor;
                    Nsiteneighbors++;
                    if(Nsiteneighbors==3*(*pmeshlink).maxreverseneighbors) my_exit("too many site neighbors in mc_mesh");
                }
            }
        }
        
        //	identify site neighbors that change heights; add change in onebody energy
        
        allocate_array(int, siteneighborsthatchangeheight, Nsiteneighbors);
        allocate_array(int, phenylsiteneighborsthatchangeheight, Nsiteneighbors);
        allocate_array(double, newheights, Nsiteneighbors);
        allocate_array(int, phenylheightchanged, Nnodes);
        allocate_array(double, newphenylheight, Nnodes);
        for(i=0;i<Nsiteneighbors;i++){
            siteneighbor=siteneighbors[i];
            
            //	check that stored height is correct
            
            meshpositions[mover]=newpos;
			if((*psolvp).interface==1) newheight=height_relative_to_interface(siteneighbor, coordarray[siteneighbor], (*psolvp), (*pmeshlink), (*pbox_dimension), meshpositions, mapmeshtoneighbors);
            meshpositions[mover]=oldpos;
            if(fabs(newheight-coordarray[siteneighbor].height)>0.0000000001){
                siteneighborsthatchangeheight[Nsiteneighborsthatchangeheight]=siteneighbor;
                newheights[Nsiteneighborsthatchangeheight]=newheight;
                Nsiteneighborsthatchangeheight++;
                
                //	expanding out chooseonebody
                
                mynodetype=coordarray[siteneighbor].nodetype;        //  expanding out chooseonebodyose                if(mynodetype==0) my_onebody_param=my_nonbonded_params.one0;
                if(mynodetype==0) my_onebody_param=my_nonbonded_params.one0;
                else if(mynodetype==1){
                    my_onebody_param=my_nonbonded_params.one1;
                    phenylsiteneighborsthatchangeheight[Nphenylsiteneighborsthatchangeheight]=siteneighbor;
                    phenylheightchanged[siteneighbor]=1;
                    newphenylheight[siteneighbor]=newheight;
                    Nphenylsiteneighborsthatchangeheight++;
                }
                else if(mynodetype==2) my_onebody_param=my_nonbonded_params.one2;
                else if(mynodetype==3) my_onebody_param=my_nonbonded_params.one3;
                else{
                    printf("don't have onebody param for type %i!", mynodetype);
                    exit(1);
                }
                
                energy+=calc_onebody_difference_heights(coordarray[siteneighbor].height, newheight, my_onebody_param, (*psolvp));
            }
        }
        free(siteneighbors);
        
        //  calculate and store two-body (phen-phen) interactions
        
        if(Nphenylsiteneighborsthatchangeheight>0){
            allocate_matrix(double, newphenphenenergies, Nphenylsiteneighborsthatchangeheight, (*pphenlink).maxneighbors);
            for(i=0;i<Nphenylsiteneighborsthatchangeheight;i++){
                siteneighbor=phenylsiteneighborsthatchangeheight[i];
                for(j=0;j<(*pphenlink).core.number_neighbors[siteneighbor];j++){
                    myphenphenneighbor=(*pphenlink).core.neighbor[siteneighbor][j];
                    newphenphenenergies[i][j]=phenylphenyl_energy_heightchanged(coordarray[siteneighbor], coordarray[myphenphenneighbor], my_nonbonded_params.p11vac, (*pbox_dimension), my_nonbonded_params.phenylcutoff2, my_nonbonded_params.solvationparams, phenylheightchanged[siteneighbor], newphenylheight[siteneighbor], phenylheightchanged[myphenphenneighbor], newphenylheight[myphenphenneighbor]);
                    energy+=newphenphenenergies[i][j];
                    energy-=(*pphenlink).core.pair_energy[siteneighbor][j];
                }
            }
        }
        free(phenylheightchanged);
        free(newphenylheight);
    }
	
	acceptance_prob=exp(-energy/temperature);
	*pmeshflag=0;
	if((acceptance_prob>1.0)||(rand_double<acceptance_prob)){
        
        //  now check for (rare) mesh overlap
        
        acceptance_prob=1;
        if(checkoverlap(meshpositions, newpos, (*psolvp).meshnumber, (*psolvp).meshsize, mover, mapmeshtoneighbors[mover], (*pbox_dimension), *pmeshmeshlink, meshtriangledata, newtriangledata, mapmeshtotriangles, maptrianglestotriangles, maptrianglestomesh)==1){
            acceptance_prob=0;
        }
        
        if(acceptance_prob>0){
            
            //	update mover, longest mesh bond
            
            for(i=0;i<6;i++) if(newbonddata[i].length>(*psolvp).longestmeshbond) (*psolvp).longestmeshbond=newbonddata[i].length;
            if((*psolvp).longestmeshbond>(*psolvp).maxmeshbond) *pmeshflag=1;
            flag=1;
            fmod_double_triple(&newpos, (*pbox_dimension));
            meshpositions[mover]=newpos;
            
            //  update triangle data
            
            for(i=0;i<6;i++){
				meshtriangledata[mapmeshtotriangles[mover][i]]=newtriangledata[i];
				meshbonddata[mapmeshtobonds[mover][i]]=newbonddata[i];
                meshbonddata[mapmeshtoouterbonds[mover][i]].meancurvature=new_outerbonds_meancurvature[i];
			}
            meshdata[mover]=newcentralvertexdata;
            for(i=0;i<6;i++){
                meshdata[mapmeshtoneighbors[mover][i]]=newoutervertexdata[i];
            }
            
            //	update heights and phenyl pair energies
            
            if(Nsiteneighborsthatchangeheight>0){
                for(i=0;i<Nsiteneighborsthatchangeheight;i++){
                    siteneighbor=siteneighborsthatchangeheight[i];
                    coordarray[siteneighbor].height=newheights[i];
                }
                if(Nphenylsiteneighborsthatchangeheight>0){
                    for(i=0;i<Nphenylsiteneighborsthatchangeheight;i++){
                        siteneighbor=phenylsiteneighborsthatchangeheight[i];
                        for(j=0;j<(*pphenlink).core.number_neighbors[siteneighbor];j++){
                            (*pphenlink).core.pair_energy[siteneighbor][j]=newphenphenenergies[i][j];
                            myphenphenneighbor=(*pphenlink).core.neighbor[siteneighbor][j];
                            k=0;
                            while(((*pphenlink).core.neighbor[myphenphenneighbor][k]!=siteneighbor)&&(k!=(*pphenlink).maxneighbors)) k++;
                            if(k==(*pphenlink).maxneighbors) my_exit("reverse phen neighbor not found in mc_mesh!");
                            (*pphenlink).core.pair_energy[myphenphenneighbor][k]=newphenphenenergies[i][j];
                        }
                    }
                }
            }
            
            //	update heights
            
            //	make sure gap between interfaces remains at least vacuumthickness
            
            if((newpos.z<0.5*my_nonbonded_params.vacuumthickness)||(newpos.z>(*pbox_dimension).z-0.5*my_nonbonded_params.vacuumthickness)){
                if(newpos.z<0.5*my_nonbonded_params.vacuumthickness){
                    if(interfaceindex==1) my_exit("upper mesh going into lower part of box!");
                    for(i=0;i<(*psolvp).meshnumber;i++) meshpositions[i].z+=(0.5*my_nonbonded_params.vacuumthickness-newpos.z);       //  shift everything
                    for(i=0;i<Nnodes;i++) coordarray[i].r.z+=(0.5*my_nonbonded_params.vacuumthickness-newpos.z);
                    (*pbox_dimension).z+=(0.5*my_nonbonded_params.vacuumthickness-newpos.z);
                }
                else if(newpos.z>(*pbox_dimension).z-0.5*my_nonbonded_params.vacuumthickness){
                    (*pbox_dimension).z+=(newpos.z-((*pbox_dimension).z-0.5*my_nonbonded_params.vacuumthickness));
                }
                flag=2;
            }
            if(flag==1){
                if(Nnodes>0){
                    updatecell_asymmetric_reverse(meshpositions[mover], &(*pmeshlink), mover);
                    updatecell(meshpositions[mover], &(*pmeshmeshlink), mover);
                }
                if(Nnodes>0){
					updateneighborlist_asymmetric_reverse(&(*pmeshlink), mover, newpos, coordarray, (*pbox_dimension));
                    updateneighborlist_positions(&(*pmeshmeshlink), mover, meshpositions, (*pbox_dimension));
                }
            }
            for(i=0;i<3;i++){
                (*psolvp).totallength[interfaceindex*3+i]=newtotallength[i];
            }
            (*psolvp).totalvolume+=volumechange;
        }
    }
    free(newtotallength);
    if(Nnodes>0){
        free(newheights);
        free(siteneighborsthatchangeheight);
        free(phenylsiteneighborsthatchangeheight);
    }
	free(newtriangledata);
    free(newbonddata);
    free(new_outerbonds_meancurvature);
    free(newoutervertexdata);
    if(Nphenylsiteneighborsthatchangeheight>0){
        free_matrix(newphenphenenergies, Nphenylsiteneighborsthatchangeheight);
    }
    return flag;
}

double calc_onebody_difference_heights(double oldheight, double newheight, onebodyparam p, solvation_parameters solvp){
    double affinityheight, interfaceheight;
    
    affinityheight=newheight-p.z0;
    interfaceheight=newheight-p.zinterface;
    
    double newaffinity=-p.solvationenergy-roomTinkcals*log(1+(exp(-p.solvationenergy/roomTinkcals)-1)/(1+exp(-4*affinityheight/solvp.interfacethickness)));
    double newinterface=-p.uinterface*exp(-0.5*pow(interfaceheight/p.sigmainterface, 2));
    
    affinityheight=oldheight-p.z0;
    interfaceheight=oldheight-p.zinterface;
    
    double oldaffinity=-p.solvationenergy-roomTinkcals*log(1+(exp(-p.solvationenergy/roomTinkcals)-1)/(1+exp(-4*affinityheight/solvp.interfacethickness)));
    double oldinterface=-p.uinterface*exp(-0.5*pow(interfaceheight/p.sigmainterface, 2));
    return (newaffinity+newinterface-oldaffinity-oldinterface);
}

void calcavmesharea_bondlengths(solvation_parameters solvp, double *avarea, double *avarea2, double *avmeancurvature, double *avmeancurvaturesquared, double *avbondlength, double *avbondlength2, double *surfacetensionenergy, double *meancurvatureenergy, double *meancurvaturesquaredenergy, double *pboundaryenergy, double *uniformbondenergy, double *pdvolume, double *pdvolumeenergy, triangledata *meshtriangledata, bonddata *meshbonddata, vertexdata *meshdata, double_triple box_dimension, double surfacepressure){
    int l, i, index;
    double meancurvaturesum, meancurvaturesquaredsum;
    declare_array(double, thisarea, 2);
    declare_array(double, thisarea2, 2);
    declare_array(double, thisbondlength, 6);
    declare_array(double, thisbondlength2, 6);
    declare_array(double, thismeancurvature, 2);
    declare_array(double, thismeancurvaturesquared, 2);
	double_triple bond;
    for(l=0;l<2;l++){
        for(i=0;i<solvp.meshnumber;i++){
            index=l*solvp.meshnumber+i;
            thisarea[l]+=meshtriangledata[index].area;
            thisarea2[l]+=pow(meshtriangledata[index].area, 2);
        }
        for(i=0;i<solvp.meshnumber*3/2;i++){
            index=solvp.meshnumber*3/2*l+i;
            thisbondlength[3*l+mod(index, 3)]+=meshbonddata[index].length;
            thisbondlength2[3*l+mod(index, 3)]+=pow(meshbonddata[index].length, 2);
            thismeancurvature[l]+=meshbonddata[index].meancurvature;
        }
        for(i=0;i<solvp.meshnumber/2;i++){
            index=solvp.meshnumber/2*l+i;
            thismeancurvaturesquared[l]+=meshdata[index].meancurvaturesquared;
        }
        for(i=0;i<3;i++){
            avbondlength[3*l+i]+=thisbondlength[3*l+i];
            avbondlength2[3*l+i]+=thisbondlength2[3*l+i];
            thisbondlength[3*l+i]/=(0.5*solvp.meshnumber);
            thisbondlength2[3*l+i]/=(0.5*solvp.meshnumber);
            uniformbondenergy[3*l+i]+=(solvp.kuniformbonds*(0.5*solvp.meshnumber)*(thisbondlength2[3*l+i]/pow(thisbondlength[3*l+i], 2)-1));
        }
        avarea[l]+=thisarea[l];
        avarea2[l]+=thisarea2[l];
        thismeancurvature[l]*=0.75;
        thismeancurvaturesquared[l]*=0.1875;
        avmeancurvature[l]+=(thismeancurvature[l]/thisarea[l]);
        avmeancurvaturesquared[l]+=(thismeancurvaturesquared[l]/thisarea[l]);
        surfacetensionenergy[l]+=(solvp.surfacetension*thisarea[l]);
        meancurvatureenergy[l]+=(-solvp.surfacetension*solvp.Tolmanlength*thismeancurvature[l]);
        meancurvaturesquaredenergy[l]+=(solvp.surfacetension*pow(solvp.meancurvaturelength, 2)*thismeancurvaturesquared[l]);
    }
    (*pboundaryenergy)+=((surfacepressure-solvp.surfacetension)*box_dimension.x*box_dimension.y);
	(*pdvolume)+=(solvp.totalvolume-solvp.targetvolume);
	(*pdvolumeenergy)+=(solvp.kvolumerestraint*pow(solvp.totalvolume-solvp.targetvolume, 2));
    free(thisarea);
    free(thisarea2);
    free(thisbondlength);
    free(thisbondlength2);
    free(thismeancurvature);
    free(thismeancurvaturesquared);
}

double_double nearest_image_3dto2d(double_triple a, double_triple b, double_triple box){
	double_double result;
    result.x=a.x-b.x;
    result.y=a.y-b.y;
    recenter(result.x, box.x);
    recenter(result.y, box.y);
    result.x+=b.x;
    result.y+=b.y;
    return result;
}

triangledata calc_triangledata(double_triple p, double_triple q, double_triple r, double_triple box_dimension){
	triangledata result;
	double a2, b2, c2;
	double_triple a, b, c;
	a=subtract_double_triple(q, p);
	b=subtract_double_triple(r, q);
	c=subtract_double_triple(p, r);
	recenter_double_triple(&a, box_dimension);
	recenter_double_triple(&b, box_dimension);
	recenter_double_triple(&c, box_dimension);
	a2=dot_product(a, a);
	b2=dot_product(b, b);
	c2=dot_product(c, c);
	result.area=sqrt(2*a2*b2+2*b2*c2+2*c2*a2-a2*a2-b2*b2-c2*c2)/4.;
	result.columnvolume=(-a.x*c.y+c.x*a.y)*(p.z+q.z+r.z)/6.;
	result.normal=cross_product(a, b);
    normalize(&(result.normal));
	return result;
}

double triangular_area(double_triple p, double_triple q, double_triple r, double_triple box_dimension){
	double a2, b2, c2;
	double_triple a, b, c;
	a=subtract_double_triple(q, p);
	b=subtract_double_triple(r, q);
	c=subtract_double_triple(p, r);
	recenter_double_triple(&a, box_dimension);
	recenter_double_triple(&b, box_dimension);
	recenter_double_triple(&c, box_dimension);
	a2=dot_product(a, a);
	b2=dot_product(b, b);
	c2=dot_product(c, c);
	return sqrt(2*a2*b2+2*b2*c2+2*c2*a2-a2*a2-b2*b2-c2*c2)/4.;
}

double projected_hexagon_area(double_triple c0, double_triple c1, double_triple c2, double_triple c3, double_triple c4, double_triple c5, double_triple box){
    double_double c1near, c2near, c3near, c4near, c5near;
    c1near=nearest_image_3dto2d(c1, c0, box);
    c2near=nearest_image_3dto2d(c2, c0, box);
    c3near=nearest_image_3dto2d(c3, c0, box);
    c4near=nearest_image_3dto2d(c4, c0, box);
    c5near=nearest_image_3dto2d(c5, c0, box);
    return 0.5*(c0.x*c1near.y+c1near.x*c2near.y+c2near.x*c3near.y+c3near.x*c4near.y+c4near.x*c5near.y+c5near.x*c0.y-c0.x*c5near.y-c5near.x*c4near.y-c4near.x*c3near.y-c3near.x*c2near.y-c2near.x*c1near.y-c1near.x*c0.y);
}

void mc_translate_single_hard(int mover, coord *coordarray, double max_translate, double_triple box_dimension, linkedlistfull *plink, bonded_params my_bonded_params, nonbonded_params my_nonbonded_params, double temperature, monomernodes **monomerid, int *chainlength, linkedlistasymmetric *pmeshlink, double_triple *meshpositions, int **mapmeshtoneighbors){
	double acceptance_prob=0, difference=0;
	coord coord_new, coord_old;
	coord_old=coord_new=coordarray[mover];
	double_triple shift=scalar_multiply_double_triple(rand_unit_cube(), (rand_int(2)*2-1)*max_translate);
	coord_new.r=add_double_triple(coordarray[mover].r, shift);
	if(my_nonbonded_params.solvationparams.interface==1) coord_new.height=height_relative_to_interface(mover, coord_new, my_nonbonded_params.solvationparams, *pmeshlink, box_dimension, meshpositions, mapmeshtoneighbors);
	difference=calc_energy_difference_changer_hard(mover, coord_old, coord_new, coordarray, box_dimension, my_bonded_params, my_nonbonded_params, monomerid, chainlength, (*plink).core);
	acceptance_prob=exp(-difference/temperature);
	if((acceptance_prob>1.0)||(rand_double<acceptance_prob)){
		fmod_double_triple(&(coord_new.r), box_dimension);
        if(my_nonbonded_params.solvationparams.interface==1) updatecell_asymmetric(coord_new.r, &(*pmeshlink), mover);
		updatecell(coord_new.r, &(*plink), mover);
		coordarray[mover]=coord_new;
        if(my_nonbonded_params.solvationparams.interface==1) updateneighborlist_asymmetric(&(*pmeshlink), mover, coord_new.r, meshpositions, box_dimension);
		updateneighborlist(&(*plink), mover, coordarray, box_dimension);
	}
}

void mc_translate_single_hard_pull(int mover, coord *coordarray, double max_translate, double_triple box_dimension, linkedlistfull *plink, bonded_params my_bonded_params, nonbonded_params my_nonbonded_params, double temperature, monomernodes **monomerid, int *chainlength, int pulledchain, double force, linkedlistasymmetric *pmeshlink, double_triple *meshpositions, int **mapmeshtoneighbors){
	double acceptance_prob, difference;
	coord coord_new, coord_old;
	coord_old=coord_new=coordarray[mover];
	double_triple shift=scalar_multiply_double_triple(rand_unit_cube(), (rand_int(2)*2-1)*max_translate);
	coord_new.r=add_double_triple(coordarray[mover].r, shift);
	if(my_nonbonded_params.solvationparams.interface==1) coord_new.height=height_relative_to_interface(mover, coord_new, my_nonbonded_params.solvationparams, (*pmeshlink), box_dimension, meshpositions, mapmeshtoneighbors);
	difference=calc_energy_difference_changer_hard(mover, coord_old, coord_new, coordarray, box_dimension, my_bonded_params, my_nonbonded_params, monomerid, chainlength, (*plink).core);
    if(coord_old.chainid==pulledchain){
        if(coord_old.monomerid==0) difference+=shift.x*force;
        else if(coord_old.monomerid==chainlength[pulledchain]-1) difference-=shift.x*force;
    }
	acceptance_prob=exp(-difference/temperature);
	if((acceptance_prob>1.0)||(rand_double<acceptance_prob)){
		fmod_double_triple(&(coord_new.r), box_dimension);
        if(my_nonbonded_params.solvationparams.interface==1) updatecell_asymmetric(coord_new.r, &(*pmeshlink), mover);
		updatecell(coord_new.r, &(*plink), mover);
		coordarray[mover]=coord_new;
        if(my_nonbonded_params.solvationparams.interface==1) updateneighborlist_asymmetric(&(*pmeshlink), mover, coord_new.r, meshpositions, box_dimension);
		updateneighborlist(&(*plink), mover, coordarray, box_dimension);
	}
}

void mc_translate_single_phenyl(int mover, coord *coordarray, double max_translate, double_triple box_dimension, linkedlistfull *plink1, linkedlistfull *plink2, bonded_params my_bonded_params, nonbonded_params my_nonbonded_params, double temperature, monomernodes **monomerid, int *chainlength, linkedlistasymmetric *pmeshlink, double_triple *meshpositions, int **mapmeshtoneighbors, int dovmmc, int *number_polymer_neighbors, int **polymer_neighbor, int *assigned, int Nchains, int max_neighbors){

    int j, target, targetpolymer, movingpolymer;
    
    double acceptance_prob, difference;
	coord coord_new, coord_old;
	coord_old=coord_new=coordarray[mover];
	double_triple shift=scalar_multiply_double_triple(rand_unit_cube(), (rand_int(2)*2-1)*max_translate);
	coord_new.r=add_double_triple(coordarray[mover].r, shift);
	if(my_nonbonded_params.solvationparams.interface==1) coord_new.height=height_relative_to_interface(mover, coord_new, my_nonbonded_params.solvationparams, (*pmeshlink), box_dimension, meshpositions, mapmeshtoneighbors);
    difference=calc_energy_difference_changer_phenyl(mover, coord_old, coord_new, coordarray, box_dimension, my_bonded_params, my_nonbonded_params, monomerid, chainlength, (*plink1).core, (*plink2).core);
    
	acceptance_prob=exp(-difference/temperature);
	if((acceptance_prob>1.0)||(rand_double<acceptance_prob)){
        
        if(dovmmc==1){
            movingpolymer=coordarray[mover].chainid;
            for(j=0;j<Nchains;j++) assigned[j]=0;
            assigned[movingpolymer]=1;
            for(j=0;j<number_polymer_neighbors[movingpolymer];j++){
                assigned[polymer_neighbor[movingpolymer][j]]=1;
            }            
        }
		fmod_double_triple(&(coord_new.r), box_dimension);
        if(my_nonbonded_params.solvationparams.interface==1) updatecell_asymmetric(coord_new.r, &(*pmeshlink), mover);
		updatecell(coord_new.r, &(*plink1), mover);
		updatecell(coord_new.r, &(*plink2), mover);
		coordarray[mover]=coord_new;
        
        if(my_nonbonded_params.solvationparams.interface==1) updateneighborlist_asymmetric(&(*pmeshlink), mover, coord_new.r, meshpositions, box_dimension);
		updateneighborlist(&(*plink1), mover, coordarray, box_dimension);
		updateneighborlist_pairenergies_phenyl(&(*plink2), mover, coordarray, box_dimension, my_nonbonded_params);
                
        if(dovmmc==1){
            for(j=0;j<(*plink2).core.number_neighbors[mover];j++){
                target=(*plink2).core.neighbor[mover][j];
                targetpolymer=coordarray[target].chainid;
                if(assigned[targetpolymer]==0){
                    
                    if((*plink2).core.pair_energy[mover][target]!=0){
                        
                        if(number_polymer_neighbors[movingpolymer]==max_neighbors) my_exit("max_neighbors not big enough (single phen self)");
                        if(number_polymer_neighbors[targetpolymer]==max_neighbors) my_exit("max_neighbors not big enough (single phen target)");
                        polymer_neighbor[movingpolymer][number_polymer_neighbors[movingpolymer]]=targetpolymer;
                        polymer_neighbor[targetpolymer][number_polymer_neighbors[targetpolymer]]=movingpolymer;
                        number_polymer_neighbors[movingpolymer]++;
                        number_polymer_neighbors[targetpolymer]++;
                        
                    }
                }
            }
        }
        
	}
}

void mc_translate_single_charged(int mover, coord *coordarray, double max_translate, double_triple box_dimension, linkedlistfull *plink1, linkedlistfull *plink2, bonded_params my_bonded_params, nonbonded_params my_nonbonded_params, double temperature, monomernodes **monomerid, int *chainlength, linkedlistasymmetric *pmeshlink, double_triple *meshpositions, int **mapmeshtoneighbors, int dovmmc, int *number_polymer_neighbors, int **polymer_neighbor, int *assigned, int Nchains, int max_neighbors){
    int j, target, targetpolymer, movingpolymer;
    double acceptance_prob, difference;
	coord coord_new, coord_old;
	coord_old=coord_new=coordarray[mover];
	double_triple shift=scalar_multiply_double_triple(rand_unit_cube(), (rand_int(2)*2-1)*max_translate);
	coord_new.r=add_double_triple(coordarray[mover].r, shift);
	if(my_nonbonded_params.solvationparams.interface==1) coord_new.height=height_relative_to_interface(mover, coord_new, my_nonbonded_params.solvationparams, (*pmeshlink), box_dimension, meshpositions, mapmeshtoneighbors);
	difference=calc_energy_difference_changer_charged(mover, coord_old, coord_new, coordarray, box_dimension, my_bonded_params, my_nonbonded_params, monomerid, chainlength, (*plink1).core, (*plink2).core);

    acceptance_prob=exp(-difference/temperature);
	if((acceptance_prob>1.0)||(rand_double<acceptance_prob)){
        if(dovmmc==1){
            movingpolymer=coordarray[mover].chainid;
            for(j=0;j<Nchains;j++) assigned[j]=0;
            assigned[movingpolymer]=1;
            for(j=0;j<number_polymer_neighbors[movingpolymer];j++){
                assigned[polymer_neighbor[movingpolymer][j]]=1;
            }
        }
        
		fmod_double_triple(&(coord_new.r), box_dimension);
        if(my_nonbonded_params.solvationparams.interface==1) updatecell_asymmetric(coord_new.r, &(*pmeshlink), mover);
		updatecell(coord_new.r, &(*plink1), mover);
		updatecell(coord_new.r, &(*plink2), mover);
		coordarray[mover]=coord_new;
        
        if(my_nonbonded_params.solvationparams.interface==1) updateneighborlist_asymmetric(&(*pmeshlink), mover, coord_new.r, meshpositions, box_dimension);
		updateneighborlist(&(*plink1), mover, coordarray, box_dimension);
		updateneighborlist_pairenergies_charged(&(*plink2), mover, coordarray, box_dimension, my_nonbonded_params);
        
        if(dovmmc==1){
            
            for(j=0;j<(*plink2).core.number_neighbors[mover];j++){
                target=(*plink2).core.neighbor[mover][j];
                targetpolymer=coordarray[target].chainid;
                if(assigned[targetpolymer]==0){
                    
                    if((*plink2).core.pair_energy[mover][target]!=0){
                        
                        if(number_polymer_neighbors[movingpolymer]==max_neighbors) my_exit("max_neighbors not big enough (single phen self)");
                        if(number_polymer_neighbors[targetpolymer]==max_neighbors) my_exit("max_neighbors not big enough (single phen target)");
                        polymer_neighbor[movingpolymer][number_polymer_neighbors[movingpolymer]]=targetpolymer;
                        polymer_neighbor[targetpolymer][number_polymer_neighbors[targetpolymer]]=movingpolymer;                        
                        number_polymer_neighbors[movingpolymer]++;
                        number_polymer_neighbors[targetpolymer]++;
                        
                    }
                }
            }
        }        
	}
}

void mc_rotate_single_hard(int mover, coord *coordarray, double max_rotate, double_triple box_dimension, linkedlist link, bonded_params my_bonded_params, nonbonded_params my_nonbonded_params, double temperature, monomernodes **monomerid, int *chainlength){
	double acceptance_prob, difference;
	coord coord_new, coord_old;
	coord_old=coord_new=coordarray[mover];
	double_triple newn=scalar_multiply_double_triple(rand_unit_cube(), max_rotate);
	newn=add_double_triple(coord_new.n, newn);
	normalize(&newn);
	coord_new.n=newn;
	difference=calc_energy_difference_changen_hard(mover, coord_old, coord_new, coordarray, box_dimension, my_bonded_params, my_nonbonded_params, monomerid, chainlength, link);
	acceptance_prob=exp(-difference/temperature);
	if((acceptance_prob>1.0)||(rand_double<acceptance_prob)){
		coordarray[mover]=coord_new;
	}
}

void mc_rotate_single_phenyl(int mover, coord *coordarray, double max_rotate, double_triple box_dimension, linkedlist link, bonded_params my_bonded_params, nonbonded_params my_nonbonded_params, double temperature, monomernodes **monomerid, int *chainlength){
	double acceptance_prob, difference;
	coord coord_new, coord_old;
	coord_old=coord_new=coordarray[mover];
	double_triple newn=scalar_multiply_double_triple(rand_unit_cube(), max_rotate);
	newn=add_double_triple(coord_new.n, newn);
	normalize(&newn);
	coord_new.n=newn;
	difference=calc_energy_difference_changen_phenyl(mover, coord_old, coord_new, coordarray, box_dimension, my_bonded_params, my_nonbonded_params, monomerid, chainlength, link);
	acceptance_prob=exp(-difference/temperature);
	if((acceptance_prob>1.0)||(rand_double<acceptance_prob)){
		coordarray[mover]=coord_new;
	}
}

void mc_rotate_single_charged(int mover, coord *coordarray, double max_rotate, double_triple box_dimension, linkedlist link, bonded_params my_bonded_params, nonbonded_params my_nonbonded_params, double temperature, monomernodes **monomerid, int *chainlength){
	double acceptance_prob, difference;
	coord coord_new, coord_old;
	coord_old=coord_new=coordarray[mover];
	double_triple newn=scalar_multiply_double_triple(rand_unit_cube(), max_rotate);
	newn=add_double_triple(coord_new.n, newn);
	normalize(&newn);
	coord_new.n=newn;
	difference=calc_energy_difference_changen_charged(mover, coord_old, coord_new, coordarray, box_dimension, my_bonded_params, my_nonbonded_params, monomerid, chainlength, link);
    acceptance_prob=exp(-difference/temperature);
	if((acceptance_prob>1.0)||(rand_double<acceptance_prob)){
		coordarray[mover]=coord_new;
	}
}

int mc_aspect_ratios_cellstruct(int Nnodes, int Nchains, coord *coordarray, coord *newcoordarray, double_triple *pbox_dimension, double maxlogaspect, bonded_params my_bonded_params, nonbonded_params *pmy_nonbonded_params, int *chainlength, monomernodes **monomerid, linkedlistset linkset, pressureparam my_pressureparam, double temperature, linkedlistasymmetric meshlink, double_triple *meshpositions, triangledata *meshtriangledata, bonddata *meshbonddata, vertexdata *meshdata, int **maptrianglestomesh, int **mapbondstomesh, int **mapmeshtoneighbors, int **mapmeshtobonds, int **mapmeshtotriangles, int **mapbondstotriangles){
	double shiftlog, factor, complement, xfactor, yfactor, zfactor, newenergy, oldenergy, acceptance_prob, surfacepressureterm, new_longestmeshbond=0, bondlength, uniformbondold=0, uniformbondnew=0, meancurvaturesum=0, meancurvatureoldsum=0, meancurvaturesquaredsum=0, meancurvaturesquaredoldsum=0, *new_totallength, surfacepressureatboundary;
	double_triple new_box_dimension, bond;
	int result=0, i, j, index, dim=rand_int(3);
    if((*pmy_nonbonded_params).solvationparams.interface==1) allocate_array(double, new_totallength, 3);
    shiftlog=(rand_double*2-1)*maxlogaspect;
	factor=exp(shiftlog);
	complement=exp(-0.5*shiftlog);
	if((*pmy_nonbonded_params).solvationparams.interface==1){
        surfacepressureatboundary=my_pressureparam.surfacepressure-(*pmy_nonbonded_params).solvationparams.surfacetension;
    }
    else{
        surfacepressureatboundary=my_pressureparam.surfacepressure;
    }
	if(my_pressureparam.pressuretype==0){			//	isotropic
		if(dim==0){
			xfactor=factor;
			yfactor=zfactor=complement;
		}
		else if(dim==1){
			yfactor=factor;
			xfactor=zfactor=complement;
		}
		else{
			zfactor=factor;
			xfactor=yfactor=complement;
		}
		surfacepressureterm=(xfactor*yfactor-1.)*(*pbox_dimension).x*(*pbox_dimension).y*surfacepressureatboundary;
	}
	else if(my_pressureparam.pressuretype==1){		//	planar (still isotropic, but only changing two sides per move)
		if(dim==0){
			xfactor=factor;
			yfactor=1./factor;
			zfactor=1.;
		}
		else if(dim==1){
			yfactor=factor;
			zfactor=1./factor;
			xfactor=1.;
		}
		else{
			zfactor=factor;
			xfactor=1./factor;
			yfactor=1.;
		}
		surfacepressureterm=(xfactor*yfactor-1.)*(*pbox_dimension).x*(*pbox_dimension).y*surfacepressureatboundary;
        surfacepressureterm+=0.5*(xfactor+1.)*0.5*(yfactor+1.)*(zfactor-1)*(*pbox_dimension).x*(*pbox_dimension).y*(*pbox_dimension).z*my_pressureparam.normalforceperunitarea;
	}
	else if(my_pressureparam.pressuretype==2){		//	uniaxial in x direction
		if(dim==0){
			xfactor=factor;
			yfactor=1./factor;
			zfactor=1.;
		}
		else if(dim==1){
			yfactor=factor;
			zfactor=1./factor;
			xfactor=1.;
		}
		else{
			zfactor=factor;
			xfactor=1./factor;
			yfactor=1.;
		}
		surfacepressureterm=(xfactor-1.)*0.5*(yfactor+1.)*(*pbox_dimension).x*(*pbox_dimension).y*surfacepressureatboundary;
	}
	else if(my_pressureparam.pressuretype==3){		//	uniaxial in y direction
		if(dim==0){
			xfactor=factor;
			yfactor=1./factor;
			zfactor=1.;
		}
		else if(dim==1){
			yfactor=factor;
			zfactor=1./factor;
			xfactor=1.;
		}
		else{
			zfactor=factor;
			xfactor=1./factor;
			yfactor=1.;
		}
		surfacepressureterm=(yfactor-1.)*0.5*(xfactor+1.)*(*pbox_dimension).x*(*pbox_dimension).y*surfacepressureatboundary;
	}
	else if(my_pressureparam.pressuretype==4){		//	fixed z, only swapping x and y
		xfactor=factor;
		yfactor=1./factor;
		zfactor=1.;
		surfacepressureterm=(xfactor*yfactor-1.)*(*pbox_dimension).x*(*pbox_dimension).y*surfacepressureatboundary;
	}
	else if(my_pressureparam.pressuretype==5){		//	fixed z, x and y fluctuating independently
		dim=rand_int(2);
		if(dim==0){
			xfactor=factor;
			yfactor=1.;
			zfactor=1.;
		}
		else if(dim==1){
			xfactor=1.;
			yfactor=factor;
			zfactor=1.;
		}
		surfacepressureterm=(xfactor*yfactor-1.)*(*pbox_dimension).x*(*pbox_dimension).y*surfacepressureatboundary;
	}
    else if(my_pressureparam.pressuretype==6){		//	both x and y, different pressures
		if(dim==0){
			xfactor=factor;
			yfactor=1./factor;
			zfactor=1.;
		}
		else if(dim==1){
			yfactor=factor;
			zfactor=1./factor;
			xfactor=1.;
		}
		else{
			zfactor=factor;
			xfactor=1./factor;
			yfactor=1.;
		}
		surfacepressureterm=(xfactor-1.)*0.5*(yfactor+1.)*(*pbox_dimension).x*(*pbox_dimension).y*surfacepressureatboundary;
		surfacepressureterm+=((yfactor-1.)*0.5*(xfactor+1.)*(*pbox_dimension).x*(*pbox_dimension).y*(my_pressureparam.ypressure-(*pmy_nonbonded_params).solvationparams.surfacetension));
	}
	else if(my_pressureparam.pressuretype==7){			//	not conserving volume; isotropic surface pressure; separate vertical pressure; changing only one dimension at a time
		if(dim==0){
			xfactor=factor;
			yfactor=zfactor=1;
		}
		else if(dim==1){
			yfactor=factor;
			xfactor=zfactor=1;
		}
		else{
			zfactor=factor;
			xfactor=yfactor=1;
		}
		surfacepressureterm=(xfactor*yfactor-1.)*(*pbox_dimension).x*(*pbox_dimension).y*surfacepressureatboundary;
        surfacepressureterm+=0.5*(xfactor+1.)*0.5*(yfactor+1.)*(zfactor-1)*(*pbox_dimension).x*(*pbox_dimension).y*(*pbox_dimension).z*my_pressureparam.normalforceperunitarea;
	}
    
	new_box_dimension.x=(*pbox_dimension).x*xfactor;
	new_box_dimension.y=(*pbox_dimension).y*yfactor;
	if((*pmy_nonbonded_params).solvationparams.interface==1){
        new_box_dimension.z=((*pbox_dimension).z-(*pmy_nonbonded_params).vacuumthickness)*zfactor+(*pmy_nonbonded_params).vacuumthickness;		//	dilate/contract only water part
    }
	else new_box_dimension.z=(*pbox_dimension).z*zfactor;
	for(i=0;i<Nnodes;i++){
		newcoordarray[i]=coordarray[i];
		newcoordarray[i].r.x*=xfactor;
		newcoordarray[i].r.y*=yfactor;
		if((*pmy_nonbonded_params).solvationparams.interface==1) newcoordarray[i].r.z=(coordarray[i].r.z-0.5*(*pbox_dimension).z)*zfactor+0.5*new_box_dimension.z;				//	dilate/contract only water part by dilating/contracting relative to center of water slab
		else newcoordarray[i].r.z*=zfactor;
		fmod_double_triple(&(newcoordarray[i].r), new_box_dimension);																											//	in case went out of bounds due to numerical drift
	}
	nonbonded_params new_nonbonded_params;
    new_nonbonded_params=(*pmy_nonbonded_params);
    double_triple *new_meshpositions;
    triangledata *new_meshtriangledata;
    bonddata *new_meshbonddata;
    vertexdata *new_meshdata;
    
    double dA=0;
	if((*pmy_nonbonded_params).solvationparams.interface==1){
        
        //  Add local surface tension term is there is an interface
        
        allocate_array(double_triple, new_meshpositions, new_nonbonded_params.solvationparams.meshnumber);
        allocate_array(triangledata, new_meshtriangledata, 2*new_nonbonded_params.solvationparams.meshnumber);
        allocate_array(bonddata, new_meshbonddata, 3*new_nonbonded_params.solvationparams.meshnumber);
        allocate_array(vertexdata, new_meshdata, 2*new_nonbonded_params.solvationparams.meshnumber);
		for(i=0;i<(*pmy_nonbonded_params).solvationparams.meshsize.x*(*pmy_nonbonded_params).solvationparams.meshsize.y;i++){
			new_meshpositions[i]=meshpositions[i];
			new_meshpositions[i].x*=xfactor;
			new_meshpositions[i].y*=yfactor;
			new_meshpositions[i].z=(meshpositions[i].z-0.5*(*pbox_dimension).z)*zfactor+0.5*new_box_dimension.z;
		}
		for(i=(*pmy_nonbonded_params).solvationparams.meshsize.x*(*pmy_nonbonded_params).solvationparams.meshsize.y;i<(*pmy_nonbonded_params).solvationparams.meshnumber;i++){
			new_meshpositions[i]=meshpositions[i];
			new_meshpositions[i].x*=xfactor;
			new_meshpositions[i].y*=yfactor;
			new_meshpositions[i].z=(meshpositions[i].z-0.5*(*pbox_dimension).z)*zfactor+0.5*new_box_dimension.z;
		}
        
        //  New triangle data for calculating new mesh energy terms, only using top interface
        
        if((new_nonbonded_params.solvationparams.meshmovecode!=1)&&(new_nonbonded_params.solvationparams.meshmovecode!=2)){
            printf("Boundary and bulk surface tension terms in mc_aspect_ratios assumes only one active interface, but mesh move code assumes two active interfaces!\n");
        }
        
        for(i=new_nonbonded_params.solvationparams.meshnumber;i<2*new_nonbonded_params.solvationparams.meshnumber;i++){
            new_meshtriangledata[i]=calc_triangledata(new_meshpositions[maptrianglestomesh[i][0]], new_meshpositions[maptrianglestomesh[i][1]], new_meshpositions[maptrianglestomesh[i][2]], (new_box_dimension));
            dA+=(new_meshtriangledata[i].area-meshtriangledata[i].area);
        }
        
        //  New bond data for calculating new uniform bond energy and curvature energy
        
        for(i=new_nonbonded_params.solvationparams.meshnumber*3/2;i<3*new_nonbonded_params.solvationparams.meshnumber;i++){
            bond=subtract_double_triple(new_meshpositions[mapbondstomesh[i][0]], new_meshpositions[mapbondstomesh[i][1]]);
            recenter_double_triple(&bond, new_box_dimension);
            bondlength=norm(bond);
            new_meshbonddata[i].length=bondlength;
            new_meshbonddata[i].meancurvature=bondlength*mysafeacos(dot_product(new_meshtriangledata[mapbondstotriangles[i][0]].normal, new_meshtriangledata[mapbondstotriangles[i][1]].normal));
            meancurvaturesum+=new_meshbonddata[i].meancurvature;
            meancurvatureoldsum+=meshbonddata[i].meancurvature;
            new_totallength[mod(i, 3)]+=bondlength;
            if(bondlength>new_longestmeshbond) new_longestmeshbond=bondlength;
        }
        
        //  New mesh data for calculating curvature energy
        
        for(i=new_nonbonded_params.solvationparams.meshnumber/2;i<new_nonbonded_params.solvationparams.meshnumber;i++){
            new_meshdata[i].area=0;
            new_meshdata[i].meancurvature=0;
            for(j=0;j<6;j++){
                new_meshdata[i].area+=new_meshtriangledata[mapmeshtotriangles[i][j]].area;
                new_meshdata[i].meancurvature+=new_meshbonddata[mapmeshtobonds[i][j]].meancurvature;
            }
            new_meshdata[i].meancurvaturesquared=pow(new_meshdata[i].meancurvature, 2)/new_meshdata[i].area;
            meancurvaturesquaredsum+=new_meshdata[i].meancurvaturesquared;
            meancurvaturesquaredoldsum+=meshdata[i].meancurvaturesquared;
        }        
        for(i=new_nonbonded_params.solvationparams.meshnumber*3/2;i<3*new_nonbonded_params.solvationparams.meshnumber;i++){
            uniformbondold+=pow((*pmy_nonbonded_params).solvationparams.meshnumber/2*meshbonddata[i].length/(*pmy_nonbonded_params).solvationparams.totallength[3+mod(i, 3)]-1, 2);
        }
        for(i=new_nonbonded_params.solvationparams.meshnumber*3/2;i<3*new_nonbonded_params.solvationparams.meshnumber;i++){
            uniformbondnew+=pow(new_nonbonded_params.solvationparams.meshnumber/2*new_meshbonddata[i].length/new_totallength[mod(i, 3)]-1, 2);
        }
        surfacepressureterm+=(*pmy_nonbonded_params).solvationparams.surfacetension*(dA-0.75*(*pmy_nonbonded_params).solvationparams.Tolmanlength*(meancurvaturesum-meancurvatureoldsum)+0.1875*pow((*pmy_nonbonded_params).solvationparams.meancurvaturelength, 2)*(meancurvaturesquaredsum-meancurvaturesquaredoldsum));
        
        for(i=0;i<Nnodes;i++){
            newcoordarray[i].height=height_relative_to_interface(i, newcoordarray[i], new_nonbonded_params.solvationparams, meshlink, new_box_dimension, new_meshpositions, mapmeshtoneighbors);
        }
    }
	newenergy=total_energy_cellstruct(newcoordarray, my_bonded_params, new_nonbonded_params, Nchains, chainlength, monomerid, new_box_dimension, linkset);			//	shouldn't need to redo cell until accepted, as long as maxlogaspect isn't too big
    if(newenergy<1.0/0.0){
		oldenergy=total_energy_cellstruct(coordarray, my_bonded_params, (*pmy_nonbonded_params), Nchains, chainlength, monomerid, *pbox_dimension, linkset);		//	shouldn't need to redo cell until accepted, as long as maxlogaspect isn't too big
		acceptance_prob=exp((oldenergy-newenergy-surfacepressureterm-(*pmy_nonbonded_params).solvationparams.kuniformbonds*(uniformbondnew-uniformbondold))/temperature);
		if((acceptance_prob>=1)||(rand_double<acceptance_prob)){
            result=1;
            (*pbox_dimension)=new_box_dimension;
            for(i=0;i<Nnodes;i++){
                coordarray[i]=newcoordarray[i];
            }
            (*pmy_nonbonded_params)=new_nonbonded_params;
            if((*pmy_nonbonded_params).solvationparams.interface==1){
                for(i=0;i<new_nonbonded_params.solvationparams.meshnumber;i++) meshpositions[i]=new_meshpositions[i];
                for(i=0;i<2*new_nonbonded_params.solvationparams.meshnumber;i++) meshtriangledata[i]=calc_triangledata(new_meshpositions[maptrianglestomesh[i][0]], new_meshpositions[maptrianglestomesh[i][1]], new_meshpositions[maptrianglestomesh[i][2]], (*pbox_dimension));
                
                for(i=new_nonbonded_params.solvationparams.meshnumber;i<2*new_nonbonded_params.solvationparams.meshnumber;i++){
                    meshtriangledata[i]=new_meshtriangledata[i];
                }
                for(i=new_nonbonded_params.solvationparams.meshnumber*3/2;i<3*new_nonbonded_params.solvationparams.meshnumber;i++){
                    meshbonddata[i]=new_meshbonddata[i];
                }
                for(i=new_nonbonded_params.solvationparams.meshnumber/2;i<new_nonbonded_params.solvationparams.meshnumber;i++){
                    meshdata[i]=new_meshdata[i];
                }
                for(i=0;i<3;i++){
                    (*pmy_nonbonded_params).solvationparams.totallength[3+i]=new_totallength[i];
                }
                (*pmy_nonbonded_params).solvationparams.longestmeshbond=new_longestmeshbond;
                if((*pmy_nonbonded_params).solvationparams.longestmeshbond>(*pmy_nonbonded_params).solvationparams.maxmeshbond) result=2;
            }
 		}
	}
    if((*pmy_nonbonded_params).solvationparams.interface==1){
        free(new_totallength);
        free(new_meshpositions);
        free(new_meshtriangledata);
        free(new_meshbonddata);
        free(new_meshdata);
    }
    if(Nnodes==0) result=0;                 //  no need to change cells if there is only a bare mesh
	return result;
}

void height_verbose_debugging(int minindex, double_triple center, double_triple *meshpositions, int *neighbors, double_triple box_dimension, double_triple mycoordr, double_triple relativepos){
    int i, sign, geometrycode, tienumber=1, allsamesign=1, found;
	double sep2, height, distance, distancealongedge, leftmag, rightleftmag;
	double_triple sep, left, right, leftdirector, rightdirector, normal, updirector, postriangleframe, nearestpoint, rightsep, leftsep, lefttriangleframe, righttriangleframe, rightleftdirector, rightleftsep, lefttriangleframeunit;
    printf("verbose debugging height\n");
    printf("site position:\t%f, %f, %f\n", mycoordr.x, mycoordr.y, mycoordr.z);
    printf("center %i:\t%f, %f, %f\n", minindex, center.x, center.y, center.z);
    printf("relative pos:\t%f, %f, %f\n", relativepos.x, relativepos.y, relativepos.z);
	declare_array(int, tieindex, 6);
	declare_array(int, tiesign, 6);
	declare_array(int, tiecode, 6);     //  0: bulk; 1 cr edge; 2 rl edge; 3 rc edge; 4 c corner; 5 r corner; 6 l corner
	declare_array(double, normaldotproductmag, 6);
	for(i=0;i<6;i++){
		right=meshpositions[neighbors[i]];
        printf("%i: right mesh point: %i: %f, %f, %f\n", i, neighbors[i], right.x, right.y, right.z);
		rightsep=subtract_double_triple(right, center);
		recenter_double_triple(&rightsep, box_dimension);
		rightdirector=rightsep;
		righttriangleframe.x=norm(rightdirector);
		rightdirector=scalar_multiply_double_triple(rightdirector, 1./righttriangleframe.x);		//	lab frame
        printf("\trighttriangleframe:\t%f, %f, %f\n", righttriangleframe.x, 0., 0.);
        printf("\trightdirector:\t%f, %f, %f\n", rightdirector.x, rightdirector.y, rightdirector.z);
        
		if(i==5) left=meshpositions[neighbors[0]];
		else left=meshpositions[neighbors[i+1]];
		leftsep=subtract_double_triple(left, center);
		recenter_double_triple(&leftsep, box_dimension);
        leftmag=norm(leftsep);
        leftdirector=scalar_multiply_double_triple(leftsep, 1./leftmag);		//	lab frame
        if(i==5) printf("\tleft mesh point: %i: %f, %f, %f\n", neighbors[0], left.x, left.y, left.z);
        else printf("\tleft mesh point: %i: %f, %f, %f\n", neighbors[i+1], left.x, left.y, left.z);
		
		normal=cross_product(rightdirector, leftdirector);
		normalize(&normal);
		updirector=cross_product(normal, rightdirector);		//	lab frame
		
		lefttriangleframe.x=dot_product(rightdirector, leftsep);
		lefttriangleframe.y=dot_product(updirector, leftsep);
		lefttriangleframe.z=0;										//	full-magnitude left vector in triangle frame
        printf("\tlefttriangleframe:\t%f, %f, %f\n", lefttriangleframe.x, lefttriangleframe.y, lefttriangleframe.z);
        printf("\tleftdirector:\t%f, %f, %f\n", leftdirector.x, leftdirector.y, leftdirector.z);
		righttriangleframe.y=0;
		righttriangleframe.z=0;
		//	calculate nearest point in triangle, in the frame of the triangle
		
		postriangleframe.x=dot_product(rightdirector, relativepos);
		postriangleframe.y=dot_product(updirector, relativepos);
		postriangleframe.z=dot_product(normal, relativepos);
        printf("\tpostriangleframe:\t%f, %f, %f\n", postriangleframe.x, postriangleframe.y, postriangleframe.z);
		
		sign=1;
		nearestpoint.z=0;
		if(postriangleframe.y<0){
			printf("\tpostriangleframe.y %f<0\n", postriangleframe.y);
			if(postriangleframe.x<0){
                printf("\t\tpostriangleframe.x<0\n");
                
                //  still need to check left edge!
                
                lefttriangleframeunit=lefttriangleframe;
                normalize(&lefttriangleframeunit);
                distancealongedge=dot_product(postriangleframe, lefttriangleframeunit);
                if(distancealongedge<0){
                    printf("\t\t\tdistance along left edge %f<0: center corner\n", distancealongedge);
                    geometrycode=4;         //  center corner
                    distance=norm(postriangleframe);
                    if(postriangleframe.z<0) sign=-1;
                    printf("\t\t\tdistance:\t%i x %.16f\n", sign, distance);
                }
                else if(distancealongedge>leftmag){
                    printf("\t\t\tdistance along left edge %f>%f: left corner\n", distancealongedge, leftmag);
                    geometrycode=6;         //  left corner
                    distance=sqrt(pow(postriangleframe.x-lefttriangleframe.x, 2)+pow(postriangleframe.y-lefttriangleframe.y, 2)+pow(postriangleframe.z, 2));		//	typo!
                    if(postriangleframe.z<0) sign=-1;
                    printf("\t\t\tdistance:\t%i x %.16f\n", sign, distance);
                }
                else{
                    printf("\t\t\tdistance along left edge %f in (0, %f): center-left edge\n", distancealongedge, leftmag);
                    geometrycode=3;         //  center-left edge
                    nearestpoint=scalar_multiply_double_triple(lefttriangleframeunit, distancealongedge);        //   do use it here
                    distance=sqrt(pow(postriangleframe.x-nearestpoint.x, 2)+pow(postriangleframe.y-nearestpoint.y, 2)+pow(postriangleframe.z, 2));
                    if(postriangleframe.z<0) sign=-1;
                    printf("\t\t\tdistance:\t%i x %.16f\n", sign, distance);
                }
			}
            else if(postriangleframe.x>righttriangleframe.x){
                printf("\t\tpostriangleframe.x %f>%f\n", postriangleframe.x, righttriangleframe.x);
                
                //  still need to check right-left edge!
                
                rightleftsep=subtract_double_triple(lefttriangleframe, righttriangleframe);
                rightleftmag=norm(rightleftsep);
                rightleftdirector=scalar_multiply_double_triple(rightleftsep, 1./rightleftmag);
                printf("\t\trightleftdirector %f, %f, %f\n", rightleftdirector.x, rightleftdirector.y, rightleftdirector.z);
                distancealongedge=dot_product(subtract_double_triple(postriangleframe, righttriangleframe), rightleftdirector);
                if(distancealongedge<0){
                    printf("\t\t\tdistance along right-left edge %f<0: right corner\n", distancealongedge);
                    geometrycode=5;         //  right corner
                    distance=sqrt(pow(postriangleframe.x-righttriangleframe.x, 2)+pow(postriangleframe.y, 2)+pow(postriangleframe.z, 2));
                    if(postriangleframe.z<0) sign=-1;
                    printf("\t\t\tdistance:\t%i x %.16f\n", sign, distance);
                }
                else if(distancealongedge>rightleftmag){
                    printf("\t\t\tdistance along right-left edge %f>%f: left corner\n", distancealongedge, rightleftmag);
                    geometrycode=6;         //  left corner
                    distance=sqrt(pow(postriangleframe.x-lefttriangleframe.x, 2)+pow(postriangleframe.y-lefttriangleframe.y, 2)+pow(postriangleframe.z, 2));
                    if(postriangleframe.z<0) sign=-1;
                    printf("\t\t\tdistance:\t%i x %.16f\n", sign, distance);
                }
                else{
                    printf("\t\t\tdistance along left edge %f in (0, %f): right-left edge\n", distancealongedge, rightleftmag);
                    geometrycode=2;         //  right-left edge
                    nearestpoint=add_double_triple(righttriangleframe, scalar_multiply_double_triple(rightleftdirector, distancealongedge));
                    distance=sqrt(pow(postriangleframe.x-nearestpoint.x, 2)+pow(postriangleframe.y-nearestpoint.y, 2)+pow(postriangleframe.z, 2));
                    if(postriangleframe.z<0) sign=-1;
                    printf("\t\t\tdistance:\t%i x %.16f\n", sign, distance);
                }
            }
			else{
                printf("\t\tpostriangleframe.x in (0, %f): center-right edge\n", righttriangleframe.x);
                geometrycode=1;         //  center-right edge
                distance=sqrt(pow(postriangleframe.y, 2)+pow(postriangleframe.z, 2));
                if(postriangleframe.z<0) sign=-1;
                printf("\t\tdistance:\t%i x %.16f\n", sign, distance);
            }
		}
        else if(lefttriangleframe.y<0) my_exit("lefttriangleframe.y<0");
        else if(postriangleframe.x<postriangleframe.y*lefttriangleframe.x/lefttriangleframe.y){		// to the left of center-left edge
            printf("\tpostriangleframe.y>0, postriangleframe.x %f < %f*%f/%f (left of center-left edge)\n", postriangleframe.x, postriangleframe.y, lefttriangleframe.x, lefttriangleframe.y);
            lefttriangleframeunit=lefttriangleframe;
            normalize(&lefttriangleframeunit);
			distancealongedge=dot_product(postriangleframe, lefttriangleframeunit);
			if(distancealongedge<0){
                printf("\t\tdistance along left edge %f < 0, center corner\n", distancealongedge);
                geometrycode=4;         //  center corner
                distance=norm(postriangleframe);
                if(postriangleframe.z<0) sign=-1;
                printf("\t\t\tdistance:\t%i x %.16f\n", sign, distance);
			}
			else if(distancealongedge>leftmag){
                printf("\t\tdistance along left edge %f > %f\n", distancealongedge, leftmag);
                
                //  still need to check right-left edge
                
                rightleftsep=subtract_double_triple(lefttriangleframe, righttriangleframe);
                rightleftmag=norm(rightleftsep);
                rightleftdirector=scalar_multiply_double_triple(rightleftsep, 1./rightleftmag);
                distancealongedge=dot_product(subtract_double_triple(postriangleframe, righttriangleframe), rightleftdirector);
                if(distancealongedge<0){
                    printf("\t\t\tdistance along right left edge %f < 0, right corner\n", distancealongedge);
                    geometrycode=5;         //  right corner
                    distance=sqrt(pow(postriangleframe.x-righttriangleframe.x, 2)+pow(postriangleframe.y, 2)+pow(postriangleframe.z, 2));
                    if(postriangleframe.z<0) sign=-1;
                    printf("\t\t\tdistance:\t%i x %.16f\n", sign, distance);
                }
                else if(distancealongedge>rightleftmag){
                    printf("\t\t\tdistance along right left edge %f > %f, left corner\n", distancealongedge, rightleftmag);
                    geometrycode=6;         //  left corner
                    distance=sqrt(pow(postriangleframe.x-lefttriangleframe.x, 2)+pow(postriangleframe.y-lefttriangleframe.y, 2)+pow(postriangleframe.z, 2));
                    if(postriangleframe.z<0) sign=-1;
                    printf("\t\t\tdistance:\t%i x %.16f\n", sign, distance);
                }
                else{
                    printf("\t\t\tdistance along right left edge %f in (0, %f), right-left edge\n", distancealongedge, rightleftmag);
                    geometrycode=2;         //  right-left edge
                    nearestpoint=add_double_triple(righttriangleframe, scalar_multiply_double_triple(rightleftdirector, distancealongedge));
                    distance=sqrt(pow(postriangleframe.x-nearestpoint.x, 2)+pow(postriangleframe.y-nearestpoint.y, 2)+pow(postriangleframe.z, 2));				
                    if(postriangleframe.z<0) sign=-1;
                    printf("\t\t\tdistance:\t%i x %.16f\n", sign, distance);
                }
			}
			else{
                printf("\t\tdistance along left edge %f in (0, %f), center-left edge\n", distancealongedge, leftmag);
                geometrycode=3;         //  center-left edge
				nearestpoint=scalar_multiply_double_triple(lefttriangleframeunit, distancealongedge);        //   do use it here
                distance=sqrt(pow(postriangleframe.x-nearestpoint.x, 2)+pow(postriangleframe.y-nearestpoint.y, 2)+pow(postriangleframe.z, 2));
                if(postriangleframe.z<0) sign=-1;
                printf("\t\t\tdistance:\t%i x %.16f\n", sign, distance);
			}
		}
        else if(postriangleframe.x>righttriangleframe.x+postriangleframe.y*(lefttriangleframe.x-righttriangleframe.x)/lefttriangleframe.y){		//	to the outside of left-right edge
            printf("\tpostriangleframe.y>0, postriangleframe.x %f > %f*%f/%f, > %f + %f*(%f-%f)/%f (inside other edges, right of right-left edge)\n", postriangleframe.x, postriangleframe.y, lefttriangleframe.x, lefttriangleframe.y, righttriangleframe.x, postriangleframe.y, lefttriangleframe.x, righttriangleframe.x, lefttriangleframe.y);
			rightleftsep=subtract_double_triple(lefttriangleframe, righttriangleframe);
            rightleftmag=norm(rightleftsep);
            rightleftdirector=scalar_multiply_double_triple(rightleftsep, 1./rightleftmag);
			distancealongedge=dot_product(subtract_double_triple(postriangleframe, righttriangleframe), rightleftdirector);
			if(distancealongedge<0){
                printf("\t\tdistance along right-left edge %f <0, right corner\n", distancealongedge);
                geometrycode=5;         //  right corner
                distance=sqrt(pow(postriangleframe.x-righttriangleframe.x, 2)+pow(postriangleframe.y, 2)+pow(postriangleframe.z, 2));
                if(postriangleframe.z<0) sign=-1;
                printf("\t\tdistance:\t%i x %.16f\n", sign, distance);
            }
			else if(distancealongedge>rightleftmag){
                printf("\t\tdistance along right-left edge %f > %f, left corner\n", distancealongedge, rightleftmag);
                geometrycode=6;         //  left corner
                distance=sqrt(pow(postriangleframe.x-lefttriangleframe.x, 2)+pow(postriangleframe.y-lefttriangleframe.y, 2)+pow(postriangleframe.z, 2));
                if(postriangleframe.z<0) sign=-1;
                printf("\t\tdistance:\t%i x %.16f\n", sign, distance);
            }
			else{
                printf("\t\tdistance along right-left edge %f in (0, %f), right-left edge\n", distancealongedge, rightleftmag);
                geometrycode=2;         //  right-left edge
                nearestpoint=add_double_triple(righttriangleframe, scalar_multiply_double_triple(rightleftdirector, distancealongedge));
                distance=sqrt(pow(postriangleframe.x-nearestpoint.x, 2)+pow(postriangleframe.y-nearestpoint.y, 2)+pow(postriangleframe.z, 2));				
                if(postriangleframe.z<0) sign=-1;
                printf("\t\tdistance:\t%i x %.16f\n", sign, distance);
            }
		}
		else{
            printf("\tinside all three, postriangleframe %f, %f, %f, direct projection\n", postriangleframe.x, postriangleframe.y, postriangleframe.z);
            geometrycode=0;         //  projects onto triangle
            distance=postriangleframe.z;
			if(postriangleframe.z<0){
				sign=-1;
				distance*=-1;
			}
            printf("\tdistance:\t%i x %.16f\n", sign, distance);
		}
		distance*=sign;
        if(i==0){
            height=distance;
			tieindex[0]=i;
			tiesign[0]=sign;
            tiecode[0]=geometrycode;
            printf("\tFIRST TRIANGLE, CLOSEST\n");
        }
		else if(fabs(distance)<fabs(height)-heightequalitytolerance){
			height=distance;
			tieindex[0]=i;
			tiesign[0]=sign;
            tiecode[0]=geometrycode;
			tienumber=1;
            printf("\tNEW CLOSEST\n");
		}
		else if(fabs(distance)<=fabs(height)+heightequalitytolerance){
			tieindex[tienumber]=i;
			tiesign[tienumber]=sign;
            tiecode[tienumber]=geometrycode;
			tienumber++;
            printf("\tTIE\n");
		}
	}
	
	declare_array(int, edgeties, 4);
	declare_array(int, nearoutsidecorner, 4);
	declare_array(int, reorder, 4);
	int *farcorner, *nextcorner, *conflictingtriangles;
    double projection, projection0, projection1;
	double_triple averagenormal, sep1, normal0, normal1;
	averagenormal.x=averagenormal.y=averagenormal.z=0;
	int exitflag=0, foundsign, firstsign, secondsign, indexdifference, sharededge, nextcornerside;
	double_double *vertices;
	if(tienumber>1){
		printf("(tieindex, tiesign, tiecode)[%i]=(%i, %i, %i)\n", 0, tieindex[0], tiesign[0], tiecode[0]);
		for(i=1;i<tienumber;i++){
			printf("(tieindex, tiesign, tiecode)[%i]=(%i, %i, %i)\n", i, tieindex[i], tiesign[i], tiecode[i]);
			if(tiesign[i]!=tiesign[0]){
				allsamesign=0;
			}
		}
		if(allsamesign==0){
			printf("not all same sign\n");
			found=0;
			allsamesign=1;
			for(i=0;i<tienumber;i++){
				if(tiecode[i]==0){
					if(found==2){
						printf("found three direct projections\n");
					}
					edgeties[found]=i;								//	using edgeties for direct projections; will reset with found=0
                    printf("direct tie: edgeties[%i]=%i\n", found, i);
					if(tiesign[i]!=tiesign[edgeties[0]]) allsamesign=0;
					found++;
				}
			}
			if(found==1){
                printf("one and only one direct projection");
				height=height/tiesign[0]*tiesign[i];                //  if one and only one direct projection, use that sign
			}
			else if(found==2){
				if(allsamesign==1){
                    printf("both direct projections have same sign\n");
					height=height/tiesign[0]*tiesign[edgeties[0]];      //  if two direct projections with the same sign, use it
				}
				else{
					allocate_array(int, farcorner, 2);
					allocate_array(int, nextcorner, 2);
					indexdifference=mod(edgeties[1]-edgeties[0], 6);
					if((indexdifference!=1)&&(indexdifference!=5)){
						printf("direct projection onto two non-adjacent triangles, %i and %i\n", edgeties[0], edgeties[1]);
					}
					if(indexdifference==1){
						sharededge=edgeties[1];
						farcorner[1]=mod(edgeties[1]+1, 6);
						farcorner[0]=edgeties[0];
						nextcorner[1]=mod(edgeties[1]+2, 6);
						nextcorner[0]=mod(edgeties[0]-1, 6);
					}
					else{
						sharededge=edgeties[0];
						farcorner[0]=mod(edgeties[0]+1, 6);
						farcorner[1]=edgeties[1];
						nextcorner[0]=mod(edgeties[0]+2, 6);
						nextcorner[1]=mod(edgeties[1]-1, 6);
					}
                    printf("sharededge=%i, far corners %i, %i, next corners %i, %i\n", sharededge, farcorner[0], farcorner[1], nextcorner[0], nextcorner[1]);
					
					//	recalculate orthonormal basis for edgeties[0] (arbitrary choice)
					
					right=meshpositions[neighbors[edgeties[0]]];
					rightsep=subtract_double_triple(right, center);
					recenter_double_triple(&rightsep, box_dimension);
					rightdirector=rightsep;
					righttriangleframe.x=norm(rightdirector);
					rightdirector=scalar_multiply_double_triple(rightdirector, 1./righttriangleframe.x);		//	lab frame
					left=meshpositions[neighbors[mod(edgeties[0]+1, 6)]];
					leftsep=subtract_double_triple(left, center);
					recenter_double_triple(&leftsep, box_dimension);
					leftmag=norm(leftsep);
					leftdirector=scalar_multiply_double_triple(leftsep, 1./leftmag);		//	lab frame
					normal=cross_product(rightdirector, leftdirector);
					normalize(&normal);
					updirector=cross_product(normal, rightdirector);		//	lab frame
                    printf("rightdirector = %f, %f, %f\n", rightdirector.x, rightdirector.y, rightdirector.z);
                    printf("updirector = %f, %f, %f\n", updirector.x, updirector.y, updirector.z);
					
					allocate_array(double_double, vertices, 3);		//	center, farcorner[0], shared corner, farcorner[1]; no handedness implied
					vertices[0].x=0;
					vertices[0].y=0;
					sep=subtract_double_triple(meshpositions[neighbors[farcorner[0]]], center);
					recenter_double_triple(&sep, box_dimension);
                    printf("sep=%f, %f, %f\n", sep.x, sep.y, sep.z);
					vertices[1].x=dot_product(sep, rightdirector);
					vertices[1].y=dot_product(sep, updirector);
                    printf("vertices[1]=%f, %f\n", vertices[1].x, vertices[1].y);
					sep=subtract_double_triple(meshpositions[neighbors[sharededge]], center);
					recenter_double_triple(&sep, box_dimension);
                    printf("sep=%f, %f, %f\n", sep.x, sep.y, sep.z);
					vertices[2].x=dot_product(sep, rightdirector);
					vertices[2].y=dot_product(sep, updirector);
                    printf("vertices[2]=%f, %f\n", vertices[2].x, vertices[2].y);
					sep=subtract_double_triple(meshpositions[neighbors[farcorner[1]]], center);
					recenter_double_triple(&sep, box_dimension);
                    printf("sep=%f, %f, %f\n", sep.x, sep.y, sep.z);
					vertices[3].x=dot_product(sep, rightdirector);
					vertices[3].y=dot_product(sep, updirector);
                    printf("vertices[3]=%f, %f\n", vertices[3].x, vertices[3].y);
					
					if(checkpointintrianglelist(vertices[1], vertices[0], vertices[2], vertices[3])==1){		//	farcorner 0 inside triangle 1
                        printf("far corner 0 inside triangle 1\n");
						sep=subtract_double_triple(meshpositions[neighbors[nextcorner[0]]], center);		//	which side is NEXT corner on?
                        recenter_double_triple(&sep, box_dimension);
                        printf("sep=%f, %f, %f\n", sep.x, sep.y, sep.z);
                        printf("normal=%f, %f, %f\n", normal.x, normal.y, normal.z);
						if(dot_product(sep, normal)>0) nextcornerside=1;
						else nextcornerside=-1;
                        printf("nextcorner side %i\n", nextcornerside);
						if(tiesign[edgeties[0]]==nextcornerside){
                            printf("same side, keep sign\n");
							height=height/tiesign[0]*tiesign[edgeties[0]];										//	farcorner on same side: keep sign
						}
						else{
							height=-height/tiesign[0]*tiesign[edgeties[0]];										//	farcorner on different side: reverse sign
                            printf("different side, reverse sign\n");
						}
					}
					else if(checkpointintrianglelist(vertices[3], vertices[0], vertices[1], vertices[2])==1){		//	farcorner 1 inside triangle 0
                        printf("far corner 1 inside triangle 0\n");
						sep=subtract_double_triple(meshpositions[neighbors[nextcorner[1]]], center);		//	which side is NEXT corner on?
                        recenter_double_triple(&sep, box_dimension);
                        printf("sep=%f, %f, %f\n", sep.x, sep.y, sep.z);
                        printf("normal=%f, %f, %f\n", normal.x, normal.y, normal.z);
						if(dot_product(sep, normal)>0) nextcornerside=1;
						else nextcornerside=-1;
                        printf("nextcorner side %i\n", nextcornerside);
						if(tiesign[edgeties[1]]==nextcornerside){
							height=height/tiesign[0]*tiesign[edgeties[1]];										//	farcorner on same side: keep sign
                            printf("same side, keep sign\n");
						}
						else{
							height=-height/tiesign[0]*tiesign[edgeties[1]];										//	farcorner on different side: reverse sign
                            printf("different side, reverse sign\n");
						}
					}
					else{
                        sep=subtract_double_triple(meshpositions[neighbors[nextcorner[0]]], center);
                        recenter_double_triple(&sep, box_dimension);
                        sep1=subtract_double_triple(meshpositions[neighbors[nextcorner[1]]], center);
                        recenter_double_triple(&sep1, box_dimension);
                        printf("sep=%f, %f, %f\n", sep.x, sep.y, sep.z);
                        printf("sep1=%f, %f, %f\n", sep1.x, sep1.y, sep1.z);
                        printf("normal=%f, %f, %f\n", normal.x, normal.y, normal.z);
                        projection0=dot_product(sep, normal);
                        projection1=dot_product(sep1, normal);
                        printf("projections %f, %f\n", projection0, projection1);
                        if(projection*projection1<0){               //  different sides; no discrepancy
                            if(projection0>0) nextcornerside=1;
                            else nextcornerside=-1;
                            printf("different sides; no discrepancy; nextcornerside=%i, tiesign[%i]=%i\n", nextcornerside, edgeties[0], tiesign[edgeties[0]]);
                            if(tiesign[edgeties[0]]==nextcornerside){
                                height=height/tiesign[0]*tiesign[edgeties[0]];										//	farcorner on same side: keep sign
                            }
                            else{
                                height=-height/tiesign[0]*tiesign[edgeties[0]];										//	farcorner on different side: reverse sign
                            }
                        }
                        else{
                            
                            //  use next corner whose triangle is more planar with original triangle; either this is the closest or two next triangles must overlap (not allowed)
                            
                            printf("check which triangle more planar\n");
                            sep=subtract_double_triple(meshpositions[neighbors[farcorner[0]]], center);
                            recenter_double_triple(&sep, box_dimension);
                            sep1=subtract_double_triple(meshpositions[neighbors[nextcorner[0]]], center);
                            recenter_double_triple(&sep, box_dimension);
                            normal0=cross_product(sep, sep1);
                            normalize(&normal0);
                            printf("normal0=%f, %f, %f\n", normal0.x, normal0.y, normal0.z);
                            sep=subtract_double_triple(meshpositions[neighbors[farcorner[1]]], center);
                            recenter_double_triple(&sep, box_dimension);
                            sep1=subtract_double_triple(meshpositions[neighbors[nextcorner[1]]], center);
                            recenter_double_triple(&sep, box_dimension);
                            normal1=cross_product(sep, sep1);
                            normalize(&normal1);
                            printf("normal1=%f, %f, %f\n", normal1.x, normal1.y, normal1.z);
                            if(fabs(dot_product(normal0, normal))>fabs(dot_product(normal1, normal))){
                                height=height/tiesign[0]*tiesign[edgeties[0]];										//	farcorner on same side: keep sign
                            }
                            else{
                                height=-height/tiesign[0]*tiesign[edgeties[0]];										//	farcorner on same side: keep sign
                            }
                        }
					}
					free(farcorner);
					free(nextcorner);
				}
			}
			else if(found>2){
				printf("found %i direct projections!\n", found);
				exit(1);
			}
			else{
				found=0;
				allsamesign=1;
				for(i=0;i<tienumber;i++){
					if(tiecode[i]<4){                                   //  looking for edges
						if(found==4){
							my_exit("site equally close to three edges");
						}
						edgeties[found]=i;
						if(tiesign[i]!=tiesign[edgeties[0]]) allsamesign=0;
						found++;
					}
				}
				if((found>0)&&(allsamesign==1)){
					printf("all edge projections have the same sign\n");
					height=height/tiesign[0]*tiesign[edgeties[0]];      //  if at least one edge site and all the same sign, use it
				}
				else{
					printf("no direct projections, edge projections with both signs\n");
					
					
					if(found==1){
						if(tiecode[edgeties[0]]!=2){
							my_exit("found only one close interior edge; should be an equivalent one");
						}
                        printf("one and only one edge projection\n");
						height=height/tiesign[0]*tiesign[edgeties[0]];      //  if one and only edge projection, use that sign
					}
					else if(found==3){
						for(i=0;i<2;i++){
							if(tiecode[edgeties[i]]==2){
								my_exit("two equally close edge projections, and one is outside edge");     //  if one outside edge, exit, since this means they are different edges
							}
							else if(tiecode[edgeties[i]]==1){                   //  right edge
								nearoutsidecorner[i]=tieindex[edgeties[i]];
							}
							else if(tiecode[edgeties[i]]==3){                   //  left edge
								if(tieindex[edgeties[i]]==5) nearoutsidecorner[i]=0;
								else nearoutsidecorner[i]=tieindex[edgeties[i]]+1;
							}
							else my_exit("far corner not defined");
						}
						found=0;
						if(nearoutsidecorner[1]==nearoutsidecorner[0]){
							reorder[0]=0;
							reorder[1]=1;
							found=2;
						}
						else if(nearoutsidecorner[2]==nearoutsidecorner[0]){
							reorder[0]=0;
							reorder[1]=2;
							found=2;
						}
						else if(nearoutsidecorner[2]==nearoutsidecorner[1]){
							reorder[0]=1;
							reorder[1]=2;
							found=2;
						}
						if(found==2){
							if(tiesign[edgeties[reorder[0]]]==tiesign[edgeties[reorder[1]]]){
								height=height/tiesign[0]*tiesign[edgeties[reorder[0]]];      //  two equally close edge projections, but same sign
							}							
							else if(normaldotproductmag[tieindex[edgeties[reorder[0]]]]>normaldotproductmag[tieindex[edgeties[reorder[1]]]]){
								height=height/tiesign[0]*tiesign[edgeties[reorder[0]]];
							}
							else height=height/tiesign[0]*tiesign[edgeties[reorder[1]]];
						}
						else{
							if(normaldotproductmag[tieindex[edgeties[0]]]>normaldotproductmag[tieindex[edgeties[1]]]){
								if(normaldotproductmag[tieindex[edgeties[0]]]>normaldotproductmag[tieindex[edgeties[2]]]){
									height=height/tiesign[0]*tiesign[edgeties[0]];
								}
								else height=height/tiesign[0]*tiesign[edgeties[2]];
							}
							else{
								if(normaldotproductmag[tieindex[edgeties[1]]]>normaldotproductmag[tieindex[edgeties[2]]]){
									height=height/tiesign[0]*tiesign[edgeties[1]];
								}
								else height=height/tiesign[0]*tiesign[edgeties[2]];
							}
						}
					}
					else if(found==2){                                          //  one pair of edges
						if(tiesign[edgeties[0]]==tiesign[edgeties[1]]){
							height=height/tiesign[0]*tiesign[edgeties[0]];      //  two equally close edge projections, but same sign
						}
						else{
							for(i=0;i<2;i++){
								if(tiecode[edgeties[i]]==2){
									my_exit("two equally close edge projections, and one is outside edge");     //  if one outside edge, exit, since this means they are different edges
								}
								else if(tiecode[edgeties[i]]==1){                   //  right edge
									nearoutsidecorner[i]=tieindex[edgeties[i]];
								}
								else if(tiecode[edgeties[i]]==3){                   //  left edge
									if(tieindex[edgeties[i]]==5) nearoutsidecorner[i]=0;
									else nearoutsidecorner[i]=tieindex[edgeties[i]]+1;
								}
								else my_exit("far corner not defined");
							}
							if(nearoutsidecorner[0]!=nearoutsidecorner[1]){
								printf("two equally close but different edges\n");
							}
							if(normaldotproductmag[tieindex[edgeties[0]]]>normaldotproductmag[tieindex[edgeties[1]]]){
								height=height/tiesign[0]*tiesign[edgeties[0]];
							}
							else height=height/tiesign[0]*tiesign[edgeties[1]];
						}
					}
					else if(found==4){                                      //  two pairs of edges
						for(i=0;i<4;i++){
							if(tiecode[edgeties[i]]==2){
								my_exit("four equally close edge projections, and one is outside edge");     //  if one outside edge, exit, since this means they are different edges
							}
							else if(tiecode[edgeties[i]]==1){                   //  right edge
								nearoutsidecorner[i]=tieindex[edgeties[i]];
							}
							else if(tiecode[edgeties[i]]==3){                   //  left edge
								if(tieindex[edgeties[i]]==5) nearoutsidecorner[i]=0;
								else nearoutsidecorner[i]=tieindex[edgeties[i]]+1;
							}
							else my_exit("far corner not defined");
						}
						reorder[0]=0;
						if(nearoutsidecorner[1]==nearoutsidecorner[0]){
							reorder[1]=1;
							reorder[2]=2;
							if(nearoutsidecorner[3]==nearoutsidecorner[2]) reorder[3]=3;
							else{
								my_exit("four edge corners don't match up");
							}
						}
						else if(nearoutsidecorner[2]==nearoutsidecorner[0]){
							reorder[1]=2;
							reorder[2]=1;
							if(nearoutsidecorner[3]==nearoutsidecorner[1]) reorder[3]=3;
							else{
								my_exit("four edge corners don't match up");
							}
						}
						else if(nearoutsidecorner[3]==nearoutsidecorner[0]){
							reorder[1]=3;
							reorder[2]=1;
							if(nearoutsidecorner[2]==nearoutsidecorner[1]) reorder[3]=2;
							else{
								my_exit("four edge corners don't match up");
							}
						}
						else{
							my_exit("four edge corners don't match up");
						}
						
						if(normaldotproductmag[tieindex[edgeties[reorder[0]]]]>normaldotproductmag[tieindex[edgeties[reorder[1]]]]){
							firstsign=tiesign[edgeties[reorder[0]]];
						}
						else firstsign=tiesign[edgeties[reorder[1]]];
						if(normaldotproductmag[tieindex[edgeties[reorder[2]]]]>normaldotproductmag[tieindex[edgeties[reorder[3]]]]){
							secondsign=tiesign[edgeties[reorder[2]]];
						}
						else secondsign=tiesign[edgeties[reorder[3]]];
						
						if(firstsign==secondsign){
							height=height/tiesign[0]*firstsign;
						}
						else{
                            							
							//	two pairs of edges with different signs
							//	if they are separated by one, assume that two intervening triangles are folded
							//	break tie as done for two direct projections
                            
							allocate_array(int, farcorner, 2);
							allocate_array(int, nextcorner, 2);
							allocate_array(int, conflictingtriangles, 2);		//	here, use conflictingtriangles NOT edgeties to index the two triangles
							if(tiecode[edgeties[reorder[0]]]==1){				//	right
								conflictingtriangles[0]=tieindex[edgeties[reorder[0]]];
							}
							else if(tiecode[edgeties[reorder[0]]]==3){
								conflictingtriangles[0]=mod(tieindex[edgeties[reorder[0]]]+1, 6);
							}
							else{
								printf("tiecode[edgeties[reorder[%i]]] = tieindex[%i] = %i not 1 or 3\n", 0, edgeties[reorder[0]], tiecode[edgeties[reorder[0]]]);
								exit(1);
							}
							if(tiecode[edgeties[reorder[1]]]==1){				//	right
								if(conflictingtriangles[0]!=tieindex[edgeties[reorder[1]]]){
									printf("tieindex[edgeties[reorder[%i]]] = tieindex[%i] = %i not matching\n", 1, edgeties[reorder[1]], tieindex[edgeties[reorder[1]]]);
									exit(1);
								}
							}
							else if(tiecode[edgeties[reorder[1]]]==3){
                                if(conflictingtriangles[0]!=mod(tieindex[edgeties[reorder[1]]]+1, 6)){
									printf("(tieindex[edgeties[reorder[%i]]] + 1) mod 6 = (tieindex[%i] + 1) mod 6 = %i not matching\n", 1, edgeties[reorder[1]], mod(tieindex[edgeties[reorder[1]]]+1, 6));
									exit(1);
								}
							}
							else{
								printf("tiecode[edgeties[reorder[%i]]] = tieindex[%i] = %i not 1 or 3\n", 1, edgeties[reorder[1]], tiecode[edgeties[reorder[1]]]);
								exit(1);
							}
							
							if(tiecode[edgeties[reorder[2]]]==1){				//	right
								conflictingtriangles[1]=tieindex[edgeties[reorder[2]]];
							}
							else if(tiecode[edgeties[reorder[2]]]==3){
								conflictingtriangles[1]=mod(tieindex[edgeties[reorder[2]]]+1, 6);
							}
							else{
								printf("tiecode[edgeties[reorder[%i]]] = tieindex[%i] = %i not 1 or 3\n", 0, edgeties[reorder[2]], tiecode[edgeties[reorder[2]]]);
								exit(1);
							}
							if(tiecode[edgeties[reorder[3]]]==1){				//	right
								if(conflictingtriangles[1]!=tieindex[edgeties[reorder[3]]]){
									printf("tieindex[edgeties[reorder[%i]]] = tieindex[%i] = %i not matching\n", 1, edgeties[reorder[3]], tieindex[edgeties[reorder[3]]]);
									exit(1);
								}
							}
							else if(tiecode[edgeties[reorder[3]]]==3){
                                if(conflictingtriangles[1]!=mod(tieindex[edgeties[reorder[3]]]+1, 6)){
									printf("(tieindex[edgeties[reorder[%i]]] + 1) mod 6 = (tieindex[%i] + 1) mod 6 = %i not matching\n", 1, edgeties[reorder[3]], mod(tieindex[edgeties[reorder[3]]]+1, 6));
									exit(1);
								}
							}
							else{
								printf("tiecode[edgeties[reorder[%i]]] = tieindex[%i] = %i not 1 or 3\n", 1, edgeties[reorder[3]], tiecode[edgeties[reorder[3]]]);
								exit(1);
							}
							
							indexdifference=mod(conflictingtriangles[1]-conflictingtriangles[0], 6);
                            if((indexdifference!=2)&&(indexdifference!=4)){
								printf("edge projection onto different edges %i and %i not separated by 2\n", conflictingtriangles[0], conflictingtriangles[1]);
							}
							if(indexdifference==2){
								sharededge=mod(conflictingtriangles[0]+1, 6);
								farcorner[1]=conflictingtriangles[1];
								farcorner[0]=conflictingtriangles[0];
								nextcorner[1]=mod(conflictingtriangles[1]+1, 6);
								nextcorner[0]=mod(conflictingtriangles[0]-1, 6);
							}
							else{
								sharededge=mod(conflictingtriangles[0]-1, 6);
								farcorner[0]=conflictingtriangles[0];
								farcorner[1]=conflictingtriangles[1];
								nextcorner[0]=mod(conflictingtriangles[0]+1, 6);
								nextcorner[1]=mod(conflictingtriangles[1]-1, 6);
							}
							
							//	recalculate orthonormal basis for conflictingtriangles[0] (arbitrary choice)
							
							right=meshpositions[neighbors[conflictingtriangles[0]]];
							rightsep=subtract_double_triple(right, center);
							recenter_double_triple(&rightsep, box_dimension);
							rightdirector=rightsep;
							righttriangleframe.x=norm(rightdirector);
							rightdirector=scalar_multiply_double_triple(rightdirector, 1./righttriangleframe.x);		//	lab frame
							left=meshpositions[neighbors[mod(conflictingtriangles[0]+1, 6)]];
							leftsep=subtract_double_triple(left, center);
							recenter_double_triple(&leftsep, box_dimension);
							leftmag=norm(leftsep);
							leftdirector=scalar_multiply_double_triple(leftsep, 1./leftmag);		//	lab frame
							normal=cross_product(rightdirector, leftdirector);
							normalize(&normal);
							updirector=cross_product(normal, rightdirector);		//	lab frame
							
							allocate_array(double_double, vertices, 4);		//	center, farcorner[0], shared corner, farcorner[1]; no handedness implied
							vertices[0].x=0;
							vertices[0].y=0;
							sep=subtract_double_triple(meshpositions[neighbors[farcorner[0]]], center);
							recenter_double_triple(&sep, box_dimension);
							vertices[1].x=dot_product(sep, rightdirector);
							vertices[1].y=dot_product(sep, updirector);
							sep=subtract_double_triple(meshpositions[neighbors[sharededge]], center);
							recenter_double_triple(&sep, box_dimension);
							vertices[2].x=dot_product(sep, rightdirector);
							vertices[2].y=dot_product(sep, updirector);
							sep=subtract_double_triple(meshpositions[neighbors[farcorner[1]]], center);
							recenter_double_triple(&sep, box_dimension);
							vertices[3].x=dot_product(sep, rightdirector);
							vertices[3].y=dot_product(sep, updirector);
							
							if(checkpointintrianglelist(vertices[1], vertices[0], vertices[2], vertices[3])==1){		//	farcorner 0 inside triangle 1
								sep=subtract_double_triple(meshpositions[neighbors[nextcorner[0]]], center);		//	which side is NEXT corner on?
								recenter_double_triple(&sep, box_dimension);
								if(dot_product(sep, normal)>0) nextcornerside=1;
								else nextcornerside=-1;
								if(tiesign[edgeties[reorder[0]]]==nextcornerside){
									height=height/tiesign[0]*tiesign[edgeties[reorder[0]]];										//	farcorner on same side: keep sign
								}
								else{
									height=-height/tiesign[0]*tiesign[edgeties[reorder[0]]];										//	farcorner on different side: reverse sign
								}
							}
							else if(checkpointintrianglelist(vertices[3], vertices[0], vertices[1], vertices[2])==1){		//	farcorner 1 inside triangle 0
								sep=subtract_double_triple(meshpositions[neighbors[nextcorner[1]]], center);		//	which side is NEXT corner on?
								recenter_double_triple(&sep, box_dimension);
								if(dot_product(sep, normal)>0) nextcornerside=1;
								else nextcornerside=-1;
								if(tiesign[edgeties[reorder[2]]]==nextcornerside){
									height=height/tiesign[0]*tiesign[edgeties[reorder[2]]];										//	farcorner on same side: keep sign
								}
								else{
									height=-height/tiesign[0]*tiesign[edgeties[reorder[2]]];										//	farcorner on different side: reverse sign
								}
							}
							else{
								sep=subtract_double_triple(meshpositions[neighbors[nextcorner[0]]], center);
								recenter_double_triple(&sep, box_dimension);
								sep1=subtract_double_triple(meshpositions[neighbors[nextcorner[1]]], center);
								recenter_double_triple(&sep1, box_dimension);
								projection0=dot_product(sep, normal);
								projection1=dot_product(sep1, normal);
								if(projection*projection1<0){               //  different sides; no discrepancy
									if(projection0>0) nextcornerside=1;
									else nextcornerside=-1;
									if(tiesign[edgeties[reorder[0]]]==nextcornerside){
										height=height/tiesign[0]*tiesign[edgeties[reorder[0]]];										//	farcorner on same side: keep sign
									}
									else{
										height=-height/tiesign[0]*tiesign[edgeties[reorder[0]]];										//	farcorner on different side: reverse sign
									}
								}
								else{
									
									//  use next corner whose triangle is more planar with original triangle; either this is the closest or two next triangles must overlap (not allowed)
									
									sep=subtract_double_triple(meshpositions[neighbors[farcorner[0]]], center);
									recenter_double_triple(&sep, box_dimension);
									sep1=subtract_double_triple(meshpositions[neighbors[nextcorner[0]]], center);
									recenter_double_triple(&sep, box_dimension);
									normal0=cross_product(sep, sep1);
									normalize(&normal0);
									sep=subtract_double_triple(meshpositions[neighbors[farcorner[1]]], center);
									recenter_double_triple(&sep, box_dimension);
									sep1=subtract_double_triple(meshpositions[neighbors[nextcorner[1]]], center);
									recenter_double_triple(&sep, box_dimension);
									normal1=cross_product(sep, sep1);
									normalize(&normal1);
									if(fabs(dot_product(normal0, normal))>fabs(dot_product(normal1, normal))){
										height=height/tiesign[0]*tiesign[edgeties[reorder[0]]];										//	farcorner on same side: keep sign
									}
									else{
										height=-height/tiesign[0]*tiesign[edgeties[reorder[0]]];										//	farcorner on same side: keep sign
									}
								}
							}
							free(farcorner);
							free(nextcorner);
							free(conflictingtriangles);
							
						}
					}
					else if(found>4){
						printf("found %i edge projections!\n", found);
						exit(1);
					}
					else{           //  all corner projections
						for(i=0;i<tienumber;i++){
							if(tiecode[i]!=4){
								//my_exit("all corner projections but different corners");
								exitflag=1;
							}
							if(exitflag==1){
								printf("all corner projections but different corners; tienumber=%i; debugging...", tienumber);
								for(i=0;i<tienumber;i++){
									printf("\t%i: %i: code %i\n", i, tieindex[i], tiecode[i]);
								}
								exit(1);
							}
							
							//  recalculate normal
							
							printf("not recalculating normal in height_verbose_debug\n");
							exit(1);
						}
						if(dot_product(relativepos, averagenormal)<0){
							height=height/tiesign[0]*(-1);
						}
						else height=height/tiesign[0];
					}
				}
			}
		}
	}
	free(nearoutsidecorner);
	free(edgeties);
	free(reorder);
	free(tieindex);
	free(tiesign);
	free(tiecode);
	free(normaldotproductmag);
	
	printf("ending here in height_verbose_debug\n\n");
}

void cornersdontmatchdebug(int *edgeties, int *tieindex, int *nearoutsidecorner, int *tiecode){
	int i;
	printf("four edge corners don't match up, edges:\n");
	for(i=0;i<4;i++){
		printf("\t%i: index %i, code %i, nearoutsidecorner %i\n", edgeties[i], tieindex[edgeties[i]], tiecode[edgeties[i]], nearoutsidecorner[i]);
	}
	printf("this is how they got these values of nearoutsidecorner:\n");
	for(i=0;i<4;i++){
		printf("\tedgeties[%i]=%i, code %i\n", i, edgeties[i], tiecode[edgeties[i]]);
		if(tiecode[edgeties[i]]==1){                   //  right edge
			printf("\t\tright edge, tieindex[%i]=%i\n", edgeties[i], tieindex[edgeties[i]]);
			nearoutsidecorner[i]=tieindex[edgeties[i]];
			printf("\t\t\tnearoutsidecorner[%i]=%i\n", i, nearoutsidecorner[i]);
		}
		else if(tiecode[edgeties[i]]==3){                   //  left edge
			printf("\t\tleft edge, tieindex[%i]=%i\n", edgeties[i], tieindex[edgeties[i]]);
			if(tieindex[edgeties[i]]==5){
				nearoutsidecorner[i]=0;
				printf("\t\t\ttieindex[%i]=5: nearoutsidecorner[%i]=%i\n", edgeties[i], i, nearoutsidecorner[i]);
			}
			else{
				nearoutsidecorner[i]=tieindex[edgeties[i]]+1;
				printf("\t\t\ttieindex[%i] not 5: nearoutsidecorner[%i]=%i\n", edgeties[i], i, nearoutsidecorner[i]);
			}
		}
	}
}


double height_relative_to_interface(int n, coord mycoord, solvation_parameters solvationparams, linkedlistasymmetric meshlink, double_triple box_dimension, double_triple *meshpositions, int **mapmeshtoneighbors){
	int i, j, minindex=-1, closesttriangle, leftneighbor, tienumber=1, sign, allsamesign=1, geometrycode, found;
	double sep2, height, distance, distancealongedge, leftmag, rightleftmag;
	double mindistance2=pow(solvationparams.meshsurfacecutoff, 2);
	double_triple sep, center, left, right, leftdirector, rightdirector, normal, relativepos, updirector, postriangleframe, nearestpoint, rightsep, leftsep, lefttriangleframe, righttriangleframe, rightleftdirector, rightleftsep, lefttriangleframeunit;
	
	//  find index of nearest meshpoint, minindex, and separation, sep=coord-meshpoint
	
	if(meshlink.core.number_neighbors[n]>0){
		for(i=0;i<meshlink.core.number_neighbors[n];i++){
			sep=subtract_double_triple(mycoord.r, meshpositions[meshlink.core.neighbor[n][i]]);
			recenter_double_triple(&sep, box_dimension);
			sep2=dot_product(sep, sep);
			if(sep2<mindistance2){
				mindistance2=sep2;
				minindex=meshlink.core.neighbor[n][i];
				relativepos=sep;
			}
		}
	}
	if(minindex==-1){						//	assign height=meshsurfacecutoff
		
		height=-solvationparams.meshsurfacecutoff;              //  ASSUME NEGATIVE, BECAUSE POSITIVE VALUE IS EXTREMELY (charged, backbone) OR VERY (phenyl) UNLIKELY
		return height;
		
	}
    
	//	project onto each of six triangles; calculate minimum distance to each triangle; keep smallest (following D. A. Simon's Ph. D. dissertation)
	
	declare_array(int, tieindex, 6);
	declare_array(int, tiesign, 6);
	declare_array(int, tiecode, 6);     //  0: bulk; 1 cr edge; 2 rl edge; 3 rc edge; 4 c corner; 5 r corner; 6 l corner
	declare_array(double, normaldotproductmag, 6);
	center=meshpositions[minindex];
	for(i=0;i<6;i++){
		right=meshpositions[mapmeshtoneighbors[minindex][i]];
		rightsep=subtract_double_triple(right, center);
		recenter_double_triple(&rightsep, box_dimension);
		rightdirector=rightsep;
		righttriangleframe.x=norm(rightdirector);
		rightdirector=scalar_multiply_double_triple(rightdirector, 1./righttriangleframe.x);		//	lab frame
		
		if(i==5) left=meshpositions[mapmeshtoneighbors[minindex][0]];
		else left=meshpositions[mapmeshtoneighbors[minindex][i+1]];
		leftsep=subtract_double_triple(left, center);
		recenter_double_triple(&leftsep, box_dimension);
		leftmag=norm(leftsep);
		leftdirector=scalar_multiply_double_triple(leftsep, 1./leftmag);		//	lab frame
		
		normal=cross_product(rightdirector, leftdirector);
		normalize(&normal);
		updirector=cross_product(normal, rightdirector);		//	lab frame
		
		lefttriangleframe.x=dot_product(rightdirector, leftsep);
		lefttriangleframe.y=dot_product(updirector, leftsep);
		lefttriangleframe.z=0;										//	full-magnitude left vector in triangle frame
		righttriangleframe.y=0;
		righttriangleframe.z=0;
        
		//	calculate nearest point in triangle, in the frame of the triangle
		
		postriangleframe.x=dot_product(rightdirector, relativepos);
		postriangleframe.y=dot_product(updirector, relativepos);
		postriangleframe.z=dot_product(normal, relativepos);
		
		sign=1;
		nearestpoint.z=0;
		if(postriangleframe.y<0){
			if(postriangleframe.x<0){
				
				//  still need to check left edge!
				
				lefttriangleframeunit=lefttriangleframe;
				normalize(&lefttriangleframeunit);
				distancealongedge=dot_product(postriangleframe, lefttriangleframeunit);
				if(distancealongedge<0){
					geometrycode=4;         //  center corner
					distance=norm(postriangleframe);
					if(postriangleframe.z<0) sign=-1;
				}
				else if(distancealongedge>leftmag){
					geometrycode=6;         //  left corner
					distance=sqrt(pow(postriangleframe.x-lefttriangleframe.x, 2)+pow(postriangleframe.y-lefttriangleframe.y, 2)+pow(postriangleframe.z, 2));
					if(postriangleframe.z<0) sign=-1;
				}
				else{
					geometrycode=3;         //  center-left edge
					nearestpoint=scalar_multiply_double_triple(lefttriangleframeunit, distancealongedge);        //   do use it here
					postriangleframe.x-=nearestpoint.x;
					postriangleframe.y-=nearestpoint.y;			//	changing postriangleframe to reference frame of nearest point
					distance=norm(postriangleframe);
					normaldotproductmag[i]=fabs(postriangleframe.z);
					if(postriangleframe.z<0) sign=-1;
				}
			}
			else if(postriangleframe.x>righttriangleframe.x){
				
				//  still need to check right-left edge!
				
				rightleftsep=subtract_double_triple(lefttriangleframe, righttriangleframe);
				rightleftmag=norm(rightleftsep);
				rightleftdirector=scalar_multiply_double_triple(rightleftsep, 1./rightleftmag);
				distancealongedge=dot_product(subtract_double_triple(postriangleframe, righttriangleframe), rightleftdirector);
				if(distancealongedge<0){
					geometrycode=5;         //  right corner
					distance=sqrt(pow(postriangleframe.x-righttriangleframe.x, 2)+pow(postriangleframe.y, 2)+pow(postriangleframe.z, 2));
					if(postriangleframe.z<0) sign=-1;
				}
				else if(distancealongedge>rightleftmag){
					geometrycode=6;         //  left corner
					distance=sqrt(pow(postriangleframe.x-lefttriangleframe.x, 2)+pow(postriangleframe.y-lefttriangleframe.y, 2)+pow(postriangleframe.z, 2));
					if(postriangleframe.z<0) sign=-1;
				}
				else{
					geometrycode=2;         //  right-left edge
					nearestpoint=add_double_triple(righttriangleframe, scalar_multiply_double_triple(rightleftdirector, distancealongedge));
					postriangleframe.x-=nearestpoint.x;
					postriangleframe.y-=nearestpoint.y;			//	changing postriangleframe to reference frame of nearest point
				
					distance=norm(postriangleframe);
					normaldotproductmag[i]=fabs(postriangleframe.z);
					if(postriangleframe.z<0) sign=-1;
				}
			}
			else{
				geometrycode=1;         //  center-right edge
				distance=sqrt(pow(postriangleframe.y, 2)+pow(postriangleframe.z, 2));
				normaldotproductmag[i]=fabs(postriangleframe.z);
				if(postriangleframe.z<0) sign=-1;
			}
		}
		else if(lefttriangleframe.y<0) my_exit("lefttriangleframe.y<0");
		else if(postriangleframe.x<postriangleframe.y*lefttriangleframe.x/lefttriangleframe.y){		// to the left of center-left edge
			lefttriangleframeunit=lefttriangleframe;
			normalize(&lefttriangleframeunit);
			distancealongedge=dot_product(postriangleframe, lefttriangleframeunit);
			if(distancealongedge<0){
				geometrycode=4;         //  center corner
				distance=norm(postriangleframe);
				if(postriangleframe.z<0) sign=-1;
			}
			else if(distancealongedge>leftmag){
				
				//  still need to check right-left edge
				
				rightleftsep=subtract_double_triple(lefttriangleframe, righttriangleframe);
				rightleftmag=norm(rightleftsep);
				rightleftdirector=scalar_multiply_double_triple(rightleftsep, 1./rightleftmag);
				distancealongedge=dot_product(subtract_double_triple(postriangleframe, righttriangleframe), rightleftdirector);
				if(distancealongedge<0){
					geometrycode=5;         //  right corner
					distance=sqrt(pow(postriangleframe.x-righttriangleframe.x, 2)+pow(postriangleframe.y, 2)+pow(postriangleframe.z, 2));
					if(postriangleframe.z<0) sign=-1;
				}
				else if(distancealongedge>rightleftmag){
					geometrycode=6;         //  left corner
					distance=sqrt(pow(postriangleframe.x-lefttriangleframe.x, 2)+pow(postriangleframe.y-lefttriangleframe.y, 2)+pow(postriangleframe.z, 2));
					if(postriangleframe.z<0) sign=-1;
				}
				else{
					geometrycode=2;         //  right-left edge
					nearestpoint=add_double_triple(righttriangleframe, scalar_multiply_double_triple(rightleftdirector, distancealongedge));
					postriangleframe.x-=nearestpoint.x;
					postriangleframe.y-=nearestpoint.y;			//	changing postriangleframe to reference frame of nearest point
	
					distance=norm(postriangleframe);
					normaldotproductmag[i]=fabs(postriangleframe.z);
					if(postriangleframe.z<0) sign=-1;
				}
			}
			else{
				geometrycode=3;         //  center-left edge
				nearestpoint=scalar_multiply_double_triple(lefttriangleframeunit, distancealongedge);        //   do use it here
				postriangleframe.x-=nearestpoint.x;
				postriangleframe.y-=nearestpoint.y;			//	changing postriangleframe to reference frame of nearest point
				distance=norm(postriangleframe);
				normaldotproductmag[i]=fabs(postriangleframe.z);
				if(postriangleframe.z<0) sign=-1;
			}
		}
		else if(postriangleframe.x>righttriangleframe.x+postriangleframe.y*(lefttriangleframe.x-righttriangleframe.x)/lefttriangleframe.y){		//	to the outside of left-right edge
			rightleftsep=subtract_double_triple(lefttriangleframe, righttriangleframe);
			rightleftmag=norm(rightleftsep);
			rightleftdirector=scalar_multiply_double_triple(rightleftsep, 1./rightleftmag);
			distancealongedge=dot_product(subtract_double_triple(postriangleframe, righttriangleframe), rightleftdirector);
			if(distancealongedge<0){
				geometrycode=5;         //  right corner
				distance=sqrt(pow(postriangleframe.x-righttriangleframe.x, 2)+pow(postriangleframe.y, 2)+pow(postriangleframe.z, 2));
				if(postriangleframe.z<0) sign=-1;
			}
			else if(distancealongedge>rightleftmag){
				geometrycode=6;         //  left corner
				distance=sqrt(pow(postriangleframe.x-lefttriangleframe.x, 2)+pow(postriangleframe.y-lefttriangleframe.y, 2)+pow(postriangleframe.z, 2));
				if(postriangleframe.z<0) sign=-1;
			}
			else{
				geometrycode=2;         //  right-left edge
				nearestpoint=add_double_triple(righttriangleframe, scalar_multiply_double_triple(rightleftdirector, distancealongedge));
				postriangleframe.x-=nearestpoint.x;
				postriangleframe.y-=nearestpoint.y;			//	changing postriangleframe to reference frame of nearest point
			
				distance=norm(postriangleframe);
				normaldotproductmag[i]=fabs(postriangleframe.z);
				if(postriangleframe.z<0) sign=-1;
			}
		}
		else{
			geometrycode=0;         //  projects onto triangle
			distance=postriangleframe.z;
			if(postriangleframe.z<0){			
				sign=-1;
				distance*=-1;
			}
		}
		distance*=sign;
		if(i==0){
			height=distance;
			tieindex[0]=i;
			tiesign[0]=sign;
			tiecode[0]=geometrycode;
		}
		else if(fabs(distance)<fabs(height)-heightequalitytolerance){
			height=distance;
			tieindex[0]=i;
			tiesign[0]=sign;
			tiecode[0]=geometrycode;
			tienumber=1;				
		}
		else if(fabs(distance)<=fabs(height)+heightequalitytolerance){
			tieindex[tienumber]=i;
			tiesign[tienumber]=sign;
			tiecode[tienumber]=geometrycode;
			tienumber++;
		}
	}
	declare_array(int, edgeties, 4);
	declare_array(int, nearoutsidecorner, 4);
	declare_array(int, reorder, 4);
	int *farcorner, *nextcorner, *conflictingtriangles;
    double projection, projection0, projection1;
	double_triple averagenormal, sep1, normal0, normal1;
	averagenormal.x=averagenormal.y=averagenormal.z=0;
	int exitflag=0, foundsign, firstsign, secondsign, indexdifference, sharededge, nextcornerside;
	double_double *vertices;
	if(tienumber>1){
		for(i=1;i<tienumber;i++){
			if(tiesign[i]!=tiesign[0]){
				allsamesign=0;
				break;
			}
		}
		if(allsamesign==0){
			found=0;
			allsamesign=1;
			for(i=0;i<tienumber;i++){
				if(tiecode[i]==0){
					if(found==2){
						printf("found three direct projections\n");
						height_verbose_debugging(minindex, center, meshpositions, mapmeshtoneighbors[minindex], box_dimension, mycoord.r, relativepos);
						exit(1);
					}
					edgeties[found]=i;								//	using edgeties for direct projections; will reset with found=0
					if(tiesign[i]!=tiesign[edgeties[0]]) allsamesign=0;
					found++;
				}
			}
			if(found==1){
				height=height/tiesign[0]*tiesign[i];                //  if one and only one direct projection, use that sign
			}
			else if(found==2){
				if(allsamesign==1){
					height=height/tiesign[0]*tiesign[edgeties[0]];      //  if two direct projections with the same sign, use it
				}
				else{
					allocate_array(int, farcorner, 2);
					allocate_array(int, nextcorner, 2);
					indexdifference=mod(edgeties[1]-edgeties[0], 6);
					if((indexdifference!=1)&&(indexdifference!=5)){
						printf("direct projection onto two non-adjacent triangles, %i and %i\n", edgeties[0], edgeties[1]);
						height_verbose_debugging(minindex, center, meshpositions, mapmeshtoneighbors[minindex], box_dimension, mycoord.r, relativepos);
						exit(1);
					}
					if(indexdifference==1){
						sharededge=edgeties[1];
						farcorner[1]=mod(edgeties[1]+1, 6);
						farcorner[0]=edgeties[0];
						nextcorner[1]=mod(edgeties[1]+2, 6);
						nextcorner[0]=mod(edgeties[0]-1, 6);
					}
					else{
						sharededge=edgeties[0];
						farcorner[0]=mod(edgeties[0]+1, 6);
						farcorner[1]=edgeties[1];
						nextcorner[0]=mod(edgeties[0]+2, 6);
						nextcorner[1]=mod(edgeties[1]-1, 6);
					}
					
					//	recalculate orthonormal basis for edgeties[0] (arbitrary choice)
					
					right=meshpositions[mapmeshtoneighbors[minindex][edgeties[0]]];
					rightsep=subtract_double_triple(right, center);
					recenter_double_triple(&rightsep, box_dimension);
					rightdirector=rightsep;
					righttriangleframe.x=norm(rightdirector);
					rightdirector=scalar_multiply_double_triple(rightdirector, 1./righttriangleframe.x);		//	lab frame
					left=meshpositions[mapmeshtoneighbors[minindex][mod(edgeties[0]+1, 6)]];
					leftsep=subtract_double_triple(left, center);
					recenter_double_triple(&leftsep, box_dimension);
					leftmag=norm(leftsep);
					leftdirector=scalar_multiply_double_triple(leftsep, 1./leftmag);		//	lab frame
					normal=cross_product(rightdirector, leftdirector);
					normalize(&normal);
					updirector=cross_product(normal, rightdirector);		//	lab frame
					
					allocate_array(double_double, vertices, 4);		//	center, farcorner[0], shared corner, farcorner[1]; no handedness implied
					vertices[0].x=0;
					vertices[0].y=0;
					sep=subtract_double_triple(meshpositions[mapmeshtoneighbors[minindex][farcorner[0]]], center);
					recenter_double_triple(&sep, box_dimension);
					vertices[1].x=dot_product(sep, rightdirector);
					vertices[1].y=dot_product(sep, updirector);
					sep=subtract_double_triple(meshpositions[mapmeshtoneighbors[minindex][sharededge]], center);
					recenter_double_triple(&sep, box_dimension);
					vertices[2].x=dot_product(sep, rightdirector);
					vertices[2].y=dot_product(sep, updirector);
					sep=subtract_double_triple(meshpositions[mapmeshtoneighbors[minindex][farcorner[1]]], center);
					recenter_double_triple(&sep, box_dimension);
					vertices[3].x=dot_product(sep, rightdirector);
					vertices[3].y=dot_product(sep, updirector);
					
					if(checkpointintrianglelist(vertices[1], vertices[0], vertices[2], vertices[3])==1){		//	farcorner 0 inside triangle 1
						sep=subtract_double_triple(meshpositions[mapmeshtoneighbors[minindex][nextcorner[0]]], center);		//	which side is NEXT corner on?
                        recenter_double_triple(&sep, box_dimension);
						if(dot_product(sep, normal)>0) nextcornerside=1;
						else nextcornerside=-1;
						if(tiesign[edgeties[0]]==nextcornerside){
							height=height/tiesign[0]*tiesign[edgeties[0]];										//	farcorner on same side: keep sign
						}
						else{
							height=-height/tiesign[0]*tiesign[edgeties[0]];										//	farcorner on different side: reverse sign
						}
					}
					else if(checkpointintrianglelist(vertices[3], vertices[0], vertices[1], vertices[2])==1){		//	farcorner 1 inside triangle 0
						sep=subtract_double_triple(meshpositions[mapmeshtoneighbors[minindex][nextcorner[1]]], center);		//	which side is NEXT corner on?
                        recenter_double_triple(&sep, box_dimension);
						if(dot_product(sep, normal)>0) nextcornerside=1;
						else nextcornerside=-1;
						if(tiesign[edgeties[1]]==nextcornerside){
							height=height/tiesign[0]*tiesign[edgeties[1]];										//	farcorner on same side: keep sign
						}
						else{
							height=-height/tiesign[0]*tiesign[edgeties[1]];										//	farcorner on different side: reverse sign
						}
					}
					else{
                        sep=subtract_double_triple(meshpositions[mapmeshtoneighbors[minindex][nextcorner[0]]], center);
                        recenter_double_triple(&sep, box_dimension);
                        sep1=subtract_double_triple(meshpositions[mapmeshtoneighbors[minindex][nextcorner[1]]], center);
                        recenter_double_triple(&sep1, box_dimension);
                        projection0=dot_product(sep, normal);
                        projection1=dot_product(sep1, normal);
                        if(projection*projection1<0){               //  different sides; no discrepancy
                            if(projection0>0) nextcornerside=1;
                            else nextcornerside=-1;
                            if(tiesign[edgeties[0]]==nextcornerside){
                                height=height/tiesign[0]*tiesign[edgeties[0]];										//	farcorner on same side: keep sign
                            }
                            else{
                                height=-height/tiesign[0]*tiesign[edgeties[0]];										//	farcorner on different side: reverse sign
                            }
                        }
                        else{
                            
                            //  use next corner whose triangle is more planar with original triangle; either this is the closest or two next triangles must overlap (not allowed)
                            
                            sep=subtract_double_triple(meshpositions[mapmeshtoneighbors[minindex][farcorner[0]]], center);
                            recenter_double_triple(&sep, box_dimension);
                            sep1=subtract_double_triple(meshpositions[mapmeshtoneighbors[minindex][nextcorner[0]]], center);
                            recenter_double_triple(&sep, box_dimension);
                            normal0=cross_product(sep, sep1);
                            normalize(&normal0);
                            sep=subtract_double_triple(meshpositions[mapmeshtoneighbors[minindex][farcorner[1]]], center);
                            recenter_double_triple(&sep, box_dimension);
                            sep1=subtract_double_triple(meshpositions[mapmeshtoneighbors[minindex][nextcorner[1]]], center);
                            recenter_double_triple(&sep, box_dimension);
                            normal1=cross_product(sep, sep1);
                            normalize(&normal1);
                            if(fabs(dot_product(normal0, normal))>fabs(dot_product(normal1, normal))){
                                height=height/tiesign[0]*tiesign[edgeties[0]];										//	farcorner on same side: keep sign
                            }
                            else{
                                height=-height/tiesign[0]*tiesign[edgeties[0]];										//	farcorner on same side: keep sign
                            }
                        }
					}
					free(farcorner);
					free(nextcorner);
				}
			}
			else if(found>2){
				printf("found %i direct projections!\n", found);
				exit(1);
			}
			else{
				found=0;
				allsamesign=1;
				for(i=0;i<tienumber;i++){
					if(tiecode[i]<4){                                   //  looking for edges
						if(found==4){
							printf("site equally close to >4 edges, verbose debugging...\n");
							height_verbose_debugging(minindex, center, meshpositions, mapmeshtoneighbors[minindex], box_dimension, mycoord.r, relativepos);
							my_exit("site equally close to three edges");
						}
						edgeties[found]=i;
						if(tiesign[i]!=tiesign[edgeties[0]]) allsamesign=0;
						found++;
					}
				}
				if((found>0)&&(allsamesign==1)){
					height=height/tiesign[0]*tiesign[edgeties[0]];      //  if at least one edge site and all the same sign, use it
				}
				else{
					if(found==1){
						if(tiecode[edgeties[0]]!=2){
							printf("found only one close interior edge; should be an equivalent one; verbose debugging\n");
							my_exit("found only one close interior edge; should be an equivalent one");
						}
						height=height/tiesign[0]*tiesign[edgeties[0]];      //  if one and only edge projection, use that sign
					}
					else if(found==3){
						for(i=0;i<2;i++){
							if(tiecode[edgeties[i]]==2){
								my_exit("two equally close edge projections, and one is outside edge");     //  if one outside edge, exit, since this means they are different edges
							}
							else if(tiecode[edgeties[i]]==1){                   //  right edge
								nearoutsidecorner[i]=tieindex[edgeties[i]];
							}
							else if(tiecode[edgeties[i]]==3){                   //  left edge
								if(tieindex[edgeties[i]]==5) nearoutsidecorner[i]=0;
								else nearoutsidecorner[i]=tieindex[edgeties[i]]+1;
							}
							else my_exit("far corner not defined");
						}
						found=0;
						if(nearoutsidecorner[1]==nearoutsidecorner[0]){
							reorder[0]=0;
							reorder[1]=1;
							found=2;
						}
						else if(nearoutsidecorner[2]==nearoutsidecorner[0]){
							reorder[0]=0;
							reorder[1]=2;
							found=2;
						}
						else if(nearoutsidecorner[2]==nearoutsidecorner[1]){
							reorder[0]=1;
							reorder[1]=2;
							found=2;
						}
						if(found==2){
							if(tiesign[edgeties[reorder[0]]]==tiesign[edgeties[reorder[1]]]){
								height=height/tiesign[0]*tiesign[edgeties[reorder[0]]];      //  two equally close edge projections, but same sign
							}							
							else if(normaldotproductmag[tieindex[edgeties[reorder[0]]]]>normaldotproductmag[tieindex[edgeties[reorder[1]]]]){
								height=height/tiesign[0]*tiesign[edgeties[reorder[0]]];
							}
							else height=height/tiesign[0]*tiesign[edgeties[reorder[1]]];
						}
						else{
							if(normaldotproductmag[tieindex[edgeties[0]]]>normaldotproductmag[tieindex[edgeties[1]]]){
								if(normaldotproductmag[tieindex[edgeties[0]]]>normaldotproductmag[tieindex[edgeties[2]]]){
									height=height/tiesign[0]*tiesign[edgeties[0]];
								}
								else height=height/tiesign[0]*tiesign[edgeties[2]];
							}
							else{
								if(normaldotproductmag[tieindex[edgeties[1]]]>normaldotproductmag[tieindex[edgeties[2]]]){
									height=height/tiesign[0]*tiesign[edgeties[1]];
								}
								else height=height/tiesign[0]*tiesign[edgeties[2]];
							}
						}						
					}
					else if(found==2){                                          //  one pair of edges
						if(tiesign[edgeties[0]]==tiesign[edgeties[1]]){
							height=height/tiesign[0]*tiesign[edgeties[0]];      //  two equally close edge projections, but same sign
						}
						else{
							for(i=0;i<2;i++){
								if(tiecode[edgeties[i]]==2){
									my_exit("two equally close edge projections, and one is outside edge");     //  if one outside edge, exit, since this means they are different edges
								}
								else if(tiecode[edgeties[i]]==1){                   //  right edge
									nearoutsidecorner[i]=tieindex[edgeties[i]];
								}
								else if(tiecode[edgeties[i]]==3){                   //  left edge
									if(tieindex[edgeties[i]]==5) nearoutsidecorner[i]=0;
									else nearoutsidecorner[i]=tieindex[edgeties[i]]+1;
								}
								else my_exit("far corner not defined");
							}
							if(nearoutsidecorner[0]!=nearoutsidecorner[1]){
								printf("two equally close but different edges\n");
							}
							if(normaldotproductmag[tieindex[edgeties[0]]]>normaldotproductmag[tieindex[edgeties[1]]]){
								height=height/tiesign[0]*tiesign[edgeties[0]];
							}
							else height=height/tiesign[0]*tiesign[edgeties[1]];
						}
					}
					else if(found==4){                                      //  two pairs of edges
						for(i=0;i<4;i++){
							if(tiecode[edgeties[i]]==2){
								my_exit("four equally close edge projections, and one is outside edge");     //  if one outside edge, exit, since this means they are different edges
							}
							else if(tiecode[edgeties[i]]==1){                   //  right edge
								nearoutsidecorner[i]=tieindex[edgeties[i]];
							}
							else if(tiecode[edgeties[i]]==3){                   //  left edge
								if(tieindex[edgeties[i]]==5) nearoutsidecorner[i]=0;
								else nearoutsidecorner[i]=tieindex[edgeties[i]]+1;
							}
							else my_exit("far corner not defined");
						}
						reorder[0]=0;
						if(nearoutsidecorner[1]==nearoutsidecorner[0]){
							reorder[1]=1;
							reorder[2]=2;
							if(nearoutsidecorner[3]==nearoutsidecorner[2]) reorder[3]=3;
							else{
								cornersdontmatchdebug(edgeties, tieindex, nearoutsidecorner, tiecode);
								height_verbose_debugging(minindex, center, meshpositions, mapmeshtoneighbors[minindex], box_dimension, mycoord.r, relativepos);
								my_exit("four edge corners don't match up");
							}
						}
						else if(nearoutsidecorner[2]==nearoutsidecorner[0]){
							reorder[1]=2;
							reorder[2]=1;
							if(nearoutsidecorner[3]==nearoutsidecorner[1]) reorder[3]=3;
							else{
								cornersdontmatchdebug(edgeties, tieindex, nearoutsidecorner, tiecode);
								height_verbose_debugging(minindex, center, meshpositions, mapmeshtoneighbors[minindex], box_dimension, mycoord.r, relativepos);
								my_exit("four edge corners don't match up");
							}
						}
						else if(nearoutsidecorner[3]==nearoutsidecorner[0]){
							reorder[1]=3;
							reorder[2]=1;
							if(nearoutsidecorner[2]==nearoutsidecorner[1]) reorder[3]=2;
							else{
								cornersdontmatchdebug(edgeties, tieindex, nearoutsidecorner, tiecode);
								height_verbose_debugging(minindex, center, meshpositions, mapmeshtoneighbors[minindex], box_dimension, mycoord.r, relativepos);
								my_exit("four edge corners don't match up");
							}
						}
						else{
							cornersdontmatchdebug(edgeties, tieindex, nearoutsidecorner, tiecode);
							height_verbose_debugging(minindex, center, meshpositions, mapmeshtoneighbors[minindex], box_dimension, mycoord.r, relativepos);
							my_exit("four edge corners don't match up");
						}
						
						if(normaldotproductmag[tieindex[edgeties[reorder[0]]]]>normaldotproductmag[tieindex[edgeties[reorder[1]]]]){
							firstsign=tiesign[edgeties[reorder[0]]];
						}
						else firstsign=tiesign[edgeties[reorder[1]]];
						if(normaldotproductmag[tieindex[edgeties[reorder[2]]]]>normaldotproductmag[tieindex[edgeties[reorder[3]]]]){
							secondsign=tiesign[edgeties[reorder[2]]];
						}
						else secondsign=tiesign[edgeties[reorder[3]]];
						
						if(firstsign==secondsign){
							height=height/tiesign[0]*firstsign;
						}
						else{
							
							//	two pairs of edges with different signs
							//	if they are separated by one, assume that two intervening triangles are folded
							//	break tie as done for two direct projections
                            
							allocate_array(int, farcorner, 2);
							allocate_array(int, nextcorner, 2);
							allocate_array(int, conflictingtriangles, 2);		//	here, use conflictingtriangles NOT edgeties to index the two triangles
							if(tiecode[edgeties[reorder[0]]]==1){				//	right
								conflictingtriangles[0]=tieindex[edgeties[reorder[0]]];
							}
							else if(tiecode[edgeties[reorder[0]]]==3){
								conflictingtriangles[0]=mod(tieindex[edgeties[reorder[0]]]+1, 6);
							}
							else{
								printf("tiecode[edgeties[reorder[%i]]] = tieindex[%i] = %i not 1 or 3\n", 0, edgeties[reorder[0]], tiecode[edgeties[reorder[0]]]);
								exit(1);
							}
							if(tiecode[edgeties[reorder[1]]]==1){				//	right
								if(conflictingtriangles[0]!=tieindex[edgeties[reorder[1]]]){
									printf("tieindex[edgeties[reorder[%i]]] = tieindex[%i] = %i not matching\n", 1, edgeties[reorder[1]], tieindex[edgeties[reorder[1]]]);
									exit(1);
								}
							}
							else if(tiecode[edgeties[reorder[1]]]==3){
                                if(conflictingtriangles[0]!=mod(tieindex[edgeties[reorder[1]]]+1, 6)){
									printf("(tieindex[edgeties[reorder[%i]]] + 1) mod 6 = (tieindex[%i] + 1) mod 6 = %i not matching\n", 1, edgeties[reorder[1]], mod(tieindex[edgeties[reorder[1]]]+1, 6));
									exit(1);
								}
							}
							else{
								printf("tiecode[edgeties[reorder[%i]]] = tieindex[%i] = %i not 1 or 3\n", 1, edgeties[reorder[1]], tiecode[edgeties[reorder[1]]]);
								exit(1);
							}
                            
							if(tiecode[edgeties[reorder[2]]]==1){				//	right
								conflictingtriangles[1]=tieindex[edgeties[reorder[2]]];
							}
							else if(tiecode[edgeties[reorder[2]]]==3){
								conflictingtriangles[1]=mod(tieindex[edgeties[reorder[2]]]+1, 6);
							}
							else{
								printf("tiecode[edgeties[reorder[%i]]] = tieindex[%i] = %i not 1 or 3\n", 0, edgeties[reorder[2]], tiecode[edgeties[reorder[2]]]);
								exit(1);
							}
							if(tiecode[edgeties[reorder[3]]]==1){				//	right
								if(conflictingtriangles[1]!=tieindex[edgeties[reorder[3]]]){
									printf("tieindex[edgeties[reorder[%i]]] = tieindex[%i] = %i not matching\n", 1, edgeties[reorder[3]], tieindex[edgeties[reorder[3]]]);
									exit(1);
								}
							}
							else if(tiecode[edgeties[reorder[3]]]==3){
                                if(conflictingtriangles[1]!=mod(tieindex[edgeties[reorder[3]]]+1, 6)){
									printf("(tieindex[edgeties[reorder[%i]]] + 1) mod 6 = (tieindex[%i] + 1) mod 6 = %i not matching\n", 1, edgeties[reorder[3]], mod(tieindex[edgeties[reorder[3]]]+1, 6));
									exit(1);
								}
							}
							else{
								printf("tiecode[edgeties[reorder[%i]]] = tieindex[%i] = %i not 1 or 3\n", 1, edgeties[reorder[3]], tiecode[edgeties[reorder[3]]]);
								exit(1);
							}
							
							indexdifference=mod(conflictingtriangles[1]-conflictingtriangles[0], 6);
                            if((indexdifference!=2)&&(indexdifference!=4)){
								printf("edge projection onto different edges %i and %i not separated by 2\n", conflictingtriangles[0], conflictingtriangles[1]);
								
								printf("deciding by normaldotproductmag, values");
								for(i=0;i<4;i++) printf(" %f", normaldotproductmag[tieindex[edgeties[reorder[i]]]]);
								printf("\n");
								if(normaldotproductmag[tieindex[edgeties[reorder[0]]]]>normaldotproductmag[tieindex[edgeties[reorder[1]]]]){
									if(normaldotproductmag[tieindex[edgeties[reorder[0]]]]>normaldotproductmag[tieindex[edgeties[reorder[2]]]]){
										if(normaldotproductmag[tieindex[edgeties[reorder[0]]]]>normaldotproductmag[tieindex[edgeties[reorder[3]]]]){
											firstsign=tiesign[edgeties[reorder[0]]];
										}
										else firstsign=tiesign[edgeties[reorder[3]]];
									}
									else{
										if(normaldotproductmag[tieindex[edgeties[reorder[2]]]]>normaldotproductmag[tieindex[edgeties[reorder[3]]]]){
											firstsign=tiesign[edgeties[reorder[2]]];
										}
										else firstsign=tiesign[edgeties[reorder[3]]];
									}
								}
								else{
									if(normaldotproductmag[tieindex[edgeties[reorder[1]]]]>normaldotproductmag[tieindex[edgeties[reorder[2]]]]){
										if(normaldotproductmag[tieindex[edgeties[reorder[1]]]]>normaldotproductmag[tieindex[edgeties[reorder[3]]]]){
											firstsign=tiesign[edgeties[reorder[1]]];
										}
										else firstsign=tiesign[edgeties[reorder[3]]];
									}
									else{
										if(normaldotproductmag[tieindex[edgeties[reorder[2]]]]>normaldotproductmag[tieindex[edgeties[reorder[3]]]]){
											firstsign=tiesign[edgeties[reorder[2]]];
										}
										else firstsign=tiesign[edgeties[reorder[3]]];
									}
								}
								height=height/tiesign[0]*firstsign;
								
								printf("debugging (without exit) so I can check that the choice is correct\n");
								height_verbose_debugging(minindex, center, meshpositions, mapmeshtoneighbors[minindex], box_dimension, mycoord.r, relativepos);
								
							}
							if(indexdifference==2){
								sharededge=mod(conflictingtriangles[0]+1, 6);
								farcorner[1]=conflictingtriangles[1];
								farcorner[0]=conflictingtriangles[0];
								nextcorner[1]=mod(conflictingtriangles[1]+1, 6);
								nextcorner[0]=mod(conflictingtriangles[0]-1, 6);
							}
							else{
								sharededge=mod(conflictingtriangles[0]-1, 6);
								farcorner[0]=conflictingtriangles[0];
								farcorner[1]=conflictingtriangles[1];
								nextcorner[0]=mod(conflictingtriangles[0]+1, 6);
								nextcorner[1]=mod(conflictingtriangles[1]-1, 6);
							}
							
							//	recalculate orthonormal basis for conflictingtriangles[0] (arbitrary choice)
							
							right=meshpositions[mapmeshtoneighbors[minindex][conflictingtriangles[0]]];
							rightsep=subtract_double_triple(right, center);
							recenter_double_triple(&rightsep, box_dimension);
							rightdirector=rightsep;
							righttriangleframe.x=norm(rightdirector);
							rightdirector=scalar_multiply_double_triple(rightdirector, 1./righttriangleframe.x);		//	lab frame
							left=meshpositions[mapmeshtoneighbors[minindex][mod(conflictingtriangles[0]+1, 6)]];
							leftsep=subtract_double_triple(left, center);
							recenter_double_triple(&leftsep, box_dimension);
							leftmag=norm(leftsep);
							leftdirector=scalar_multiply_double_triple(leftsep, 1./leftmag);		//	lab frame
							normal=cross_product(rightdirector, leftdirector);
							normalize(&normal);
							updirector=cross_product(normal, rightdirector);		//	lab frame
							
							allocate_array(double_double, vertices, 4);		//	center, farcorner[0], shared corner, farcorner[1]; no handedness implied
							vertices[0].x=0;
							vertices[0].y=0;
							sep=subtract_double_triple(meshpositions[mapmeshtoneighbors[minindex][farcorner[0]]], center);
							recenter_double_triple(&sep, box_dimension);
							vertices[1].x=dot_product(sep, rightdirector);
							vertices[1].y=dot_product(sep, updirector);
							sep=subtract_double_triple(meshpositions[mapmeshtoneighbors[minindex][sharededge]], center);
							recenter_double_triple(&sep, box_dimension);
							vertices[2].x=dot_product(sep, rightdirector);
							vertices[2].y=dot_product(sep, updirector);
							sep=subtract_double_triple(meshpositions[mapmeshtoneighbors[minindex][farcorner[1]]], center);
							recenter_double_triple(&sep, box_dimension);
							vertices[3].x=dot_product(sep, rightdirector);
							vertices[3].y=dot_product(sep, updirector);
							
							if(checkpointintrianglelist(vertices[1], vertices[0], vertices[2], vertices[3])==1){		//	farcorner 0 inside triangle 1
								sep=subtract_double_triple(meshpositions[mapmeshtoneighbors[minindex][nextcorner[0]]], center);		//	which side is NEXT corner on?
								recenter_double_triple(&sep, box_dimension);
								if(dot_product(sep, normal)>0) nextcornerside=1;
								else nextcornerside=-1;
								if(tiesign[edgeties[reorder[0]]]==nextcornerside){
									height=height/tiesign[0]*tiesign[edgeties[reorder[0]]];										//	farcorner on same side: keep sign
								}
								else{
									height=-height/tiesign[0]*tiesign[edgeties[reorder[0]]];										//	farcorner on different side: reverse sign
								}
							}
							else if(checkpointintrianglelist(vertices[3], vertices[0], vertices[1], vertices[2])==1){		//	farcorner 1 inside triangle 0
								sep=subtract_double_triple(meshpositions[mapmeshtoneighbors[minindex][nextcorner[1]]], center);		//	which side is NEXT corner on?
								recenter_double_triple(&sep, box_dimension);
								if(dot_product(sep, normal)>0) nextcornerside=1;
								else nextcornerside=-1;
								if(tiesign[edgeties[reorder[2]]]==nextcornerside){
									height=height/tiesign[0]*tiesign[edgeties[reorder[2]]];										//	farcorner on same side: keep sign
								}
								else{
									height=-height/tiesign[0]*tiesign[edgeties[reorder[2]]];										//	farcorner on different side: reverse sign
								}
							}
							else{
								sep=subtract_double_triple(meshpositions[mapmeshtoneighbors[minindex][nextcorner[0]]], center);
								recenter_double_triple(&sep, box_dimension);
								sep1=subtract_double_triple(meshpositions[mapmeshtoneighbors[minindex][nextcorner[1]]], center);
								recenter_double_triple(&sep1, box_dimension);
								projection0=dot_product(sep, normal);
								projection1=dot_product(sep1, normal);
								if(projection*projection1<0){               //  different sides; no discrepancy
									if(projection0>0) nextcornerside=1;
									else nextcornerside=-1;
									if(tiesign[edgeties[reorder[0]]]==nextcornerside){
										height=height/tiesign[0]*tiesign[edgeties[reorder[0]]];										//	farcorner on same side: keep sign
									}
									else{
										height=-height/tiesign[0]*tiesign[edgeties[reorder[0]]];										//	farcorner on different side: reverse sign
									}
								}
								else{
									
									//  use next corner whose triangle is more planar with original triangle; either this is the closest or two next triangles must overlap (not allowed)
									
									sep=subtract_double_triple(meshpositions[mapmeshtoneighbors[minindex][farcorner[0]]], center);
									recenter_double_triple(&sep, box_dimension);
									sep1=subtract_double_triple(meshpositions[mapmeshtoneighbors[minindex][nextcorner[0]]], center);
									recenter_double_triple(&sep, box_dimension);
									normal0=cross_product(sep, sep1);
									normalize(&normal0);
									sep=subtract_double_triple(meshpositions[mapmeshtoneighbors[minindex][farcorner[1]]], center);
									recenter_double_triple(&sep, box_dimension);
									sep1=subtract_double_triple(meshpositions[mapmeshtoneighbors[minindex][nextcorner[1]]], center);
									recenter_double_triple(&sep, box_dimension);
									normal1=cross_product(sep, sep1);
									normalize(&normal1);
									if(fabs(dot_product(normal0, normal))>fabs(dot_product(normal1, normal))){
										height=height/tiesign[0]*tiesign[edgeties[reorder[0]]];										//	farcorner on same side: keep sign
									}
									else{
										height=-height/tiesign[0]*tiesign[edgeties[reorder[0]]];										//	farcorner on same side: keep sign
									}
								}
							}
							free(farcorner);
							free(nextcorner);
							free(conflictingtriangles);
                            
                            
						}
					}
					else if(found>4){
						printf("found %i edge projections!\n", found);
						exit(1);
					}
					else{           //  all corner projections
						for(i=0;i<tienumber;i++){
							if(tiecode[i]!=4){
								//my_exit("all corner projections but different corners");
								exitflag=1;
							}
							if(exitflag==1){
								printf("all corner projections but different corners; tienumber=%i; debugging...", tienumber);
								for(i=0;i<tienumber;i++){
									printf("\t%i: %i: code %i\n", i, tieindex[i], tiecode[i]);
								}
								exit(1);
							}
							
							//  recalculate normal
							
							right=meshpositions[mapmeshtoneighbors[minindex][i]];
							rightsep=subtract_double_triple(right, center);
							recenter_double_triple(&rightsep, box_dimension);
							rightdirector=rightsep;
							rightdirector=scalar_multiply_double_triple(rightdirector, 1./norm(rightdirector));
							if(i==5) left=meshpositions[mapmeshtoneighbors[minindex][0]];
							else left=meshpositions[mapmeshtoneighbors[minindex][i+1]];
							leftsep=subtract_double_triple(left, center);
							recenter_double_triple(&leftsep, box_dimension);
							leftmag=norm(leftsep);
							leftdirector=scalar_multiply_double_triple(leftsep, 1./leftmag);
							normal=cross_product(rightdirector, leftdirector);
							normalize(&normal);
							
							averagenormal=add_double_triple(averagenormal, normal);     //  sign based on average normal
							
						}
						if(dot_product(relativepos, averagenormal)<0){
							height=height/tiesign[0]*(-1);
						}
						else height=height/tiesign[0];
					}
				}
			}
		}
	}
	
	free(nearoutsidecorner);
	free(edgeties);
	free(reorder);
	free(tieindex);
	free(tiesign);
	free(tiecode);
	free(normaldotproductmag);
	
	if(minindex<solvationparams.meshsize.x*solvationparams.meshsize.y){
		height*=-1;
	}
	if(height>solvationparams.meshsurfacecutoff) height=solvationparams.meshsurfacecutoff;
	else if(height<-solvationparams.meshsurfacecutoff) height=-solvationparams.meshsurfacecutoff;
    
    if(isnan(height)>0){
        printf("height=%f!!", height);
        exit(1);
    }
	return height;
}

int change_cells_doublylinked(linkedlistfull *plink, double_triple box_dimension){
	int flag=0, i, j;
	int_triple newcellsperside;
	newcellsperside.x=floor(box_dimension.x/(*plink).mincellwidth)+1;
	newcellsperside.y=floor(box_dimension.y/(*plink).mincellwidth)+1;
	newcellsperside.z=floor(box_dimension.z/(*plink).mincellwidth)+1;
	if(newcellsperside.x!=(*plink).core.cellsperside.x) flag=1;
	if(newcellsperside.y!=(*plink).core.cellsperside.y) flag=1;
	if(newcellsperside.z!=(*plink).core.cellsperside.z) flag=1;
	if((newcellsperside.x<3)||(newcellsperside.y<3)||(newcellsperside.y<3)) my_exit("code not written for small cellsperside (change_cells_doublylinked).");
	if(flag==1){
		for(i=0;i<(*plink).core.cellsperside.x;i++){
			for(j=0;j<(*plink).core.cellsperside.y;j++){
				free((*plink).core.head[i][j]);
			}
			free((*plink).core.head[i]);
		}
		free((*plink).core.head);
		(*plink).core.head=xcalloc(newcellsperside.x, sizeof(int **));
		for(i=0;i<newcellsperside.x;i++){
			(*plink).core.head[i]=xcalloc(newcellsperside.y, sizeof(int *));
			for(j=0;j<newcellsperside.y;j++){
				(*plink).core.head[i][j]=xcalloc(newcellsperside.z, sizeof(int));
			}
		}
	}
	(*plink).core.cellsperside=newcellsperside;
	(*plink).cellwidth.x=box_dimension.x/(*plink).core.cellsperside.x;
	(*plink).cellwidth.y=box_dimension.y/(*plink).core.cellsperside.y;
	(*plink).cellwidth.z=box_dimension.z/(*plink).core.cellsperside.z;
	return flag;
}

int change_cells_asymmetric(linkedlistasymmetric *plink, double_triple box_dimension){
	int flag=0, i, j;
	int_triple newcellsperside;
	newcellsperside.x=floor(box_dimension.x/(*plink).mincellwidth)+1;
	newcellsperside.y=floor(box_dimension.y/(*plink).mincellwidth)+1;
	newcellsperside.z=floor(box_dimension.z/(*plink).mincellwidth)+1;
	if(newcellsperside.x!=(*plink).core.cellsperside.x) flag=1;
	if(newcellsperside.y!=(*plink).core.cellsperside.y) flag=1;
	if(newcellsperside.z!=(*plink).core.cellsperside.z) flag=1;
	if((newcellsperside.x<3)||(newcellsperside.y<3)||(newcellsperside.y<3)) my_exit("code not written for small cellsperside (change_cells_asymmetric).");
	if(flag==1){
		for(i=0;i<(*plink).core.cellsperside.x;i++){
			for(j=0;j<(*plink).core.cellsperside.y;j++){
				free((*plink).core.head[i][j]);
				free((*plink).headreverse[i][j]);
			}
			free((*plink).core.head[i]);
			free((*plink).headreverse[i]);
		}
		free((*plink).core.head);
		free((*plink).headreverse);
		(*plink).core.head=xcalloc(newcellsperside.x, sizeof(int **));
		(*plink).headreverse=xcalloc(newcellsperside.x, sizeof(int **));
		for(i=0;i<newcellsperside.x;i++){
			(*plink).core.head[i]=xcalloc(newcellsperside.y, sizeof(int *));
			(*plink).headreverse[i]=xcalloc(newcellsperside.y, sizeof(int *));
			for(j=0;j<newcellsperside.y;j++){
				(*plink).core.head[i][j]=xcalloc(newcellsperside.z, sizeof(int));
				(*plink).headreverse[i][j]=xcalloc(newcellsperside.z, sizeof(int));
			}
		}
	}
	(*plink).core.cellsperside=newcellsperside;
	(*plink).cellwidth.x=box_dimension.x/(*plink).core.cellsperside.x;
	(*plink).cellwidth.y=box_dimension.y/(*plink).core.cellsperside.y;
	(*plink).cellwidth.z=box_dimension.z/(*plink).core.cellsperside.z;
	return flag;
}

void vmmc_translate_wholepolymer(int Nchains, int *chainlength, coord *newcoordarray, coord *reversecoordarray, coord *coordarray, double max_translate, double_triple box_dimension, linkedlistset *plinkset, bonded_params my_bonded_params, nonbonded_params my_nonbonded_params, double temperature, monomernodes **monomerid, int *number_neighbors, int **neighbor, int max_frustrated_links, int max_neighbors, int *interior_member, int *exterior_member, int *in_cluster, int_double *frustrated_link, int *assigned, int move3d, linkedlistasymmetric *pmeshlink, double_triple *meshpositions, int **mapmeshtoneighbors){
	double individual_energy, collective_energy, probability, reverse_individual_energy, reverse_probability;
	int seed=rand_int(Nchains), cluster_complete=0, number_interior_members=0, number_exterior_members=1, cluster_number=1, i, j, k, l, m, candidate, index, firstnode, lastnode, exterior_member_rank=0, number_unchecked_neighbors, neighbor_rank, frustrated_link_index=0, temp, firstnodecandidate, lastnodecandidate, found;
	
	for(i=0;i<Nchains;i++) in_cluster[i]=0;
	in_cluster[seed]=1;
    double_triple shift;
    if(move3d==0){
        shift.x=(rand_int(2)*2-1)*max_translate*rand_double;
        shift.y=(rand_int(2)*2-1)*max_translate*rand_double;
        shift.z=0;
    }
	else{
        shift=scalar_multiply_double_triple(rand_unit_cube(), (rand_int(2)*2-1)*max_translate);
    }
    
	index=seed;
    
	exterior_member[0]=index;
	firstnode=monomerid[index][0].backbone;
	lastnode=monomerid[index][chainlength[index]-1].sidechain;
	if(coordarray[lastnode].nodetype<0) lastnode=monomerid[index][chainlength[index]-1].backbone;						//	in case long terminus
	for(i=firstnode;i<=lastnode;i++){
 		newcoordarray[i]=coordarray[i];
		reversecoordarray[i]=coordarray[i];
		newcoordarray[i].r=add_double_triple(coordarray[i].r, shift);
		reversecoordarray[i].r=subtract_double_triple(coordarray[i].r, shift);											//	don't need to apply b.c. because calc_pair computes minimum image
		if(my_nonbonded_params.solvationparams.interface==1){
            newcoordarray[i].height=height_relative_to_interface(i, newcoordarray[i], my_nonbonded_params.solvationparams, (*pmeshlink), box_dimension, meshpositions, mapmeshtoneighbors);
            reversecoordarray[i].height=height_relative_to_interface(i, reversecoordarray[i], my_nonbonded_params.solvationparams, (*pmeshlink), box_dimension, meshpositions, mapmeshtoneighbors);
        }
	}
	while(cluster_complete==0){
		number_unchecked_neighbors=number_neighbors[index];
		while((number_unchecked_neighbors>0)&&(cluster_complete==0)){
			neighbor_rank=rand_int(number_unchecked_neighbors);						//	choose randomly from neighbors
			candidate=neighbor[index][neighbor_rank];
			for(i=neighbor_rank;i<number_unchecked_neighbors-1;i++){				//	shift neighbor list to close gap
				neighbor[index][i]=neighbor[index][i+1];
			}
			neighbor[index][i]=candidate;
			number_unchecked_neighbors--;
			if(in_cluster[candidate]==0){											//	not already in cluster
				individual_energy=collective_energy=0;
                
				firstnodecandidate=monomerid[candidate][0].backbone;
				lastnodecandidate=monomerid[candidate][chainlength[index]-1].sidechain;
				if(coordarray[lastnodecandidate].nodetype<0) lastnodecandidate=monomerid[candidate][chainlength[index]-1].backbone;

				for(i=firstnode;i<=lastnode;i++){
					individual_energy+=calc_hard_energy_cellstruct_givenpolymer((*plinkset).shortrange.core.number_neighbors[i], (*plinkset).shortrange.core.neighbor[i], newcoordarray[i], coordarray, box_dimension, my_nonbonded_params, candidate);
					collective_energy+=calc_hard_energy_cellstruct_givenpolymer((*plinkset).shortrange.core.number_neighbors[i], (*plinkset).shortrange.core.neighbor[i], coordarray[i], coordarray, box_dimension, my_nonbonded_params, candidate);
                    
					if(newcoordarray[i].nodetype==1){
						individual_energy+=calc_phenyl_energy_cellstruct_givenpolymer((*plinkset).phenylphenyl.core.number_neighbors[i], (*plinkset).phenylphenyl.core.neighbor[i], newcoordarray[i], coordarray, box_dimension, my_nonbonded_params, candidate);
						collective_energy+=calc_phenyl_energy_cellstruct_givenpolymer((*plinkset).phenylphenyl.core.number_neighbors[i], (*plinkset).phenylphenyl.core.neighbor[i], coordarray[i], coordarray, box_dimension, my_nonbonded_params, candidate);
                    }
					else if((newcoordarray[i].nodetype==2)||(newcoordarray[i].nodetype==3)){
						individual_energy+=calc_charged_energy_cellstruct_givenpolymer((*plinkset).chargedcharged.core.number_neighbors[i], (*plinkset).chargedcharged.core.neighbor[i], newcoordarray[i], coordarray, box_dimension, my_nonbonded_params, candidate);
						collective_energy+=calc_charged_energy_cellstruct_givenpolymer((*plinkset).chargedcharged.core.number_neighbors[i], (*plinkset).chargedcharged.core.neighbor[i], coordarray[i], coordarray, box_dimension, my_nonbonded_params, candidate);
					}
				}
                
                if(collective_energy==0){
                    
                    //  remove candidate from index's neighbor list (currently at position number_unchecked_neighbors
                    
                    for(i=number_unchecked_neighbors;i<number_neighbors[index]-1;i++){
                        neighbor[index][i]=neighbor[index][i+1];
                    }
                    number_neighbors[index]--;
                    
                    //  remove index from candidate's neighbor list
                    
                    found=0;
                    for(i=0;i<number_neighbors[candidate];i++){
                        if(neighbor[candidate][i]==index){
                            for(j=i;j<number_neighbors[candidate]-1;j++){
                                neighbor[candidate][j]=neighbor[candidate][j+1];
                            }
                            number_neighbors[candidate]--;
                            found=1;
                            break;
                        }
                    }
                    if(found==0){
                        printf("reverse map not found!!\n");
                        exit(1);
                    }
				}
				else{           
                    probability=1.0-exp(-(individual_energy-collective_energy)/temperature);
                    if(probability>0){
                        if((rand_double<probability)||(probability==1)){						//	pre-link
                            reverse_individual_energy=0;
                            for(i=firstnode;i<=lastnode;i++){
                                reverse_individual_energy+=calc_hard_energy_cellstruct_givenpolymer((*plinkset).shortrange.core.number_neighbors[i], (*plinkset).shortrange.core.neighbor[i], reversecoordarray[i], coordarray, box_dimension, my_nonbonded_params, candidate);
                                if(newcoordarray[i].nodetype==1){
                                    reverse_individual_energy+=calc_phenyl_energy_cellstruct_givenpolymer((*plinkset).phenylphenyl.core.number_neighbors[i], (*plinkset).phenylphenyl.core.neighbor[i], reversecoordarray[i], coordarray, box_dimension, my_nonbonded_params, candidate);
                                }
                                else if((newcoordarray[i].nodetype==2)||(newcoordarray[i].nodetype==3)){
                                    reverse_individual_energy+=calc_charged_energy_cellstruct_givenpolymer((*plinkset).chargedcharged.core.number_neighbors[i], (*plinkset).chargedcharged.core.neighbor[i], reversecoordarray[i], coordarray, box_dimension, my_nonbonded_params, candidate);
                                }
                            }
                            reverse_probability=1.0-exp(-(reverse_individual_energy-collective_energy)/temperature);
                            
                            if((reverse_probability>probability)||(rand_double<reverse_probability/probability)){		//	form link
                                cluster_number++;
                                exterior_member[number_exterior_members]=candidate;
                                number_exterior_members++;
                                in_cluster[candidate]=1;
                                for(i=firstnodecandidate;i<=lastnodecandidate;i++){
                                    newcoordarray[i].r=add_double_triple(coordarray[i].r, shift);
                                    reversecoordarray[i].r=subtract_double_triple(coordarray[i].r, shift);					//	don't need to apply b.c. because calc_pair computes minimum image
                                    if(my_nonbonded_params.solvationparams.interface==1){
                                        newcoordarray[i].height=height_relative_to_interface(i, newcoordarray[i], my_nonbonded_params.solvationparams, (*pmeshlink), box_dimension, meshpositions, mapmeshtoneighbors);
                                        reversecoordarray[i].height=height_relative_to_interface(i, reversecoordarray[i], my_nonbonded_params.solvationparams, (*pmeshlink), box_dimension, meshpositions, mapmeshtoneighbors);
                                    }
                                }
                            }
                            else{
                                frustrated_link[frustrated_link_index].x=index;				//	mark LINK as frustrated
                                frustrated_link[frustrated_link_index].y=candidate;
                                frustrated_link_index++;
                                if(frustrated_link_index==max_frustrated_links){
                                    printf("max_frustrated_links=%i not big enough.\n", max_frustrated_links);
                                    exit(1);
                                }
                            }
                        }
                    }
				}
			}
		}
		interior_member[number_interior_members]=index;								//	searched all neighbors; shift index from exterior to interior
		number_interior_members++;
		for(i=exterior_member_rank;i<number_exterior_members-1;i++){
			exterior_member[i]=exterior_member[i+1];
		}
		number_exterior_members--;
		if(number_exterior_members==0) cluster_complete=1;
		else{																		//	assign new cluster member to check neighbors of
			exterior_member_rank=rand_int(number_exterior_members);
			index=exterior_member[exterior_member_rank];
			firstnode=monomerid[index][0].backbone;
			lastnode=monomerid[index][chainlength[index]-1].sidechain;
			if(coordarray[lastnode].nodetype<0) lastnode=monomerid[index][chainlength[index]-1].backbone;
		}
	}
	for(i=0;i<frustrated_link_index;i++){
		if((in_cluster[frustrated_link[i].x]==0)||(in_cluster[frustrated_link[i].y]==0)){
			return;																	//	if there are any frustrated LINKS external to the cluster, reject
		}
	}
	
	//	Acceptance prob is min(1, x), where x is a product over pairs {i, j}
	//	where i is in cluster and j is outside
	//	energy either
	//		initially overlapping (positive, neighbor) and finally noninteracting (0)
	//		or initially noninteracting (0, not a neighbor but check in case it is a neighbor) and finally overlapping (positive)
	//	Note that if I use D(C) not equal to 1, then that becomes the acceptance probability. But for now I am putting the size dependence in the distribution of max_number.
    
	double sum_energy=0, initial_energy, final_energy;
	declare_array(int, neighborofindex, Nchains);
	declare_array(double, initiallynoninteractingenergiesbypolymer, Nchains);
	for(l=0;l<number_interior_members;l++){
		index=interior_member[l];
		firstnode=monomerid[index][0].backbone;
		lastnode=monomerid[index][chainlength[index]-1].sidechain;
		if(coordarray[lastnode].nodetype<0) lastnode=monomerid[index][chainlength[index]-1].backbone;
		for(m=0;m<Nchains;m++){
			neighborofindex[m]=0;
			initiallynoninteractingenergiesbypolymer[m]=0;
		}
		neighborofindex[index]=1;		//	self
		for(m=0;m<number_neighbors[index];m++){
			candidate=neighbor[index][m];
			neighborofindex[candidate]=1;
			if(in_cluster[candidate]==0){
				firstnodecandidate=monomerid[candidate][0].backbone;
				lastnodecandidate=monomerid[candidate][chainlength[index]-1].sidechain;
				if(coordarray[lastnode].nodetype<0) lastnodecandidate=monomerid[candidate][chainlength[index]-1].backbone;
				initial_energy=0;
                
				for(i=firstnode;i<=lastnode;i++){
					for(j=firstnodecandidate;j<=lastnodecandidate;j++){
						initial_energy+=calc_energy_pair(coordarray[i], coordarray[j], my_nonbonded_params, box_dimension);
					}
				}
				
				if(initial_energy>0){
					final_energy=0;
					for(i=firstnode;i<=lastnode;i++){
						final_energy+=calc_hard_energy_cellstruct_givenpolymer((*plinkset).shortrange.core.number_neighbors[i], (*plinkset).shortrange.core.neighbor[i], newcoordarray[i], coordarray, box_dimension, my_nonbonded_params, candidate);
						if(newcoordarray[i].nodetype==1){
							final_energy+=calc_phenyl_energy_cellstruct_givenpolymer((*plinkset).phenylphenyl.core.number_neighbors[i], (*plinkset).phenylphenyl.core.neighbor[i], newcoordarray[i], coordarray, box_dimension, my_nonbonded_params, candidate);
						}
						else if((newcoordarray[i].nodetype==2)||(newcoordarray[i].nodetype==3)){
							final_energy+=calc_charged_energy_cellstruct_givenpolymer((*plinkset).chargedcharged.core.number_neighbors[i], (*plinkset).chargedcharged.core.neighbor[i], newcoordarray[i], coordarray, box_dimension, my_nonbonded_params, candidate);
						}
					}
					
					if(final_energy==0) sum_energy-=initial_energy;
				}
				else if(initial_energy==0){
					final_energy=0;
					for(i=firstnode;i<=lastnode;i++){
						final_energy+=calc_hard_energy_cellstruct_givenpolymer((*plinkset).shortrange.core.number_neighbors[i], (*plinkset).shortrange.core.neighbor[i], newcoordarray[i], coordarray, box_dimension, my_nonbonded_params, candidate);
						if(newcoordarray[i].nodetype==1){
							final_energy+=calc_phenyl_energy_cellstruct_givenpolymer((*plinkset).phenylphenyl.core.number_neighbors[i], (*plinkset).phenylphenyl.core.neighbor[i], newcoordarray[i], coordarray, box_dimension, my_nonbonded_params, candidate);
						}
						else if((newcoordarray[i].nodetype==2)||(newcoordarray[i].nodetype==3)){
							final_energy+=calc_charged_energy_cellstruct_givenpolymer((*plinkset).chargedcharged.core.number_neighbors[i], (*plinkset).chargedcharged.core.neighbor[i], newcoordarray[i], coordarray, box_dimension, my_nonbonded_params, candidate);
						}
					}
					if(final_energy>0) sum_energy+=final_energy;
				}
			}
		}
		//	next, add up energies for any polymers that were not initially neighbors
		//	bin the energies by polymers, because only positive pairwise polymer-polymer energies contribute to acceptance criterion
		
		for(i=firstnode;i<lastnode;i++){
			fmod_double_triple(&(newcoordarray[i].r), box_dimension);                   //  need to apply boundary conditions because I look up cells, not just neighbors, here
            add_hardenergy_if_polymer_not_in_list(newcoordarray[i], (*plinkset).shortrange, coordarray, neighborofindex, initiallynoninteractingenergiesbypolymer, my_nonbonded_params, box_dimension);
			if(coordarray[i].nodetype==1) add_phenylphenylenergy_if_polymer_not_in_list(newcoordarray[i], (*plinkset).phenylphenyl, coordarray, neighborofindex, initiallynoninteractingenergiesbypolymer, my_nonbonded_params, box_dimension);
			else if((coordarray[i].nodetype==2)||(coordarray[i].nodetype==3)){
                add_chargedchargedenergy_if_polymer_not_in_list(newcoordarray[i], (*plinkset).chargedcharged, coordarray, neighborofindex, initiallynoninteractingenergiesbypolymer, my_nonbonded_params, box_dimension);
            }
		}
		for(m=0;m<Nchains;m++){
			if(initiallynoninteractingenergiesbypolymer[m]>0) sum_energy+=initiallynoninteractingenergiesbypolymer[m];
		}
	}
	free(neighborofindex);
	free(initiallynoninteractingenergiesbypolymer);
	
	//	external potential
	
	for(l=0;l<number_interior_members;l++){
		index=interior_member[l];
		firstnode=monomerid[index][0].backbone;
		lastnode=monomerid[index][chainlength[index]-1].sidechain;
		if(coordarray[lastnode].nodetype<0) lastnode=monomerid[index][chainlength[index]-1].backbone;
		for(i=firstnode;i<=lastnode;i++){
			sum_energy+=calc_onebody_difference(coordarray[i], newcoordarray[i], chooseonebody(newcoordarray[i].nodetype, my_nonbonded_params), my_nonbonded_params.solvationparams);
            
            //  Treat change in phenyl-phenyl interaction due to change in height as external potential
            
            if(my_nonbonded_params.solvationparams.interface==1){
                sum_energy-=calc_phenyl_energy_cellstruct_polymerlist_singlecount((*plinkset).phenylphenyl.core.number_neighbors[i], (*plinkset).phenylphenyl.core.neighbor[i], coordarray[i], coordarray, box_dimension, my_nonbonded_params, in_cluster);
                sum_energy+=calc_phenyl_energy_cellstruct_polymerlist_singlecount((*plinkset).phenylphenyl.core.number_neighbors[i], (*plinkset).phenylphenyl.core.neighbor[i], newcoordarray[i], newcoordarray, box_dimension, my_nonbonded_params, in_cluster);
            }
            
		}
	}
	probability=exp(-sum_energy/temperature);
	if((probability<1.0)&&(rand_double>probability)){			//	reject
		return;
	}
    
	//	otherwise accept
    
	//	update coordinates and cells
    
	for(l=0;l<number_interior_members;l++){
		index=interior_member[l];
		firstnode=monomerid[index][0].backbone;
		lastnode=monomerid[index][chainlength[index]-1].sidechain;
		if(coordarray[lastnode].nodetype<0) lastnode=monomerid[index][chainlength[index]-1].backbone;
		for(i=firstnode;i<=lastnode;i++){
			fmod_double_triple(&(newcoordarray[i].r), box_dimension);
            if(my_nonbonded_params.solvationparams.interface==1) updatecell_asymmetric(newcoordarray[i].r, &(*pmeshlink), i);
			updatecell(newcoordarray[i].r, &((*plinkset).shortrange), i);
			if(newcoordarray[i].nodetype==1) updatecell(newcoordarray[i].r, &((*plinkset).phenylphenyl), i);
			else if((newcoordarray[i].nodetype==2)||(newcoordarray[i].nodetype==3)) updatecell(newcoordarray[i].r, &((*plinkset).chargedcharged), i);
			coordarray[i]=newcoordarray[i];
		}
	}
    
	//	update site-site neighbor list
	
	for(l=0;l<number_interior_members;l++){
		index=interior_member[l];
		firstnode=monomerid[index][0].backbone;
		lastnode=monomerid[index][chainlength[index]-1].sidechain;
		if(coordarray[lastnode].nodetype<0) lastnode=monomerid[index][chainlength[index]-1].backbone;
		for(i=firstnode;i<=lastnode;i++){
            if(my_nonbonded_params.solvationparams.interface==1) updateneighborlist_asymmetric(&(*pmeshlink), i, coordarray[i].r, meshpositions, box_dimension);
			updateneighborlist(&((*plinkset).shortrange), i, coordarray, box_dimension);
			if(newcoordarray[i].nodetype==1) updateneighborlist_pairenergies_phenyl(&((*plinkset).phenylphenyl), i, coordarray, box_dimension, my_nonbonded_params);
			else if((newcoordarray[i].nodetype==2)||(newcoordarray[i].nodetype==3)) updateneighborlist_pairenergies_charged(&((*plinkset).chargedcharged), i, coordarray, box_dimension, my_nonbonded_params);
		}
	}
    
	//	Update polymer-polymer neighbor lists. First, wipe out neighbor lists starting and ending in the cluster
	
	for(i=0;i<number_interior_members;i++){
		index=interior_member[i];
		for(j=0;j<number_neighbors[index];j++){
			temp=neighbor[index][j];
			if(in_cluster[temp]==0){
				for(k=0;k<number_neighbors[temp];k++){
					if(neighbor[temp][k]==index) break;								//	find the reverse neighbor map; can I do this more efficiently?
				}
				if(k==number_neighbors[temp]){
					printf("Error: did not find reverse map.\n");
					exit(1);
				}
				number_neighbors[temp]--;											//	remove reverse neighbor map
				for(l=k;l<number_neighbors[temp];l++){
					neighbor[temp][l]=neighbor[temp][l+1];
				}
			}
		}
		number_neighbors[index]=0;
	}
    
	//	Now update neighbor lists
    
	int target, targetpolymer;
	double rsq;
	double_triple r;
	for(l=0;l<number_interior_members;l++){
		index=interior_member[l];
		firstnode=monomerid[index][0].backbone;
		lastnode=monomerid[index][chainlength[index]-1].sidechain;
		if(coordarray[lastnode].nodetype<0) lastnode=monomerid[index][chainlength[index]-1].backbone;
		for(i=0;i<Nchains;i++) assigned[i]=0;
		assigned[index]=1;
        
        for(i=0;i<number_neighbors[index];i++) assigned[neighbor[index][i]]=1;
		
        for(i=firstnode;i<=lastnode;i++){
			
			//	neighbor if interacting (but don't need to check hard core interactions, since these infinite interactions would not have been allowed)
			
			if(coordarray[i].nodetype==1){
                for(j=0;j<(*plinkset).phenylphenyl.core.number_neighbors[i];j++){
					target=(*plinkset).phenylphenyl.core.neighbor[i][j];
					targetpolymer=coordarray[target].chainid;
					if(assigned[targetpolymer]==0){
						r=subtract_double_triple(coordarray[target].r, coordarray[i].r);
						rsq=dot_product(r, r);
						if(rsq<my_nonbonded_params.phenylcutoff2){				//	interacting
							assigned[targetpolymer]=1;
							neighbor[index][number_neighbors[index]]=targetpolymer;
							number_neighbors[index]++;
							if(number_neighbors[index]==max_neighbors) my_exit("max_neighbors not big enough (a)");
							if(in_cluster[targetpolymer]==0){					//	safe to add reverse neighbor link without double counting
								neighbor[targetpolymer][number_neighbors[targetpolymer]]=index;
								number_neighbors[targetpolymer]++;
								if(number_neighbors[targetpolymer]==max_neighbors) my_exit("max_neighbors not big enough (b)");
							}
							else{
								found=0;
								for(temp=0;temp<number_neighbors[targetpolymer];temp++){			//	only add link if target->index doesn't already exist
									if(neighbor[targetpolymer][temp]==index) found=1;
								}
								if(found==0){
									neighbor[targetpolymer][number_neighbors[targetpolymer]]=index;
									number_neighbors[targetpolymer]++;
									if(number_neighbors[targetpolymer]==max_neighbors) my_exit("max_neighbors not big enough (c)");
								}
							}
						}
					}
				}
			}
			else if((coordarray[i].nodetype==2)||(coordarray[i].nodetype==3)){
                for(j=0;j<(*plinkset).chargedcharged.core.number_neighbors[i];j++){
					target=(*plinkset).chargedcharged.core.neighbor[i][j];
					targetpolymer=coordarray[target].chainid;
					if(assigned[targetpolymer]==0){
						r=subtract_double_triple(coordarray[target].r, coordarray[i].r);
						rsq=dot_product(r, r);
						if(rsq<my_nonbonded_params.cutoff2){				//	interacting
							assigned[targetpolymer]=1;
							neighbor[index][number_neighbors[index]]=targetpolymer;
							number_neighbors[index]++;
							if(number_neighbors[index]==max_neighbors) my_exit("max_neighbors not big enough (d)");
							if(in_cluster[targetpolymer]==0){					//	safe to add reverse neighbor link without double counting
								neighbor[targetpolymer][number_neighbors[targetpolymer]]=index;
								number_neighbors[targetpolymer]++;
 								if(number_neighbors[targetpolymer]==max_neighbors) my_exit("max_neighbors not big enough (e)");
							}
							else{
								found=0;
								for(temp=0;temp<number_neighbors[targetpolymer];temp++){			//	only add link if target->index doesn't already exist
									if(neighbor[targetpolymer][temp]==index) found=1;
								}
								if(found==0){
									neighbor[targetpolymer][number_neighbors[targetpolymer]]=index;
									number_neighbors[targetpolymer]++;
 									if(number_neighbors[targetpolymer]==max_neighbors) my_exit("max_neighbors not big enough (f)");
								}
							}
						}
					}
				}
			}
		}
	}    
}		

double_triple com(int firstnode, int lastnode, coord *coordarray, double_triple box_dimension){
	int realnodes=0, i;
	double_triple result, sep, pos;
	for(i=firstnode;i<=lastnode;i++){
        if(coordarray[i].nodetype>=0){
            realnodes++;
            if(i==firstnode){
                pos=coordarray[i].r;
                result=pos;
            }
            else{
                sep=subtract_double_triple(coordarray[i].r, coordarray[i-1].r);
                recenter_double_triple(&sep, box_dimension);
                pos=add_double_triple(pos, sep);
                result=add_double_triple(result, pos);
            }
        }
	}
	result=scalar_multiply_double_triple(result, (1./realnodes));
    fmod_double_triple(&result, box_dimension);
	return result;
}

int calccomcheckr(double_triple *pcom, int firstnode, int lastnode, coord *coordarray, double_triple box_dimension, double *pmaxr){
	int realnodes=0, i;
	double_triple sep, pos;
	for(i=firstnode;i<=lastnode;i++){
        if(coordarray[i].nodetype>=0){
            realnodes++;
            if(i==firstnode){
                pos=coordarray[i].r;
                (*pcom)=pos;
            }
            else{
                sep=subtract_double_triple(coordarray[i].r, coordarray[i-1].r);
                recenter_double_triple(&sep, box_dimension);
                pos=add_double_triple(pos, sep);
                (*pcom)=add_double_triple((*pcom), pos);
            }
        }
	}
	(*pcom)=scalar_multiply_double_triple((*pcom), (1./realnodes));
    fmod_double_triple(&(*pcom), box_dimension);
	for(i=firstnode;i<=lastnode;i++){
        if(coordarray[i].nodetype>=0){
			sep=subtract_double_triple(coordarray[i].r, (*pcom));
			recenter_double_triple(&sep, box_dimension);
			if(norm(sep)>(*pmaxr)){
				//printf("site distance from com=%f!\n", norm(sep));
                printf("increasing maxrforvmmc from %f to %f\n", (*pmaxr), norm(sep)+maxrforvmmcbuffer);
                (*pmaxr)=norm(sep)+maxrforvmmcbuffer;
				return 1;
			}
		}
	}
	return 0;
}

int_triple cellofpos(double_triple pos, double_triple cellwidth){
	int_triple new;
	new.x=(int) floor(pos.x/cellwidth.x);
	new.y=(int) floor(pos.y/cellwidth.y);
	new.z=(int) floor(pos.z/cellwidth.z);
	return new;
}

void vmmc_rotate_wholepolymer(int Nchains, int *chainlength, coord *newcoordarray, coord *reversecoordarray, coord *coordarray, double max_rotate, double_triple box_dimension, linkedlistset *plinkset, bonded_params my_bonded_params, nonbonded_params my_nonbonded_params, double temperature, monomernodes **monomerid, int *number_neighbors, int **neighbor, int max_frustrated_links, int max_neighbors, int *interior_member, int *exterior_member, int *in_cluster, int_double *frustrated_link, int *assigned, int move3d, double *pmaxr, linkedlistasymmetric *pmeshlink, double_triple *meshpositions, int **mapmeshtoneighbors){
	double individual_energy, collective_energy, probability, reverse_individual_energy, reverse_probability;
	int checkrflag, seed=rand_int(Nchains), cluster_complete=0, number_interior_members=0, number_exterior_members=1, cluster_number=1, i, j, k, l, m, candidate, index, firstnode, lastnode, exterior_member_rank=0, number_unchecked_neighbors, neighbor_rank, frustrated_link_index=0, temp, firstnodecandidate, lastnodecandidate, checkselfoverlap=0, found;
	for(i=0;i<Nchains;i++) in_cluster[i]=0;
	in_cluster[seed]=1;
    double_triple axis_vector, sep;
	double angle=(2.0*rand_double-1.0)*max_rotate;
	declare_matrix(double, forward_matrix, 3, 3);
	declare_matrix(double, backward_matrix, 3, 3);
	declare_array_nozero(double_triple, mycom, Nchains);
    if(move3d==0){
		axis_vector.z=1;
		axis_vector.x=0;
		axis_vector.y=0;
    }
	else axis_vector=rand_unit_sphere();
	forward_and_backward_matrix(angle, axis_vector, forward_matrix, backward_matrix);
    
	double minbox=box_dimension.x;
	if(box_dimension.y<minbox) minbox=box_dimension.y;
	if(box_dimension.z<minbox) minbox=box_dimension.z;
	
	index=seed;
    
	exterior_member[0]=index;
	firstnode=monomerid[index][0].backbone;
	lastnode=monomerid[index][chainlength[index]-1].sidechain;
	if(coordarray[lastnode].nodetype<0) lastnode=monomerid[index][chainlength[index]-1].backbone;						//	in case long terminus
	
	checkrflag=calccomcheckr(&(mycom[index]), firstnode, lastnode, coordarray, box_dimension, pmaxr);										//	self-consistently, make sure no site is too far from com or could defeat selfoverlap check
	if(checkrflag==1){
		free_matrix(forward_matrix, 3);
		free_matrix(backward_matrix, 3);
		free(mycom);
        //printf("\tout, maxr\n");
		return;																										
	}
	
	for(i=firstnode;i<=lastnode;i++){
 		newcoordarray[i]=coordarray[i];
		reversecoordarray[i]=coordarray[i];
		sep=subtract_double_triple(coordarray[i].r, mycom[index]);
		recenter_double_triple(&sep, box_dimension);
		newcoordarray[i].r=add_double_triple(mycom[seed], rotate_by_matrix(sep, forward_matrix));
		reversecoordarray[i].r=add_double_triple(mycom[seed], rotate_by_matrix(sep, backward_matrix));
		newcoordarray[i].n=rotate_by_matrix(coordarray[i].n, forward_matrix);
		reversecoordarray[i].n=rotate_by_matrix(coordarray[i].n, backward_matrix);
		if(my_nonbonded_params.solvationparams.interface==1){
            newcoordarray[i].height=height_relative_to_interface_recreatecell(newcoordarray[i], my_nonbonded_params.solvationparams, i, (*pmeshlink), box_dimension, meshpositions, mapmeshtoneighbors);
            reversecoordarray[i].height=height_relative_to_interface_recreatecell(reversecoordarray[i], my_nonbonded_params.solvationparams, i, (*pmeshlink), box_dimension, meshpositions, mapmeshtoneighbors);
        }
	}
	while(cluster_complete==0){
		number_unchecked_neighbors=number_neighbors[index];
		while((number_unchecked_neighbors>0)&&(cluster_complete==0)){
			neighbor_rank=rand_int(number_unchecked_neighbors);						//	choose randomly from neighbors
			candidate=neighbor[index][neighbor_rank];
			for(i=neighbor_rank;i<number_unchecked_neighbors-1;i++){				//	shift neighbor list to close gap
				neighbor[index][i]=neighbor[index][i+1];
			}
			neighbor[index][i]=candidate;
			number_unchecked_neighbors--;
			if(in_cluster[candidate]==0){											//	not already in cluster
				individual_energy=collective_energy=0;				
				firstnodecandidate=monomerid[candidate][0].backbone;
				lastnodecandidate=monomerid[candidate][chainlength[index]-1].sidechain;
				if(coordarray[lastnodecandidate].nodetype<0) lastnodecandidate=monomerid[candidate][chainlength[index]-1].backbone;
				for(i=firstnode;i<=lastnode;i++){
					individual_energy+=calc_hard_energy_cellstruct_givenpolymer_recreatecell((*plinkset).shortrange, newcoordarray[i], coordarray, box_dimension, my_nonbonded_params, candidate);
					collective_energy+=calc_hard_energy_cellstruct_givenpolymer((*plinkset).shortrange.core.number_neighbors[i], (*plinkset).shortrange.core.neighbor[i], coordarray[i], coordarray, box_dimension, my_nonbonded_params, candidate);
					if(newcoordarray[i].nodetype==1){
						individual_energy+=calc_phenyl_energy_cellstruct_givenpolymer_recreatecell((*plinkset).phenylphenyl, newcoordarray[i], coordarray, box_dimension, my_nonbonded_params, candidate);
						collective_energy+=calc_phenyl_energy_cellstruct_givenpolymer((*plinkset).phenylphenyl.core.number_neighbors[i], (*plinkset).phenylphenyl.core.neighbor[i], coordarray[i], coordarray, box_dimension, my_nonbonded_params, candidate);
					}
					else if((newcoordarray[i].nodetype==2)||(newcoordarray[i].nodetype==3)){
						individual_energy+=calc_charged_energy_cellstruct_givenpolymer_recreatecell((*plinkset).chargedcharged, newcoordarray[i], coordarray, box_dimension, my_nonbonded_params, candidate);
						collective_energy+=calc_charged_energy_cellstruct_givenpolymer((*plinkset).chargedcharged.core.number_neighbors[i], (*plinkset).chargedcharged.core.neighbor[i], coordarray[i], coordarray, box_dimension, my_nonbonded_params, candidate);
					}
				}
                if(collective_energy==0){
                    
                    //  remove candidate from index's neighbor list (currently at position number_unchecked_neighbors
                    
                    for(i=number_unchecked_neighbors;i<number_neighbors[index]-1;i++){
                        neighbor[index][i]=neighbor[index][i+1];
                    }
                    number_neighbors[index]--;
                    
                    //  remove index from candidate's neighbor list
                    
                    found=0;
                    for(i=0;i<number_neighbors[candidate];i++){
                        if(neighbor[candidate][i]==index){
                            for(j=i;j<number_neighbors[candidate]-1;j++){
                                neighbor[candidate][j]=neighbor[candidate][j+1];
                            }
                            number_neighbors[candidate]--;
                            found=1;
                            break;
                        }
                    }
                    if(found==0){
                        printf("reverse map not found!!\n");
                        exit(1);
                    }
				}
				else{           
                    probability=1.0-exp(-(individual_energy-collective_energy)/temperature);
                    if(probability>0){
                        if((rand_double<probability)||(probability==1)){						//	pre-link
                            reverse_individual_energy=0;
                            for(i=firstnode;i<=lastnode;i++){
                                reverse_individual_energy+=calc_hard_energy_cellstruct_givenpolymer_recreatecell((*plinkset).shortrange, reversecoordarray[i], coordarray, box_dimension, my_nonbonded_params, candidate);
                                if(newcoordarray[i].nodetype==1){
                                    reverse_individual_energy+=calc_phenyl_energy_cellstruct_givenpolymer_recreatecell((*plinkset).phenylphenyl, reversecoordarray[i], coordarray, box_dimension, my_nonbonded_params, candidate);
                                }
                                else if((newcoordarray[i].nodetype==2)||(newcoordarray[i].nodetype==3)){
                                    reverse_individual_energy+=calc_charged_energy_cellstruct_givenpolymer_recreatecell((*plinkset).chargedcharged, reversecoordarray[i], coordarray, box_dimension, my_nonbonded_params, candidate);
                                }
                            }
                            reverse_probability=1.0-exp(-(reverse_individual_energy-collective_energy)/temperature);
                            if((reverse_probability>probability)||(rand_double<reverse_probability/probability)){		//	form link
                                cluster_number++;
                                exterior_member[number_exterior_members]=candidate;
                                number_exterior_members++;
                                in_cluster[candidate]=1;
                                checkrflag=calccomcheckr(&(mycom[candidate]), firstnodecandidate, lastnodecandidate, coordarray, box_dimension, pmaxr);
                                if(checkrflag==1){
                                    free_matrix(forward_matrix, 3);
                                    free_matrix(backward_matrix, 3);
                                    free(mycom);
                                    return;
                                }							
                                mycom[candidate]=subtract_double_triple(mycom[candidate], mycom[index]);
                                recenter_double_triple(&(mycom[candidate]), box_dimension);
                                mycom[candidate]=add_double_triple(mycom[candidate], mycom[index]);
                                for(i=firstnodecandidate;i<=lastnodecandidate;i++){
                                    sep=subtract_double_triple(coordarray[i].r, mycom[candidate]);
                                    recenter_double_triple(&sep, box_dimension);
                                    sep=subtract_double_triple(add_double_triple(sep, mycom[candidate]), mycom[seed]);
                                    newcoordarray[i].r=add_double_triple(mycom[seed], rotate_by_matrix(sep, forward_matrix));
                                    reversecoordarray[i].r=add_double_triple(mycom[seed], rotate_by_matrix(sep, backward_matrix));
                                    newcoordarray[i].n=rotate_by_matrix(coordarray[i].n, forward_matrix);
                                    reversecoordarray[i].n=rotate_by_matrix(coordarray[i].n, backward_matrix);
                                    if(my_nonbonded_params.solvationparams.interface==1){
                                        newcoordarray[i].height=height_relative_to_interface_recreatecell(newcoordarray[i], my_nonbonded_params.solvationparams, i, (*pmeshlink), box_dimension, meshpositions, mapmeshtoneighbors);
                                        reversecoordarray[i].height=height_relative_to_interface_recreatecell(reversecoordarray[i], my_nonbonded_params.solvationparams, i, (*pmeshlink), box_dimension, meshpositions, mapmeshtoneighbors);
                                    }
                                }
                                if(checkselfoverlap==0){
                                    if(norm(mycom[candidate])+2*(*pmaxr)+sqrt(my_nonbonded_params.cutoff2)>minbox) checkselfoverlap=1;
                                }
                            }
                            else{
                                frustrated_link[frustrated_link_index].x=index;				//	mark LINK as frustrated
                                frustrated_link[frustrated_link_index].y=candidate;
                                frustrated_link_index++;
                                if(frustrated_link_index==max_frustrated_links){
                                    printf("max_frustrated_links=%i not big enough.\n", max_frustrated_links);
                                    exit(1);
                                }
                            }
                        }
                    }
				}
			}
		}
		interior_member[number_interior_members]=index;								//	searched all neighbors; shift index from exterior to interior
		number_interior_members++;
		for(i=exterior_member_rank;i<number_exterior_members-1;i++){
			exterior_member[i]=exterior_member[i+1];
		}
		number_exterior_members--;
		if(number_exterior_members==0) cluster_complete=1;
		else{																		//	assign new cluster member to check neighbors of
			exterior_member_rank=rand_int(number_exterior_members);
			index=exterior_member[exterior_member_rank];
			firstnode=monomerid[index][0].backbone;
			lastnode=monomerid[index][chainlength[index]-1].sidechain;
			if(coordarray[lastnode].nodetype<0) lastnode=monomerid[index][chainlength[index]-1].backbone;
		}
	}
	for(i=0;i<frustrated_link_index;i++){
		if((in_cluster[frustrated_link[i].x]==0)||(in_cluster[frustrated_link[i].y]==0)){
			free_matrix(forward_matrix, 3);
			free_matrix(backward_matrix, 3);
			free(mycom);
			return;																	//	if there are any frustrated LINKS external to the cluster, reject
		}
	}
	
	//	Acceptance prob is min(1, x), where x is a product over pairs {i, j}
	//	where i is in cluster and j is outside
	//	energy either
	//		initially overlapping (positive) and finally noninteracting (0)
	//		or initially noninteracting (0) and finally overlapping (positive)
	//	Note that if I use D(C) not equal to 1, then that becomes the acceptance probability. But for now I am putting the size dependence in the distribution of max_number.
    
	//	BECAUSE OF POSSIBLE LARGE DISPLACEMENTS IN ROTATIONS, NEED TO ACCOUNT FOR POSSIBILITY THAT FINALLY POSITIVE INTERACTING PAIRS WERE NOT INITIALLY NEIGHBORS
    
	double sum_energy=0, initial_energy, final_energy;
	declare_array(int, neighborofindex, Nchains);
	declare_array(double, initiallynoninteractingenergiesbypolymer, Nchains);
	for(l=0;l<number_interior_members;l++){
		index=interior_member[l];
		firstnode=monomerid[index][0].backbone;
		lastnode=monomerid[index][chainlength[index]-1].sidechain;
		if(coordarray[lastnode].nodetype<0) lastnode=monomerid[index][chainlength[index]-1].backbone;
        
		//	first, do regular sum over neighbors of polymer #index
		
		for(m=0;m<Nchains;m++){
			neighborofindex[m]=0;
			initiallynoninteractingenergiesbypolymer[m]=0;
		}
		neighborofindex[index]=1;		//	self
		for(m=0;m<number_neighbors[index];m++){
			candidate=neighbor[index][m];
			neighborofindex[candidate]=1;
			
			//	calculating contribution from initial neighbors index and candidate
			
			if(in_cluster[candidate]==0){
				firstnodecandidate=monomerid[candidate][0].backbone;
				lastnodecandidate=monomerid[candidate][chainlength[index]-1].sidechain;
				if(coordarray[lastnode].nodetype<0) lastnodecandidate=monomerid[candidate][chainlength[index]-1].backbone;
				initial_energy=0;
				
				for(i=firstnode;i<=lastnode;i++){
					for(j=firstnodecandidate;j<=lastnodecandidate;j++){
						initial_energy+=calc_energy_pair(coordarray[i], coordarray[j], my_nonbonded_params, box_dimension);
					}
				}
				
				if(initial_energy>0){
					final_energy=0;
					for(i=firstnode;i<=lastnode;i++){
						final_energy+=calc_hard_energy_cellstruct_givenpolymer((*plinkset).shortrange.core.number_neighbors[i], (*plinkset).shortrange.core.neighbor[i], newcoordarray[i], coordarray, box_dimension, my_nonbonded_params, candidate);
						if(newcoordarray[i].nodetype==1){
							final_energy+=calc_phenyl_energy_cellstruct_givenpolymer((*plinkset).phenylphenyl.core.number_neighbors[i], (*plinkset).phenylphenyl.core.neighbor[i], newcoordarray[i], coordarray, box_dimension, my_nonbonded_params, candidate);
						}
						else if((newcoordarray[i].nodetype==2)||(newcoordarray[i].nodetype==3)){
							final_energy+=calc_charged_energy_cellstruct_givenpolymer((*plinkset).chargedcharged.core.number_neighbors[i], (*plinkset).chargedcharged.core.neighbor[i], newcoordarray[i], coordarray, box_dimension, my_nonbonded_params, candidate);
						}
					}
					
					if(final_energy==0) sum_energy-=initial_energy;
				}
				else if(initial_energy==0){
					final_energy=0;
					for(i=firstnode;i<=lastnode;i++){
						final_energy+=calc_hard_energy_cellstruct_givenpolymer((*plinkset).shortrange.core.number_neighbors[i], (*plinkset).shortrange.core.neighbor[i], newcoordarray[i], coordarray, box_dimension, my_nonbonded_params, candidate);
						if(newcoordarray[i].nodetype==1){
							final_energy+=calc_phenyl_energy_cellstruct_givenpolymer((*plinkset).phenylphenyl.core.number_neighbors[i], (*plinkset).phenylphenyl.core.neighbor[i], newcoordarray[i], coordarray, box_dimension, my_nonbonded_params, candidate);
						}
						else if((newcoordarray[i].nodetype==2)||(newcoordarray[i].nodetype==3)){
							final_energy+=calc_charged_energy_cellstruct_givenpolymer((*plinkset).chargedcharged.core.number_neighbors[i], (*plinkset).chargedcharged.core.neighbor[i], newcoordarray[i], coordarray, box_dimension, my_nonbonded_params, candidate);
						}
					}
					if(final_energy>0) sum_energy+=final_energy;
				}
			}
		}
		
		//	next, add up energies for any polymers that were not initially neighbors
		//	bin the energies by polymers, because only positive pairwise polymer-polymer energies contribute to acceptance criterion
		
		for(i=firstnode;i<lastnode;i++){
			fmod_double_triple(&(newcoordarray[i].r), box_dimension);                   //  need to apply boundary conditions because I look up cells, not just neighbors, here
            add_hardenergy_if_polymer_not_in_list(newcoordarray[i], (*plinkset).shortrange, coordarray, neighborofindex, initiallynoninteractingenergiesbypolymer, my_nonbonded_params, box_dimension);
			if(coordarray[i].nodetype==1) add_phenylphenylenergy_if_polymer_not_in_list(newcoordarray[i], (*plinkset).phenylphenyl, coordarray, neighborofindex, initiallynoninteractingenergiesbypolymer, my_nonbonded_params, box_dimension);
			else if((coordarray[i].nodetype==2)||(coordarray[i].nodetype==3)){
                add_chargedchargedenergy_if_polymer_not_in_list(newcoordarray[i], (*plinkset).chargedcharged, coordarray, neighborofindex, initiallynoninteractingenergiesbypolymer, my_nonbonded_params, box_dimension);
            }
		}
		for(m=0;m<Nchains;m++){
			if(initiallynoninteractingenergiesbypolymer[m]>0) sum_energy+=initiallynoninteractingenergiesbypolymer[m];
		}
	}
	free(neighborofindex);
	free(initiallynoninteractingenergiesbypolymer);
	
	//	external potential
	
    for(l=0;l<number_interior_members;l++){
		index=interior_member[l];
		firstnode=monomerid[index][0].backbone;
		lastnode=monomerid[index][chainlength[index]-1].sidechain;
		if(coordarray[lastnode].nodetype<0) lastnode=monomerid[index][chainlength[index]-1].backbone;
		for(i=firstnode;i<=lastnode;i++){
            sum_energy-=calc_onebody_energy(coordarray[i], chooseonebody(newcoordarray[i].nodetype, my_nonbonded_params), my_nonbonded_params.solvationparams);
			sum_energy+=calc_onebody_energy(newcoordarray[i], chooseonebody(newcoordarray[i].nodetype, my_nonbonded_params), my_nonbonded_params.solvationparams);
            
            //  Treat change in phenyl-phenyl interaction due to change in height as external potential
            
            if(my_nonbonded_params.solvationparams.interface==1){
                sum_energy-=calc_phenyl_energy_cellstruct_polymerlist_singlecount((*plinkset).phenylphenyl.core.number_neighbors[i], (*plinkset).phenylphenyl.core.neighbor[i], coordarray[i], coordarray, box_dimension, my_nonbonded_params, in_cluster);
                sum_energy+=calc_phenyl_energy_cellstruct_polymerlist_singlecount((*plinkset).phenylphenyl.core.number_neighbors[i], (*plinkset).phenylphenyl.core.neighbor[i], newcoordarray[i], newcoordarray, box_dimension, my_nonbonded_params, in_cluster);
            }
		}
	}
    
	probability=exp(-sum_energy/temperature);
	if((probability<1.0)&&(rand_double>probability)){			//	reject
		free_matrix(forward_matrix, 3);
		free_matrix(backward_matrix, 3);
		free(mycom);
        return;
	}
    
	//	otherwise accept	
	    
	//	update coordinates and cells
    
	for(l=0;l<number_interior_members;l++){
		index=interior_member[l];
		firstnode=monomerid[index][0].backbone;
		lastnode=monomerid[index][chainlength[index]-1].sidechain;
		if(coordarray[lastnode].nodetype<0) lastnode=monomerid[index][chainlength[index]-1].backbone;
		for(i=firstnode;i<=lastnode;i++){
			fmod_double_triple(&(newcoordarray[i].r), box_dimension);
            if(my_nonbonded_params.solvationparams.interface==1)updatecell_asymmetric(newcoordarray[i].r, &(*pmeshlink), i);
			updatecell(newcoordarray[i].r, &((*plinkset).shortrange), i);
			if(newcoordarray[i].nodetype==1) updatecell(newcoordarray[i].r, &((*plinkset).phenylphenyl), i);
			else if((newcoordarray[i].nodetype==2)||(newcoordarray[i].nodetype==3)) updatecell(newcoordarray[i].r, &((*plinkset).chargedcharged), i);
			coordarray[i]=newcoordarray[i];
		}
	}
    
	//	update site-site neighbor list
    	
	for(l=0;l<number_interior_members;l++){
		index=interior_member[l];
		firstnode=monomerid[index][0].backbone;
		lastnode=monomerid[index][chainlength[index]-1].sidechain;
		if(coordarray[lastnode].nodetype<0) lastnode=monomerid[index][chainlength[index]-1].backbone;
		for(i=firstnode;i<=lastnode;i++){
            if(my_nonbonded_params.solvationparams.interface==1) updateneighborlist_asymmetric(&(*pmeshlink), i, coordarray[i].r, meshpositions, box_dimension);
			updateneighborlist(&((*plinkset).shortrange), i, coordarray, box_dimension);
			if(newcoordarray[i].nodetype==1) updateneighborlist_pairenergies_phenyl(&((*plinkset).phenylphenyl), i, coordarray, box_dimension, my_nonbonded_params);
			else if((newcoordarray[i].nodetype==2)||(newcoordarray[i].nodetype==3)) updateneighborlist_pairenergies_charged(&((*plinkset).chargedcharged), i, coordarray, box_dimension, my_nonbonded_params);
		}
	}
    
	//	Update polymer-polymer neighbor lists. First, wipe out neighbor lists starting and ending in the cluster
	
	for(i=0;i<number_interior_members;i++){
		index=interior_member[i];
		for(j=0;j<number_neighbors[index];j++){
			temp=neighbor[index][j];
			if(in_cluster[temp]==0){
				for(k=0;k<number_neighbors[temp];k++){
					if(neighbor[temp][k]==index) break;								//	find the reverse neighbor map; can I do this more efficiently?
				}
				if(k==number_neighbors[temp]){
					printf("Error: did not find reverse map.\n");
					exit(1);
				}
				number_neighbors[temp]--;											//	remove reverse neighbor map
				for(l=k;l<number_neighbors[temp];l++){
					neighbor[temp][l]=neighbor[temp][l+1];
				}
			}
		}
		number_neighbors[index]=0;                                                  //  remove all forward maps
	}
    
	//	Now update neighbor lists
    
	int target, targetpolymer;
	double rsq;
	double_triple r;
	for(l=0;l<number_interior_members;l++){
		index=interior_member[l];
		firstnode=monomerid[index][0].backbone;
		lastnode=monomerid[index][chainlength[index]-1].sidechain;
		if(coordarray[lastnode].nodetype<0) lastnode=monomerid[index][chainlength[index]-1].backbone;
		for(i=0;i<Nchains;i++) assigned[i]=0;
		assigned[index]=1;
        
        for(i=0;i<number_neighbors[index];i++) assigned[neighbor[index][i]]=1;
		
        for(i=firstnode;i<=lastnode;i++){
			
			//	neighbor if interacting (but don't need to check hard core interactions, since these infinite interactions would not have been allowed)
			
			if(coordarray[i].nodetype==1){
                for(j=0;j<(*plinkset).phenylphenyl.core.number_neighbors[i];j++){
					target=(*plinkset).phenylphenyl.core.neighbor[i][j];
					targetpolymer=coordarray[target].chainid;
					if(assigned[targetpolymer]==0){
						r=subtract_double_triple(coordarray[target].r, coordarray[i].r);
						rsq=dot_product(r, r);
						if(rsq<my_nonbonded_params.phenylcutoff2){				//	interacting
							assigned[targetpolymer]=1;
							neighbor[index][number_neighbors[index]]=targetpolymer;
							number_neighbors[index]++;
							if(number_neighbors[index]==max_neighbors) my_exit("max_neighbors not big enough (a)");
							if(in_cluster[targetpolymer]==0){					//	safe to add reverse neighbor link without double counting
                                
                                //  okay because in_cluster equivalent to interior_member, since there should be no exterior_members
                                
								neighbor[targetpolymer][number_neighbors[targetpolymer]]=index;
								number_neighbors[targetpolymer]++;
								if(number_neighbors[targetpolymer]==max_neighbors) my_exit("max_neighbors not big enough (b)");
							}
							else{
								found=0;
								for(temp=0;temp<number_neighbors[targetpolymer];temp++){			//	only add link if target->index doesn't already exist
									if(neighbor[targetpolymer][temp]==index) found=1;
								}
								if(found==0){
									neighbor[targetpolymer][number_neighbors[targetpolymer]]=index;
									number_neighbors[targetpolymer]++;
									if(number_neighbors[targetpolymer]==max_neighbors) my_exit("max_neighbors not big enough (c)");
								}
							}
						}
					}
				}
			}
			else if((coordarray[i].nodetype==2)||(coordarray[i].nodetype==3)){
                for(j=0;j<(*plinkset).chargedcharged.core.number_neighbors[i];j++){
					target=(*plinkset).chargedcharged.core.neighbor[i][j];
					targetpolymer=coordarray[target].chainid;
					if(assigned[targetpolymer]==0){
						r=subtract_double_triple(coordarray[target].r, coordarray[i].r);
						rsq=dot_product(r, r);
						if(rsq<my_nonbonded_params.cutoff2){				//	interacting
							assigned[targetpolymer]=1;
							neighbor[index][number_neighbors[index]]=targetpolymer;
							number_neighbors[index]++;
							if(number_neighbors[index]==max_neighbors) my_exit("max_neighbors not big enough (d)");
							if(in_cluster[targetpolymer]==0){					//	safe to add reverse neighbor link without double counting
								neighbor[targetpolymer][number_neighbors[targetpolymer]]=index;
								number_neighbors[targetpolymer]++;
 								if(number_neighbors[targetpolymer]==max_neighbors) my_exit("max_neighbors not big enough (e)");
							}
							else{
								found=0;
								for(temp=0;temp<number_neighbors[targetpolymer];temp++){			//	only add link if target->index doesn't already exist
									if(neighbor[targetpolymer][temp]==index) found=1;
								}
								if(found==0){
									neighbor[targetpolymer][number_neighbors[targetpolymer]]=index;
									number_neighbors[targetpolymer]++;
 									if(number_neighbors[targetpolymer]==max_neighbors) my_exit("max_neighbors not big enough (f)");
								}
							}
						}
					}
				}
			}
		}
	}
	free_matrix(forward_matrix, 3);
	free_matrix(backward_matrix, 3);
	free(mycom);
}

double calc_hard_energy_cellstruct_givenpolymer_recreatecell(linkedlistfull mylinkedlistfull, coord coord_new, coord *coordarray, double_triple box_dimension, nonbonded_params my_nonbonded_params, int polymer){
    fmod_double_triple(&(coord_new.r), box_dimension);
	int_triple newcell=cellofpos(coord_new.r, mylinkedlistfull.cellwidth), imagecell;
	int trialneighbor, a, b, c;
    double result=0;
	for(a=-1;a<=1;a++){
		if((newcell.x==0)&&(a==-1)) imagecell.x=mylinkedlistfull.core.cellsperside.x-1;
		else if((newcell.x==mylinkedlistfull.core.cellsperside.x-1)&&(a==1)) imagecell.x=0;
		else imagecell.x=newcell.x+a;
		for(b=-1;b<=1;b++){
			if((newcell.y==0)&&(b==-1)) imagecell.y=mylinkedlistfull.core.cellsperside.y-1;
			else if((newcell.y==mylinkedlistfull.core.cellsperside.y-1)&&(b==1)) imagecell.y=0;
			else imagecell.y=newcell.y+b;
			for(c=-1;c<=1;c++){
				if((newcell.z==0)&&(c==-1)) imagecell.z=mylinkedlistfull.core.cellsperside.z-1;
				else if((newcell.z==mylinkedlistfull.core.cellsperside.z-1)&&(c==1)) imagecell.z=0;
				else imagecell.z=newcell.z+c;
				trialneighbor=mylinkedlistfull.core.head[imagecell.x][imagecell.y][imagecell.z];
				while(trialneighbor>=0){
					if(coordarray[trialneighbor].chainid==polymer){
                        
						if(coord_new.nodetype==0){
							if(coordarray[trialneighbor].nodetype==0){
								result+=hard_energy_sumradii(coord_new.r, coordarray[trialneighbor].r, 2.*my_nonbonded_params.rhard0, box_dimension);
							}
							else if(coordarray[trialneighbor].nodetype==1){
								result+=hard_energy_sumradii(coord_new.r, coordarray[trialneighbor].r, my_nonbonded_params.rhard0+my_nonbonded_params.rhard1, box_dimension);
							}
							else if(coordarray[trialneighbor].nodetype==2){
								result+=hard_energy_sumradii(coord_new.r, coordarray[trialneighbor].r, my_nonbonded_params.rhard0+my_nonbonded_params.p2.rhard, box_dimension);
							}
							else if(coordarray[trialneighbor].nodetype==3){
								result+=hard_energy_sumradii(coord_new.r, coordarray[trialneighbor].r, my_nonbonded_params.rhard0+my_nonbonded_params.p3.rhard, box_dimension);
							}
						}
						else if(coord_new.nodetype==1){
							if(coordarray[trialneighbor].nodetype==0){
								result+=hard_energy_sumradii(coord_new.r, coordarray[trialneighbor].r, my_nonbonded_params.rhard1+my_nonbonded_params.rhard0, box_dimension);
							}
							else if(coordarray[trialneighbor].nodetype==2){
								result+=hard_energy_sumradii(coord_new.r, coordarray[trialneighbor].r, my_nonbonded_params.rhard1+my_nonbonded_params.p2.rhard, box_dimension);
							}
							else if(coordarray[trialneighbor].nodetype==3){
								result+=hard_energy_sumradii(coord_new.r, coordarray[trialneighbor].r, my_nonbonded_params.rhard1+my_nonbonded_params.p3.rhard, box_dimension);
							}
						}
						else if(coord_new.nodetype==2){
							if(coordarray[trialneighbor].nodetype==0){
								result+=hard_energy_sumradii(coord_new.r, coordarray[trialneighbor].r, my_nonbonded_params.p2.rhard+my_nonbonded_params.rhard0, box_dimension);
							}
							else if(coordarray[trialneighbor].nodetype==1){
								result+=hard_energy_sumradii(coord_new.r, coordarray[trialneighbor].r, my_nonbonded_params.p2.rhard+my_nonbonded_params.rhard1, box_dimension);
							}
						}
						else if(coord_new.nodetype==3){
							if(coordarray[trialneighbor].nodetype==0){
								result+=hard_energy_sumradii(coord_new.r, coordarray[trialneighbor].r, my_nonbonded_params.p3.rhard+my_nonbonded_params.rhard0, box_dimension);
							}
							else if(coordarray[trialneighbor].nodetype==1){
								result+=hard_energy_sumradii(coord_new.r, coordarray[trialneighbor].r, my_nonbonded_params.p3.rhard+my_nonbonded_params.rhard1, box_dimension);
							}
						}
						
					}
					trialneighbor=mylinkedlistfull.core.list[trialneighbor];
                    
				}
			}
		}
	}
    return result;
}

double calc_phenyl_energy_cellstruct_givenpolymer_recreatecell(linkedlistfull mylinkedlistfull, coord coord_new, coord *coordarray, double_triple box_dimension, nonbonded_params my_nonbonded_params, int polymer){
    fmod_double_triple(&(coord_new.r), box_dimension);
	int_triple newcell=cellofpos(coord_new.r, mylinkedlistfull.cellwidth), imagecell;
	int trialneighbor, a, b, c;
    double result=0;
	for(a=-1;a<=1;a++){
		if((newcell.x==0)&&(a==-1)) imagecell.x=mylinkedlistfull.core.cellsperside.x-1;
		else if((newcell.x==mylinkedlistfull.core.cellsperside.x-1)&&(a==1)) imagecell.x=0;
		else imagecell.x=newcell.x+a;
		for(b=-1;b<=1;b++){
			if((newcell.y==0)&&(b==-1)) imagecell.y=mylinkedlistfull.core.cellsperside.y-1;
			else if((newcell.y==mylinkedlistfull.core.cellsperside.y-1)&&(b==1)) imagecell.y=0;
			else imagecell.y=newcell.y+b;
			for(c=-1;c<=1;c++){
				if((newcell.z==0)&&(c==-1)) imagecell.z=mylinkedlistfull.core.cellsperside.z-1;
				else if((newcell.z==mylinkedlistfull.core.cellsperside.z-1)&&(c==1)) imagecell.z=0;
				else imagecell.z=newcell.z+c;
				trialneighbor=mylinkedlistfull.core.head[imagecell.x][imagecell.y][imagecell.z];
				while(trialneighbor>=0){
					if(coordarray[trialneighbor].chainid==polymer){
						result+=phenylphenyl_energy(coord_new, coordarray[trialneighbor], my_nonbonded_params.p11vac, box_dimension, my_nonbonded_params.phenylcutoff2, my_nonbonded_params.solvationparams);
					}
					trialneighbor=mylinkedlistfull.core.list[trialneighbor];
					
				}
			}
		}
	}
    return result;
}

double calc_charged_energy_cellstruct_givenpolymer_recreatecell(linkedlistfull mylinkedlistfull, coord coord_new, coord *coordarray, double_triple box_dimension, nonbonded_params my_nonbonded_params, int polymer){
    fmod_double_triple(&(coord_new.r), box_dimension);
	int_triple newcell=cellofpos(coord_new.r, mylinkedlistfull.cellwidth), imagecell;
	int trialneighbor, a, b, c;
    double result=0;
	for(a=-1;a<=1;a++){
		if((newcell.x==0)&&(a==-1)) imagecell.x=mylinkedlistfull.core.cellsperside.x-1;
		else if((newcell.x==mylinkedlistfull.core.cellsperside.x-1)&&(a==1)) imagecell.x=0;
		else imagecell.x=newcell.x+a;
		for(b=-1;b<=1;b++){
			if((newcell.y==0)&&(b==-1)) imagecell.y=mylinkedlistfull.core.cellsperside.y-1;
			else if((newcell.y==mylinkedlistfull.core.cellsperside.y-1)&&(b==1)) imagecell.y=0;
			else imagecell.y=newcell.y+b;
			for(c=-1;c<=1;c++){
				if((newcell.z==0)&&(c==-1)) imagecell.z=mylinkedlistfull.core.cellsperside.z-1;
				else if((newcell.z==mylinkedlistfull.core.cellsperside.z-1)&&(c==1)) imagecell.z=0;
				else imagecell.z=newcell.z+c;
				trialneighbor=mylinkedlistfull.core.head[imagecell.x][imagecell.y][imagecell.z];
				while(trialneighbor>=0){
					if(coordarray[trialneighbor].chainid==polymer){
						if(coord_new.nodetype==2){
							if(coordarray[trialneighbor].nodetype==2){
								result+=electrostatic_energy(coord_new, coordarray[trialneighbor], my_nonbonded_params.solvationparams, my_nonbonded_params.p2, my_nonbonded_params.p2, box_dimension, my_nonbonded_params.cutoff2);
							}
							else if(coordarray[trialneighbor].nodetype==3){
								result+=electrostatic_energy(coord_new, coordarray[trialneighbor], my_nonbonded_params.solvationparams, my_nonbonded_params.p2, my_nonbonded_params.p3, box_dimension, my_nonbonded_params.cutoff2);
							}
						}
						else if(coord_new.nodetype==3){
							if(coordarray[trialneighbor].nodetype==2){
								result+=electrostatic_energy(coord_new, coordarray[trialneighbor], my_nonbonded_params.solvationparams, my_nonbonded_params.p3, my_nonbonded_params.p2, box_dimension, my_nonbonded_params.cutoff2);
							}
							else if(coordarray[trialneighbor].nodetype==3){
								result+=electrostatic_energy(coord_new, coordarray[trialneighbor], my_nonbonded_params.solvationparams, my_nonbonded_params.p3, my_nonbonded_params.p3, box_dimension, my_nonbonded_params.cutoff2);
							}
						}
					}
                    trialneighbor=mylinkedlistfull.core.list[trialneighbor];
				}
			}
		}
	}
    return result;
}

void add_hardenergy_if_polymer_not_in_list(coord newcoord, linkedlistfull mylinkedlistfull, coord *coordarray, int *list, double *energy, nonbonded_params my_nonbonded_params, double_triple box_dimension){
	int_triple newcell=cellofpos(newcoord.r, mylinkedlistfull.cellwidth), imagecell;
	int trialneighbor, trialpolymer, a, b, c;
	for(a=-1;a<=1;a++){
		if((newcell.x==0)&&(a==-1)) imagecell.x=mylinkedlistfull.core.cellsperside.x-1;
		else if((newcell.x==mylinkedlistfull.core.cellsperside.x-1)&&(a==1)) imagecell.x=0;
		else imagecell.x=newcell.x+a;
		for(b=-1;b<=1;b++){
			if((newcell.y==0)&&(b==-1)) imagecell.y=mylinkedlistfull.core.cellsperside.y-1;
			else if((newcell.y==mylinkedlistfull.core.cellsperside.y-1)&&(b==1)) imagecell.y=0;
			else imagecell.y=newcell.y+b;
			for(c=-1;c<=1;c++){
				if((newcell.z==0)&&(c==-1)) imagecell.z=mylinkedlistfull.core.cellsperside.z-1;
				else if((newcell.z==mylinkedlistfull.core.cellsperside.z-1)&&(c==1)) imagecell.z=0;
				else imagecell.z=newcell.z+c;
				trialneighbor=mylinkedlistfull.core.head[imagecell.x][imagecell.y][imagecell.z];
				while(trialneighbor>=0){
					trialpolymer=coordarray[trialneighbor].chainid;
					if(list[trialpolymer]==0){
                        
						if(newcoord.nodetype==0){
							if(coordarray[trialneighbor].nodetype==0){
								energy[trialpolymer]+=hard_energy_sumradii(newcoord.r, coordarray[trialneighbor].r, 2.*my_nonbonded_params.rhard0, box_dimension);
							}
							else if(coordarray[trialneighbor].nodetype==1){
								energy[trialpolymer]+=hard_energy_sumradii(newcoord.r, coordarray[trialneighbor].r, my_nonbonded_params.rhard0+my_nonbonded_params.rhard1, box_dimension);
							}
							else if(coordarray[trialneighbor].nodetype==2){
								energy[trialpolymer]+=hard_energy_sumradii(newcoord.r, coordarray[trialneighbor].r, my_nonbonded_params.rhard0+my_nonbonded_params.p2.rhard, box_dimension);
							}
							else if(coordarray[trialneighbor].nodetype==3){
								energy[trialpolymer]+=hard_energy_sumradii(newcoord.r, coordarray[trialneighbor].r, my_nonbonded_params.rhard0+my_nonbonded_params.p3.rhard, box_dimension);
							}
						}
						else if(newcoord.nodetype==1){
							if(coordarray[trialneighbor].nodetype==0){
								energy[trialpolymer]+=hard_energy_sumradii(newcoord.r, coordarray[trialneighbor].r, my_nonbonded_params.rhard1+my_nonbonded_params.rhard0, box_dimension);
							}
							else if(coordarray[trialneighbor].nodetype==2){
								energy[trialpolymer]+=hard_energy_sumradii(newcoord.r, coordarray[trialneighbor].r, my_nonbonded_params.rhard1+my_nonbonded_params.p2.rhard, box_dimension);
							}
							else if(coordarray[trialneighbor].nodetype==3){
								energy[trialpolymer]+=hard_energy_sumradii(newcoord.r, coordarray[trialneighbor].r, my_nonbonded_params.rhard1+my_nonbonded_params.p3.rhard, box_dimension);
							}
						}
						else if(newcoord.nodetype==2){
							if(coordarray[trialneighbor].nodetype==0){
								energy[trialpolymer]+=hard_energy_sumradii(newcoord.r, coordarray[trialneighbor].r, my_nonbonded_params.p2.rhard+my_nonbonded_params.rhard0, box_dimension);
							}
							else if(coordarray[trialneighbor].nodetype==1){
								energy[trialpolymer]+=hard_energy_sumradii(newcoord.r, coordarray[trialneighbor].r, my_nonbonded_params.p2.rhard+my_nonbonded_params.rhard1, box_dimension);
							}
						}
						else if(newcoord.nodetype==3){
							if(coordarray[trialneighbor].nodetype==0){
								energy[trialpolymer]+=hard_energy_sumradii(newcoord.r, coordarray[trialneighbor].r, my_nonbonded_params.p3.rhard+my_nonbonded_params.rhard0, box_dimension);
							}
							else if(coordarray[trialneighbor].nodetype==1){
								energy[trialpolymer]+=hard_energy_sumradii(newcoord.r, coordarray[trialneighbor].r, my_nonbonded_params.p3.rhard+my_nonbonded_params.rhard1, box_dimension);
							}
						}
						
					}
					trialneighbor=mylinkedlistfull.core.list[trialneighbor];
                    
				}
			}
		}
	}
}

void add_phenylphenylenergy_if_polymer_not_in_list(coord newcoord, linkedlistfull mylinkedlistfull, coord *coordarray, int *list, double *energy, nonbonded_params my_nonbonded_params, double_triple box_dimension){
	int_triple newcell=cellofpos(newcoord.r, mylinkedlistfull.cellwidth), imagecell;
	int trialneighbor, trialpolymer, a, b, c;
	for(a=-1;a<=1;a++){
		if((newcell.x==0)&&(a==-1)) imagecell.x=mylinkedlistfull.core.cellsperside.x-1;
		else if((newcell.x==mylinkedlistfull.core.cellsperside.x-1)&&(a==1)) imagecell.x=0;
		else imagecell.x=newcell.x+a;
		for(b=-1;b<=1;b++){
			if((newcell.y==0)&&(b==-1)) imagecell.y=mylinkedlistfull.core.cellsperside.y-1;
			else if((newcell.y==mylinkedlistfull.core.cellsperside.y-1)&&(b==1)) imagecell.y=0;
			else imagecell.y=newcell.y+b;
			for(c=-1;c<=1;c++){
				if((newcell.z==0)&&(c==-1)) imagecell.z=mylinkedlistfull.core.cellsperside.z-1;
				else if((newcell.z==mylinkedlistfull.core.cellsperside.z-1)&&(c==1)) imagecell.z=0;
				else imagecell.z=newcell.z+c;
				trialneighbor=mylinkedlistfull.core.head[imagecell.x][imagecell.y][imagecell.z];
				while(trialneighbor>=0){
					trialpolymer=coordarray[trialneighbor].chainid;
					if(list[trialpolymer]==0){						
						energy[trialpolymer]+=phenylphenyl_energy(newcoord, coordarray[trialneighbor], my_nonbonded_params.p11vac, box_dimension, my_nonbonded_params.phenylcutoff2, my_nonbonded_params.solvationparams);
					}
					trialneighbor=mylinkedlistfull.core.list[trialneighbor];
					
				}
			}
		}
	}
}

double calc_phenyl_energy_cellstruct_polymerlist_singlecount(int numberneighbors, int *neighborlist, coord coord_new, coord *coordarray, double_triple box_dimension, nonbonded_params my_nonbonded_params, int *polymerlist){
	int i, index, selfpolymer=coord_new.chainid, otherpolymer;
	double result=0;
	for(i=0;i<numberneighbors;i++){
		index=neighborlist[i];
        otherpolymer=coordarray[index].chainid;
        if(((otherpolymer==selfpolymer)&&(coordarray[index].monomerid>coord_new.monomerid))||((polymerlist[otherpolymer]==1)&&(otherpolymer>selfpolymer))){
 			
			// should only be in linked list if both types are phenyl
			
			result+=phenylphenyl_energy(coord_new, coordarray[index], my_nonbonded_params.p11vac, box_dimension, my_nonbonded_params.phenylcutoff2, my_nonbonded_params.solvationparams);
		}
	}
	return result;
}

void add_chargedchargedenergy_if_polymer_not_in_list(coord newcoord, linkedlistfull mylinkedlistfull, coord *coordarray, int *list, double *energy, nonbonded_params my_nonbonded_params, double_triple box_dimension){
	int_triple newcell=cellofpos(newcoord.r, mylinkedlistfull.cellwidth), imagecell;
	int trialneighbor, trialpolymer, a, b, c;
	for(a=-1;a<=1;a++){
		if((newcell.x==0)&&(a==-1)) imagecell.x=mylinkedlistfull.core.cellsperside.x-1;
		else if((newcell.x==mylinkedlistfull.core.cellsperside.x-1)&&(a==1)) imagecell.x=0;
		else imagecell.x=newcell.x+a;
		for(b=-1;b<=1;b++){
			if((newcell.y==0)&&(b==-1)) imagecell.y=mylinkedlistfull.core.cellsperside.y-1;
			else if((newcell.y==mylinkedlistfull.core.cellsperside.y-1)&&(b==1)) imagecell.y=0;
			else imagecell.y=newcell.y+b;
			for(c=-1;c<=1;c++){
				if((newcell.z==0)&&(c==-1)) imagecell.z=mylinkedlistfull.core.cellsperside.z-1;
				else if((newcell.z==mylinkedlistfull.core.cellsperside.z-1)&&(c==1)) imagecell.z=0;
				else imagecell.z=newcell.z+c;
				trialneighbor=mylinkedlistfull.core.head[imagecell.x][imagecell.y][imagecell.z];
				while(trialneighbor>=0){
					trialpolymer=coordarray[trialneighbor].chainid;
					if(list[trialpolymer]==0){		
						if(newcoord.nodetype==2){
							if(coordarray[trialneighbor].nodetype==2){
								energy[trialpolymer]+=electrostatic_energy(newcoord, coordarray[trialneighbor], my_nonbonded_params.solvationparams, my_nonbonded_params.p2, my_nonbonded_params.p2, box_dimension, my_nonbonded_params.cutoff2);
							}
							else if(coordarray[trialneighbor].nodetype==3){
								energy[trialpolymer]+=electrostatic_energy(newcoord, coordarray[trialneighbor], my_nonbonded_params.solvationparams, my_nonbonded_params.p2, my_nonbonded_params.p3, box_dimension, my_nonbonded_params.cutoff2);
							}
						}
						else if(newcoord.nodetype==3){
							if(coordarray[trialneighbor].nodetype==2){
								energy[trialpolymer]+=electrostatic_energy(newcoord, coordarray[trialneighbor], my_nonbonded_params.solvationparams, my_nonbonded_params.p3, my_nonbonded_params.p2, box_dimension, my_nonbonded_params.cutoff2);
							}
							else if(coordarray[trialneighbor].nodetype==3){
								energy[trialpolymer]+=electrostatic_energy(newcoord, coordarray[trialneighbor], my_nonbonded_params.solvationparams, my_nonbonded_params.p3, my_nonbonded_params.p3, box_dimension, my_nonbonded_params.cutoff2);
							}
						}
					}
                    trialneighbor=mylinkedlistfull.core.list[trialneighbor];					
				}
			}
		}
	}
}

void mc_translate_wholepolymer(int Nchains, int *chainlength, coord *newcoordarray, coord *coordarray, double max_translate, double_triple box_dimension, linkedlistset *plinkset, bonded_params my_bonded_params, nonbonded_params my_nonbonded_params, double temperature, monomernodes **monomerid, linkedlistasymmetric *pmeshlink, double_triple *meshpositions, int **mapmeshtoneighbors){
	int movingpolymer=rand_int(Nchains), firstnode, lastnode, i;
	firstnode=monomerid[movingpolymer][0].backbone;
	lastnode=monomerid[movingpolymer][chainlength[movingpolymer]-1].sidechain;
	if(coordarray[lastnode].nodetype<0) lastnode=monomerid[movingpolymer][chainlength[movingpolymer]-1].backbone;
	double acceptance_prob, difference=0;
	double_triple shift=scalar_multiply_double_triple(rand_unit_cube(), (rand_int(2)*2-1)*max_translate);
	for(i=firstnode;i<=lastnode;i++){
		newcoordarray[i]=coordarray[i];
		newcoordarray[i].r=add_double_triple(coordarray[i].r, shift);
		if(my_nonbonded_params.solvationparams.interface==1) newcoordarray[i].height=height_relative_to_interface(i, newcoordarray[i], my_nonbonded_params.solvationparams, (*pmeshlink), box_dimension, meshpositions, mapmeshtoneighbors);
		difference+=calc_hard_energy_cellstruct_otherpolymer((*plinkset).shortrange.core.number_neighbors[i], (*plinkset).shortrange.core.neighbor[i], newcoordarray[i], coordarray, box_dimension, my_nonbonded_params);
		difference-=calc_hard_energy_cellstruct_otherpolymer((*plinkset).shortrange.core.number_neighbors[i], (*plinkset).shortrange.core.neighbor[i], coordarray[i], coordarray, box_dimension, my_nonbonded_params);
        if(newcoordarray[i].nodetype==1){
			difference+=calc_phenyl_energy_cellstruct_otherpolymer((*plinkset).phenylphenyl.core.number_neighbors[i], (*plinkset).phenylphenyl.core.neighbor[i], newcoordarray[i], coordarray, box_dimension, my_nonbonded_params);
			difference-=calc_phenyl_energy_cellstruct_otherpolymer((*plinkset).phenylphenyl.core.number_neighbors[i], (*plinkset).phenylphenyl.core.neighbor[i], coordarray[i], coordarray, box_dimension, my_nonbonded_params);
        }
		else if((newcoordarray[i].nodetype==2)||(newcoordarray[i].nodetype==3)){
			difference+=calc_charged_energy_cellstruct_otherpolymer((*plinkset).chargedcharged.core.number_neighbors[i], (*plinkset).chargedcharged.core.neighbor[i], newcoordarray[i], coordarray, box_dimension, my_nonbonded_params);
			difference-=calc_charged_energy_cellstruct_otherpolymer((*plinkset).chargedcharged.core.number_neighbors[i], (*plinkset).chargedcharged.core.neighbor[i], coordarray[i], coordarray, box_dimension, my_nonbonded_params);
		}
		difference+=calc_onebody_difference(coordarray[i], newcoordarray[i], chooseonebody(newcoordarray[i].nodetype, my_nonbonded_params), my_nonbonded_params.solvationparams);
	}
	acceptance_prob=exp(-difference/temperature);
	if((acceptance_prob>1.0)||(rand_double<acceptance_prob)){
        for(i=firstnode;i<=lastnode;i++){
            fmod_double_triple(&(newcoordarray[i].r), box_dimension);
            if(my_nonbonded_params.solvationparams.interface==1) updatecell_asymmetric(newcoordarray[i].r, &(*pmeshlink), i);
            updatecell(newcoordarray[i].r, &((*plinkset).shortrange), i);
            if(newcoordarray[i].nodetype==1) updatecell(newcoordarray[i].r, &((*plinkset).phenylphenyl), i);
            else if((newcoordarray[i].nodetype==2)||(newcoordarray[i].nodetype==3)) updatecell(newcoordarray[i].r, &((*plinkset).chargedcharged), i);
            coordarray[i]=newcoordarray[i];
        }
        for(i=firstnode;i<=lastnode;i++){
            if(my_nonbonded_params.solvationparams.interface==1) updateneighborlist_asymmetric(&(*pmeshlink), i, coordarray[i].r, meshpositions, box_dimension);
            updateneighborlist(&((*plinkset).shortrange), i, coordarray, box_dimension);
            if(newcoordarray[i].nodetype==1) updateneighborlist_pairenergies_phenyl(&((*plinkset).phenylphenyl), i, coordarray, box_dimension, my_nonbonded_params);
            else if((newcoordarray[i].nodetype==2)||(newcoordarray[i].nodetype==3)) updateneighborlist_pairenergies_charged(&((*plinkset).chargedcharged), i, coordarray, box_dimension, my_nonbonded_params);
        }
	}
}

void mc_translate_wholebilayer(int Nnodes, int Nsheets, int *chainlength, coord *newcoordarray, coord *coordarray, double max_translate, double_triple box_dimension, linkedlistset *plinkset, bonded_params my_bonded_params, nonbonded_params my_nonbonded_params, double temperature, monomernodes **monomerid, linkedlistasymmetric *pmeshlink, double_triple *meshpositions, int **mapmeshtoneighbors){
	int movingsheet=rand_int(Nsheets), firstnode, lastnode, i, foundfirst=0;
    for(i=0;i<Nnodes;i++){
        if(coordarray[i].leafid/2==movingsheet){
            if(foundfirst==0){
                firstnode=i;
                foundfirst=1;
            }
            lastnode=i;
        }
    }
	if(coordarray[lastnode].nodetype<0) lastnode--;                         //  accounting for possibility of long right terminus
	double acceptance_prob, difference=0;
	double_triple shift=scalar_multiply_double_triple(rand_unit_cube(), (rand_int(2)*2-1)*max_translate);
	for(i=firstnode;i<=lastnode;i++){
		newcoordarray[i]=coordarray[i];
		newcoordarray[i].r=add_double_triple(coordarray[i].r, shift);
		if(my_nonbonded_params.solvationparams.interface==1) newcoordarray[i].height=height_relative_to_interface(i, newcoordarray[i], my_nonbonded_params.solvationparams, (*pmeshlink), box_dimension, meshpositions, mapmeshtoneighbors);
		difference+=calc_hard_energy_cellstruct_otherbilayer((*plinkset).shortrange.core.number_neighbors[i], (*plinkset).shortrange.core.neighbor[i], newcoordarray[i], coordarray, box_dimension, my_nonbonded_params);
		difference-=calc_hard_energy_cellstruct_otherbilayer((*plinkset).shortrange.core.number_neighbors[i], (*plinkset).shortrange.core.neighbor[i], coordarray[i], coordarray, box_dimension, my_nonbonded_params);
        if(newcoordarray[i].nodetype==1){
			difference+=calc_phenyl_energy_cellstruct_otherbilayer((*plinkset).phenylphenyl.core.number_neighbors[i], (*plinkset).phenylphenyl.core.neighbor[i], newcoordarray[i], coordarray, box_dimension, my_nonbonded_params);
			difference-=calc_phenyl_energy_cellstruct_otherbilayer((*plinkset).phenylphenyl.core.number_neighbors[i], (*plinkset).phenylphenyl.core.neighbor[i], coordarray[i], coordarray, box_dimension, my_nonbonded_params);
        }
		else if((newcoordarray[i].nodetype==2)||(newcoordarray[i].nodetype==3)){
			difference+=calc_charged_energy_cellstruct_otherbilayer((*plinkset).chargedcharged.core.number_neighbors[i], (*plinkset).chargedcharged.core.neighbor[i], newcoordarray[i], coordarray, box_dimension, my_nonbonded_params);
			difference-=calc_charged_energy_cellstruct_otherbilayer((*plinkset).chargedcharged.core.number_neighbors[i], (*plinkset).chargedcharged.core.neighbor[i], coordarray[i], coordarray, box_dimension, my_nonbonded_params);
		}
		difference+=calc_onebody_difference(coordarray[i], newcoordarray[i], chooseonebody(newcoordarray[i].nodetype, my_nonbonded_params), my_nonbonded_params.solvationparams);
	}
	acceptance_prob=exp(-difference/temperature);
	if((acceptance_prob>1.0)||(rand_double<acceptance_prob)){
        for(i=firstnode;i<=lastnode;i++){
            fmod_double_triple(&(newcoordarray[i].r), box_dimension);
            if(my_nonbonded_params.solvationparams.interface==1) updatecell_asymmetric(newcoordarray[i].r, &(*pmeshlink), i);
            updatecell(newcoordarray[i].r, &((*plinkset).shortrange), i);
            if(newcoordarray[i].nodetype==1) updatecell(newcoordarray[i].r, &((*plinkset).phenylphenyl), i);
            else if((newcoordarray[i].nodetype==2)||(newcoordarray[i].nodetype==3)) updatecell(newcoordarray[i].r, &((*plinkset).chargedcharged), i);
            coordarray[i]=newcoordarray[i];
        }
        for(i=firstnode;i<=lastnode;i++){
            if(my_nonbonded_params.solvationparams.interface==1) updateneighborlist_asymmetric(&(*pmeshlink), i, coordarray[i].r, meshpositions, box_dimension);
            updateneighborlist(&((*plinkset).shortrange), i, coordarray, box_dimension);
            if(newcoordarray[i].nodetype==1) updateneighborlist_pairenergies_phenyl(&((*plinkset).phenylphenyl), i, coordarray, box_dimension, my_nonbonded_params);
            else if((newcoordarray[i].nodetype==2)||(newcoordarray[i].nodetype==3)) updateneighborlist_pairenergies_charged(&((*plinkset).chargedcharged), i, coordarray, box_dimension, my_nonbonded_params);
        }
	}
}

int mc_shiftbilayergap(int Nnodes, int Nsheets, int *chainlength, coord *newcoordarray, coord *coordarray, double max_translate, double_triple *pbox_dimension, linkedlistset *plinkset, bonded_params my_bonded_params, nonbonded_params my_nonbonded_params, double temperature, monomernodes **monomerid, double normalforceperunitarea, int sitespersheet){
	int changinginterval=rand_int(Nsheets), result=0, lowermover, uppermover, i, r;
    lowermover=changinginterval;
    uppermover=(changinginterval+1);                                            //  don't modulate this for purposes of shifting; only calculating energy
	double sep;
    int countlower=0, countupper=0;
	declare_array(int, count, Nsheets);
	declare_array(double, com, Nsheets);
	declare_array(double, first, Nsheets);
    
    for(r=0;r<Nsheets;r++){
        first[r]=coordarray[r*sitespersheet].r.z;
        for(i=0;i<sitespersheet;i++){
            if(coordarray[i].nodetype>=0){										//  accounting for possibility of long right terminus
                sep=coordarray[r*sitespersheet+i].r.z-first[r];
                recenter((sep), ((*pbox_dimension).z));
                com[r]+=sep;
                count[r]++;
            }
        }
        com[r]=com[r]/(1.*count[r])+first[r];
        fmod((com[r]), ((*pbox_dimension).z));
    }
    for(r=0;r<lowermover;r++){
        sep=com[lowermover]-com[r];
        fmod((sep), ((*pbox_dimension).z));
        com[r]=com[lowermover]-sep;                                         //  no boundary should separate sheets 0, ..., lowermover-1 from lowermover
    }
    for(r=uppermover;r<Nsheets;r++){
        sep=com[r]-com[lowermover];
        fmod((sep), ((*pbox_dimension).z));
        com[r]=com[lowermover]+sep;                                         //  no boundary should separate sheets uppermover, ..., Nsheets-1 from lowermover
    }
    for(r=0;r<Nsheets;r++){
        for(i=0;i<sitespersheet;i++){
            newcoordarray[r*sitespersheet+i]=coordarray[r*sitespersheet+i];
            sep=coordarray[r*sitespersheet+i].r.z-com[r];
            recenter(sep, (*pbox_dimension).z);
            newcoordarray[r*sitespersheet+i].r.z=com[r]+sep;
        }
    }                                                                       //  boundary is between sheet Nsheets-1 and 0
    
    double acceptance_prob, difference=0;
    double shift=rand_double*2-1;
    double_triple new_box_dimension=(*pbox_dimension);
    new_box_dimension.z+=shift;
    for(i=uppermover*sitespersheet;i<Nnodes;i++){
		newcoordarray[i].r.z+=shift;    //  move upper sheets and sheets above it by shift
    }
    
    uppermover=uppermover%Nsheets;                                          //  now shift uppermover from Nsheets to 0 if equal to Nsheets
    
    for(i=lowermover*sitespersheet;i<(lowermover+1)*sitespersheet;i++){
		difference+=calc_hard_energy_cellstruct_specificbilayer((*plinkset).shortrange.core.number_neighbors[i], (*plinkset).shortrange.core.neighbor[i], newcoordarray[i], newcoordarray, new_box_dimension, my_nonbonded_params, uppermover);
		difference-=calc_hard_energy_cellstruct_specificbilayer((*plinkset).shortrange.core.number_neighbors[i], (*plinkset).shortrange.core.neighbor[i], coordarray[i], coordarray, (*pbox_dimension), my_nonbonded_params, uppermover);
        if(newcoordarray[i].nodetype==1){
			difference+=calc_phenyl_energy_cellstruct_specificbilayer((*plinkset).phenylphenyl.core.number_neighbors[i], (*plinkset).phenylphenyl.core.neighbor[i], newcoordarray[i], newcoordarray, new_box_dimension, my_nonbonded_params, uppermover);
			difference-=calc_phenyl_energy_cellstruct_specificbilayer((*plinkset).phenylphenyl.core.number_neighbors[i], (*plinkset).phenylphenyl.core.neighbor[i], coordarray[i], coordarray, (*pbox_dimension), my_nonbonded_params, uppermover);
        }
		else if((newcoordarray[i].nodetype==2)||(newcoordarray[i].nodetype==3)){
			difference+=calc_charged_energy_cellstruct_specificbilayer((*plinkset).chargedcharged.core.number_neighbors[i], (*plinkset).chargedcharged.core.neighbor[i], newcoordarray[i], newcoordarray, new_box_dimension, my_nonbonded_params, uppermover);
			difference-=calc_charged_energy_cellstruct_specificbilayer((*plinkset).chargedcharged.core.number_neighbors[i], (*plinkset).chargedcharged.core.neighbor[i], coordarray[i], coordarray, (*pbox_dimension), my_nonbonded_params, uppermover);
		}        
        if(my_nonbonded_params.solvationparams.interface!=0) my_exit("mc_translate_wholebilyaer_changez not written for interface!\n");
	}
	acceptance_prob=exp(-(difference+normalforceperunitarea*(new_box_dimension.z-(*pbox_dimension).z)*(*pbox_dimension).x*(*pbox_dimension).y)/temperature);
	if((acceptance_prob>1.0)||(rand_double<acceptance_prob)){
        result=1;
		printf("box from %f to ", (*pbox_dimension).z);
        (*pbox_dimension)=new_box_dimension;
		printf("%f\n", (*pbox_dimension).z);
		for(i=0;i<Nnodes;i++){
            fmod_double_triple(&(newcoordarray[i].r), (*pbox_dimension));
            coordarray[i]=newcoordarray[i];
        }
	}
	free(count);
	free(com);
	free(first);
    return result;
}

void calc_cellstruct_leafid(coord *coordarray, bonded_params my_bonded_params, nonbonded_params my_nonbonded_params, int Nchains, int *chainlength, monomernodes **monomerid, double_triple box_dimension, energycomponents *pmy_energycomponents, linkedlistset linkset, double_triple *pavboxdimension){
	int i, j, left=-1, twoleft=-1, threeleft=-1, fourleft=-1, backboneid, sidechainid;
	(*pavboxdimension).x+=box_dimension.x;
	(*pavboxdimension).y+=box_dimension.y;
	(*pavboxdimension).z+=box_dimension.z;
	for(i=0;i<Nchains;i++){
		if(chainlength[i]>2){
			for(j=0;j<chainlength[i];j++){				
                backboneid=monomerid[i][j].backbone;
				if(my_nonbonded_params.solvationparams.interface==1) calc_onebody_energy_components(coordarray[backboneid], chooseonebody(coordarray[backboneid].nodetype, my_nonbonded_params), my_nonbonded_params.solvationparams, pmy_energycomponents);
				
				//	using old sidechainid:
				//	backbone only has short-ranged interactions (zero)
                
				sidechainid=monomerid[i][j].sidechain;
                if(coordarray[sidechainid].nodetype>=0){
                    if(my_nonbonded_params.solvationparams.interface==1) calc_onebody_energy_components(coordarray[sidechainid], chooseonebody(coordarray[sidechainid].nodetype, my_nonbonded_params), my_nonbonded_params.solvationparams, pmy_energycomponents);
                    if(coordarray[sidechainid].nodetype==1){
						calc_nonbonded_energy_components_cellstruct_leafid(sidechainid, backboneid, -1, -1, linkset.phenylphenyl.core, coordarray, box_dimension, my_nonbonded_params, pmy_energycomponents);
                    }
					if(coordarray[sidechainid].nodetype>1){
						calc_nonbonded_energy_components_cellstruct_leafid(sidechainid, backboneid, -1, -1, linkset.chargedcharged.core, coordarray, box_dimension, my_nonbonded_params, pmy_energycomponents);
					}
                    if(coordarray[sidechainid].nodetype>=0) (*pmy_energycomponents).sidechain+=calc_sidechain_bonded_energy(my_bonded_params.sidechain[coordarray[sidechainid].nodetype], coordarray[backboneid], coordarray[sidechainid], box_dimension);
                }
				if(j==2){
                    (*pmy_energycomponents).backbone+=calc_backbone_bonded_energy(my_bonded_params, coordarray[left], coordarray[twoleft], coordarray[backboneid], box_dimension);
                }
				else if(j>=3){
                    (*pmy_energycomponents).backbone+=calc_backbone_bonded_energy_norim1i(my_bonded_params, coordarray[left], coordarray[twoleft], coordarray[backboneid], box_dimension);
                }
				fourleft=threeleft;
				threeleft=twoleft;
				twoleft=left;
				left=backboneid;
			}
            
			//	backbone only has short-ranged interactions (zero)
            
		}
		else if(chainlength[i]==2){
			for(j=0;j<chainlength[i];j++){
				backboneid=monomerid[i][j].backbone;
				if(my_nonbonded_params.solvationparams.interface==1) calc_onebody_energy_components(coordarray[backboneid], chooseonebody(coordarray[backboneid].nodetype, my_nonbonded_params), my_nonbonded_params.solvationparams, pmy_energycomponents);
				sidechainid=monomerid[i][j].sidechain;
                if(coordarray[sidechainid].nodetype>=0){
                    if(my_nonbonded_params.solvationparams.interface==1) calc_onebody_energy_components(coordarray[sidechainid], chooseonebody(coordarray[sidechainid].nodetype, my_nonbonded_params), my_nonbonded_params.solvationparams, pmy_energycomponents);
                    
                    if(coordarray[sidechainid].nodetype==1) calc_nonbonded_energy_components_cellstruct_leafid(sidechainid, backboneid, -1, -1, linkset.phenylphenyl.core, coordarray, box_dimension, my_nonbonded_params, pmy_energycomponents);
                    if(coordarray[sidechainid].nodetype>1) calc_nonbonded_energy_components_cellstruct_leafid(sidechainid, backboneid, -1, -1, linkset.chargedcharged.core, coordarray, box_dimension, my_nonbonded_params, pmy_energycomponents);
                    
                    if(coordarray[sidechainid].nodetype>=0) (*pmy_energycomponents).sidechain+=calc_sidechain_bonded_energy(my_bonded_params.sidechain[coordarray[sidechainid].nodetype], coordarray[backboneid], coordarray[sidechainid], box_dimension);
                }
				if(j==1){
					
					//	backbone only has short-ranged interactions (zero)
                    
					(*pmy_energycomponents).backbone+=calc_backbone_bonded_energy_onlyleft(my_bonded_params, coordarray[backboneid], coordarray[left], box_dimension);
				}
				left=backboneid;
			}
			
			//	backbone only has short-ranged interactions (zero)
            
		}
		else{
			for(j=0;j<chainlength[i];j++){
				backboneid=monomerid[i][j].backbone;
				if(my_nonbonded_params.solvationparams.interface==1) calc_onebody_energy_components(coordarray[backboneid], chooseonebody(coordarray[backboneid].nodetype, my_nonbonded_params), my_nonbonded_params.solvationparams, pmy_energycomponents);
				sidechainid=monomerid[i][j].sidechain;
                if(coordarray[sidechainid].nodetype>=0){
                    if(my_nonbonded_params.solvationparams.interface==1) calc_onebody_energy_components(coordarray[sidechainid], chooseonebody(coordarray[sidechainid].nodetype, my_nonbonded_params), my_nonbonded_params.solvationparams, pmy_energycomponents);
                    
                    if(coordarray[sidechainid].nodetype==1) calc_nonbonded_energy_components_cellstruct_leafid(sidechainid, backboneid, -1, -1, linkset.phenylphenyl.core, coordarray, box_dimension, my_nonbonded_params, pmy_energycomponents);
                    if(coordarray[sidechainid].nodetype>1) calc_nonbonded_energy_components_cellstruct_leafid(sidechainid, backboneid, -1, -1, linkset.chargedcharged.core, coordarray, box_dimension, my_nonbonded_params, pmy_energycomponents);
                    
                    if(coordarray[sidechainid].nodetype>=0) (*pmy_energycomponents).sidechain+=calc_sidechain_bonded_energy(my_bonded_params.sidechain[coordarray[sidechainid].nodetype], coordarray[backboneid], coordarray[sidechainid], box_dimension);
                }
			}
        }
	}
}

void calc_nonbonded_energy_pair_components_leafid_cluster(coord a, coord b, double_triple box_dimension, nonbonded_params params, energycomponents *pmy_energycomponents, int *cluster_index_of_chain, int *size_of_cluster, int *head_of_cluster, int *list_of_chain, int *pNclusters){
    double val;
    int chaina, chainb, clustera, clusterb, mychain, nextchain;
    if(a.leafid==b.leafid){
		if(a.chainid==b.chainid){
			if(a.nodetype==0){
				if(b.nodetype==0){
					(*pmy_energycomponents).nnsamepoly+=hard_energy(a.r, b.r, params.rhard0, params.rhard0, box_dimension);
					return;
				}
				else if(b.nodetype==1){
					(*pmy_energycomponents).nnsamepoly+=hard_energy(a.r, b.r, params.rhard0, params.rhard1, box_dimension);
					return;
				}
				else if(b.nodetype==2){
					(*pmy_energycomponents).cnsamepoly+=hard_energy(a.r, b.r, params.rhard0, params.p2.rhard, box_dimension);
					return;
				}
				else if(b.nodetype==3){
					(*pmy_energycomponents).cnsamepoly+=hard_energy(a.r, b.r, params.rhard0, params.p3.rhard, box_dimension);
					return;
				}
				else my_exit("don't have energy functions for these node types!");
			}
			else if(a.nodetype==1){
				if(b.nodetype==0){
					(*pmy_energycomponents).nnsamepoly+=hard_energy(a.r, b.r, params.rhard1, params.rhard0, box_dimension);
					return;
				}
				else if(b.nodetype==1){
					(*pmy_energycomponents).nnsamepoly+=phenylphenyl_energy(a, b, params.p11vac, box_dimension, params.phenylcutoff2, params.solvationparams);
					return;
				}
				else if(b.nodetype==2){
					(*pmy_energycomponents).cnsamepoly+=hard_energy(a.r, b.r, params.rhard1, params.p2.rhard, box_dimension);
					return;
				}
				else if(b.nodetype==3){
					(*pmy_energycomponents).cnsamepoly+=hard_energy(a.r, b.r, params.rhard1, params.p3.rhard, box_dimension);
					return;
				}
				else my_exit("don't have energy functions for these node types!");
			}
			else if(a.nodetype==2){
				if(b.nodetype==0){
					(*pmy_energycomponents).cnsamepoly+=hard_energy(a.r, b.r, params.p2.rhard, params.rhard0, box_dimension);
					return;
				}
				else if(b.nodetype==1){
					(*pmy_energycomponents).cnsamepoly+=hard_energy(a.r, b.r, params.p2.rhard, params.rhard1, box_dimension);
					return;
				}
				else if(b.nodetype==2){
					(*pmy_energycomponents).cclikesamepoly+=electrostatic_energy(a, b, params.solvationparams, params.p2, params.p2, box_dimension, params.cutoff2);
					return;
				}
				else if(b.nodetype==3){
					(*pmy_energycomponents).ccunlikesamepoly+=electrostatic_energy(a, b, params.solvationparams, params.p2, params.p3, box_dimension, params.cutoff2);
					return;
				}
				else my_exit("don't have energy functions for these node types!");
			}
			else if(a.nodetype==3){
				if(b.nodetype==0){
					(*pmy_energycomponents).cnsamepoly+=hard_energy(a.r, b.r, params.p3.rhard, params.rhard0, box_dimension);
					return;
				}
				else if(b.nodetype==1){
					(*pmy_energycomponents).cnsamepoly+=hard_energy(a.r, b.r, params.p3.rhard, params.rhard1, box_dimension);
					return;
				}
				else if(b.nodetype==2){
					(*pmy_energycomponents).ccunlikesamepoly+=electrostatic_energy(a, b, params.solvationparams, params.p3, params.p2, box_dimension, params.cutoff2);
					return;
				}
				else if(b.nodetype==3){
					(*pmy_energycomponents).cclikesamepoly+=electrostatic_energy(a, b, params.solvationparams, params.p3, params.p3, box_dimension, params.cutoff2);
					return;
				}
				else my_exit("don't have energy functions for these node types!");
			}
			else my_exit("don't have energy functions for these node types!");
		}
		else{
			if(a.nodetype==0){
				if(b.nodetype==0){
					(*pmy_energycomponents).nn+=hard_energy(a.r, b.r, params.rhard0, params.rhard0, box_dimension);
					return;
				}
				else if(b.nodetype==1){
					(*pmy_energycomponents).nn+=hard_energy(a.r, b.r, params.rhard0, params.rhard1, box_dimension);
					return;
				}
				else if(b.nodetype==2){
					(*pmy_energycomponents).cn+=hard_energy(a.r, b.r, params.rhard0, params.p2.rhard, box_dimension);
					return;
				}
				else if(b.nodetype==3){
					(*pmy_energycomponents).cn+=hard_energy(a.r, b.r, params.rhard0, params.p3.rhard, box_dimension);
					return;
				}
				else my_exit("don't have energy functions for these node types!");
			}
			else if(a.nodetype==1){
				if(b.nodetype==0){
					(*pmy_energycomponents).nn+=hard_energy(a.r, b.r, params.rhard1, params.rhard0, box_dimension);
					return;
				}
				else if(b.nodetype==1){
                    val=phenylphenyl_energy(a, b, params.p11vac, box_dimension, params.phenylcutoff2, params.solvationparams);
                    
                    if(val<0){
                        chaina=a.chainid;
                        chainb=b.chainid;
                        if(chaina!=chainb){
                            clustera=cluster_index_of_chain[chaina];
                            clusterb=cluster_index_of_chain[chainb];
                            if(clustera==-1){
                                if(clusterb==-1){
                                    cluster_index_of_chain[chaina]=*pNclusters;
                                    cluster_index_of_chain[chainb]=*pNclusters;
                                    size_of_cluster[*pNclusters]=2;
                                    head_of_cluster[*pNclusters]=chaina;
                                    list_of_chain[chaina]=chainb;
                                    (*pNclusters)++;
                                }
                                else{
                                    cluster_index_of_chain[chaina]=clusterb;
                                    size_of_cluster[clusterb]++;
                                    list_of_chain[chaina]=head_of_cluster[clusterb];
                                    head_of_cluster[clusterb]=chaina;
                                }
                            }
                            else{
                                if(clusterb==-1){
                                    cluster_index_of_chain[chainb]=clustera;
                                    size_of_cluster[clustera]++;
                                    list_of_chain[chainb]=head_of_cluster[clustera];
                                    head_of_cluster[clustera]=chainb;
                                }
                                else if(clustera!=clusterb){
                                    size_of_cluster[clustera]+=size_of_cluster[clusterb];
                                    size_of_cluster[clusterb]=0;        //  will still be in list of clusters (don't want to have to shift list)
                                    mychain=head_of_cluster[clusterb];
                                    cluster_index_of_chain[mychain]=clustera;
                                    nextchain=list_of_chain[mychain];
                                    while(nextchain!=-1){
                                        mychain=nextchain;
                                        cluster_index_of_chain[nextchain]=clustera;
                                        nextchain=list_of_chain[mychain];
                                    }
                                    list_of_chain[mychain]=head_of_cluster[clustera];
                                    head_of_cluster[clustera]=head_of_cluster[clusterb];
                                }
                            }
                        }   
                    }   
                    
					(*pmy_energycomponents).nn+=val;
					return;
				}
				else if(b.nodetype==2){
					(*pmy_energycomponents).cn+=hard_energy(a.r, b.r, params.rhard1, params.p2.rhard, box_dimension);
					return;
				}
				else if(b.nodetype==3){
					(*pmy_energycomponents).cn+=hard_energy(a.r, b.r, params.rhard1, params.p3.rhard, box_dimension);
					return;
				}
				else my_exit("don't have energy functions for these node types!");
			}
			else if(a.nodetype==2){
				if(b.nodetype==0){
					(*pmy_energycomponents).cn+=hard_energy(a.r, b.r, params.p2.rhard, params.rhard0, box_dimension);
					return;
				}
				else if(b.nodetype==1){
					(*pmy_energycomponents).cn+=hard_energy(a.r, b.r, params.p2.rhard, params.rhard1, box_dimension);
					return;
				}
				else if(b.nodetype==2){
					(*pmy_energycomponents).cclike+=electrostatic_energy(a, b, params.solvationparams, params.p2, params.p2, box_dimension, params.cutoff2);
					return;
				}
				else if(b.nodetype==3){
                    val=electrostatic_energy(a, b, params.solvationparams, params.p2, params.p3, box_dimension, params.cutoff2);
                    
                    if(val<0){
                        chaina=a.chainid;
                        chainb=b.chainid;
                        if(chaina!=chainb){
                            clustera=cluster_index_of_chain[chaina];
                            clusterb=cluster_index_of_chain[chainb];
                            if(clustera==-1){
                                if(clusterb==-1){
                                    cluster_index_of_chain[chaina]=*pNclusters;
                                    cluster_index_of_chain[chainb]=*pNclusters;
                                    size_of_cluster[*pNclusters]=2;
                                    head_of_cluster[*pNclusters]=chaina;
                                    list_of_chain[chaina]=chainb;
                                    (*pNclusters)++;
                                }
                                else{
                                    cluster_index_of_chain[chaina]=clusterb;
                                    size_of_cluster[clusterb]++;
                                    list_of_chain[chaina]=head_of_cluster[clusterb];
                                    head_of_cluster[clusterb]=chaina;
                                }
                            }
                            else{
                                if(clusterb==-1){
                                    cluster_index_of_chain[chainb]=clustera;
                                    size_of_cluster[clustera]++;
                                    list_of_chain[chainb]=head_of_cluster[clustera];
                                    head_of_cluster[clustera]=chainb;
                                }
                                else if(clustera!=clusterb){
                                    size_of_cluster[clustera]+=size_of_cluster[clusterb];
                                    size_of_cluster[clusterb]=0;        //  will still be in list of clusters (don't want to have to shift list)
                                    mychain=head_of_cluster[clusterb];
                                    cluster_index_of_chain[mychain]=clustera;
                                    nextchain=list_of_chain[mychain];
                                    while(nextchain!=-1){
                                        mychain=nextchain;
                                        cluster_index_of_chain[nextchain]=clustera;
                                        nextchain=list_of_chain[mychain];
                                    }
                                    list_of_chain[mychain]=head_of_cluster[clustera];
                                    head_of_cluster[clustera]=head_of_cluster[clusterb];
                                }
                            }
                        }   
                    }   
                    
					(*pmy_energycomponents).ccunlike+=val;
					return;
				}
				else my_exit("don't have energy functions for these node types!");
			}
			else if(a.nodetype==3){
				if(b.nodetype==0){
					(*pmy_energycomponents).cn+=hard_energy(a.r, b.r, params.p3.rhard, params.rhard0, box_dimension);
					return;
				}
				else if(b.nodetype==1){
					(*pmy_energycomponents).cn+=hard_energy(a.r, b.r, params.p3.rhard, params.rhard1, box_dimension);
					return;
				}
				else if(b.nodetype==2){
                    val=electrostatic_energy(a, b, params.solvationparams, params.p3, params.p2, box_dimension, params.cutoff2);
                    
                    if(val<0){
                        chaina=a.chainid;
                        chainb=b.chainid;
                        if(chaina!=chainb){
                            clustera=cluster_index_of_chain[chaina];
                            clusterb=cluster_index_of_chain[chainb];
                            if(clustera==-1){
                                if(clusterb==-1){
                                    cluster_index_of_chain[chaina]=*pNclusters;
                                    cluster_index_of_chain[chainb]=*pNclusters;
                                    size_of_cluster[*pNclusters]=2;
                                    head_of_cluster[*pNclusters]=chaina;
                                    list_of_chain[chaina]=chainb;
                                    (*pNclusters)++;
                                }
                                else{
                                    cluster_index_of_chain[chaina]=clusterb;
                                    size_of_cluster[clusterb]++;
                                    list_of_chain[chaina]=head_of_cluster[clusterb];
                                    head_of_cluster[clusterb]=chaina;
                                }
                            }
                            else{
                                if(clusterb==-1){
                                    cluster_index_of_chain[chainb]=clustera;
                                    size_of_cluster[clustera]++;
                                    list_of_chain[chainb]=head_of_cluster[clustera];
                                    head_of_cluster[clustera]=chainb;
                                }
                                else if(clustera!=clusterb){
                                    size_of_cluster[clustera]+=size_of_cluster[clusterb];
                                    size_of_cluster[clusterb]=0;        //  will still be in list of clusters (don't want to have to shift list)
                                    mychain=head_of_cluster[clusterb];
                                    cluster_index_of_chain[mychain]=clustera;
                                    nextchain=list_of_chain[mychain];
                                    while(nextchain!=-1){
                                        mychain=nextchain;
                                        cluster_index_of_chain[nextchain]=clustera;
                                        nextchain=list_of_chain[mychain];
                                    }
                                    list_of_chain[mychain]=head_of_cluster[clustera];
                                    head_of_cluster[clustera]=head_of_cluster[clusterb];
                                }
                            }
                        }   
                    }   
                    
					(*pmy_energycomponents).ccunlike+=val;
					return;
				}
				else if(b.nodetype==3){
					(*pmy_energycomponents).cclike+=electrostatic_energy(a, b, params.solvationparams, params.p3, params.p3, box_dimension, params.cutoff2);
					return;
				}
				else my_exit("don't have energy functions for these node types!");
			}
			else my_exit("don't have energy functions for these node types!");
		}
    }
    else if((a.leafid)/2==(b.leafid)/2){
        if(a.nodetype==0){
            if(b.nodetype==0){
                (*pmy_energycomponents).nncross+=hard_energy(a.r, b.r, params.rhard0, params.rhard0, box_dimension);
                return;
            }
            else if(b.nodetype==1){
                (*pmy_energycomponents).nncross+=hard_energy(a.r, b.r, params.rhard0, params.rhard1, box_dimension);
                return;
            }
            else if(b.nodetype==2){
                (*pmy_energycomponents).cncross+=hard_energy(a.r, b.r, params.rhard0, params.p2.rhard, box_dimension);
                return;
            }
            else if(b.nodetype==3){
                (*pmy_energycomponents).cncross+=hard_energy(a.r, b.r, params.rhard0, params.p3.rhard, box_dimension);
                return;
            }
            else my_exit("don't have energy functions for these node types!");
        }
        else if(a.nodetype==1){
            if(b.nodetype==0){
                (*pmy_energycomponents).nncross+=hard_energy(a.r, b.r, params.rhard1, params.rhard0, box_dimension);
                return;
            }
            else if(b.nodetype==1){
				(*pmy_energycomponents).nncross+=phenylphenyl_energy(a, b, params.p11vac, box_dimension, params.phenylcutoff2, params.solvationparams);                
				return;
            }
            else if(b.nodetype==2){
                (*pmy_energycomponents).cncross+=hard_energy(a.r, b.r, params.rhard1, params.p2.rhard, box_dimension);
                return;
            }
            else if(b.nodetype==3){
                (*pmy_energycomponents).cncross+=hard_energy(a.r, b.r, params.rhard1, params.p3.rhard, box_dimension);
                return;
            }
            else my_exit("don't have energy functions for these node types!");
        }
        else if(a.nodetype==2){
            if(b.nodetype==0){
                (*pmy_energycomponents).cncross+=hard_energy(a.r, b.r, params.p2.rhard, params.rhard0, box_dimension);
                return;
            }
            else if(b.nodetype==1){
                (*pmy_energycomponents).cncross+=hard_energy(a.r, b.r, params.p2.rhard, params.rhard1, box_dimension);
                return;
            }
            else if(b.nodetype==2){
                (*pmy_energycomponents).cclikecross+=electrostatic_energy(a, b, params.solvationparams, params.p2, params.p2, box_dimension, params.cutoff2);
                return;
            }
            else if(b.nodetype==3){
                (*pmy_energycomponents).ccunlikecross+=electrostatic_energy(a, b, params.solvationparams, params.p2, params.p3, box_dimension, params.cutoff2);
                return;
            }
            else my_exit("don't have energy functions for these node types!");
        }
        else if(a.nodetype==3){
            if(b.nodetype==0){
                (*pmy_energycomponents).cncross+=hard_energy(a.r, b.r, params.p3.rhard, params.rhard0, box_dimension);
                return;
            }
            else if(b.nodetype==1){
                (*pmy_energycomponents).cncross+=hard_energy(a.r, b.r, params.p3.rhard, params.rhard1, box_dimension);
                return;
            }
            else if(b.nodetype==2){
                (*pmy_energycomponents).ccunlikecross+=electrostatic_energy(a, b, params.solvationparams, params.p3, params.p2, box_dimension, params.cutoff2);
                return;
            }
            else if(b.nodetype==3){
                (*pmy_energycomponents).cclikecross+=electrostatic_energy(a, b, params.solvationparams, params.p3, params.p3, box_dimension, params.cutoff2);
                return;
            }
            else my_exit("don't have energy functions for these node types!");
        }
        else my_exit("don't have energy functions for these node types!");
    }
    else{
        if(a.nodetype==0){
            if(b.nodetype==0){
                (*pmy_energycomponents).nndifferent+=hard_energy(a.r, b.r, params.rhard0, params.rhard0, box_dimension);
                return;
            }
            else if(b.nodetype==1){
                (*pmy_energycomponents).nndifferent+=hard_energy(a.r, b.r, params.rhard0, params.rhard1, box_dimension);
                return;
            }
            else if(b.nodetype==2){
                (*pmy_energycomponents).cndifferent+=hard_energy(a.r, b.r, params.rhard0, params.p2.rhard, box_dimension);
                return;
            }
            else if(b.nodetype==3){
                (*pmy_energycomponents).cndifferent+=hard_energy(a.r, b.r, params.rhard0, params.p3.rhard, box_dimension);
                return;
            }
            else my_exit("don't have energy functions for these node types!");
        }
        else if(a.nodetype==1){
            if(b.nodetype==0){
                (*pmy_energycomponents).nndifferent+=hard_energy(a.r, b.r, params.rhard1, params.rhard0, box_dimension);
                return;
            }
            else if(b.nodetype==1){
				(*pmy_energycomponents).nndifferent+=phenylphenyl_energy(a, b, params.p11vac, box_dimension, params.phenylcutoff2, params.solvationparams);
				return;
            }
            else if(b.nodetype==2){
                (*pmy_energycomponents).cndifferent+=hard_energy(a.r, b.r, params.rhard1, params.p2.rhard, box_dimension);
                return;
            }
            else if(b.nodetype==3){
                (*pmy_energycomponents).cndifferent+=hard_energy(a.r, b.r, params.rhard1, params.p3.rhard, box_dimension);
                return;
            }
            else my_exit("don't have energy functions for these node types!");
        }
        else if(a.nodetype==2){
            if(b.nodetype==0){
                (*pmy_energycomponents).cndifferent+=hard_energy(a.r, b.r, params.p2.rhard, params.rhard0, box_dimension);
                return;
            }
            else if(b.nodetype==1){
                (*pmy_energycomponents).cndifferent+=hard_energy(a.r, b.r, params.p2.rhard, params.rhard1, box_dimension);
                return;
            }
            else if(b.nodetype==2){
                (*pmy_energycomponents).cclikedifferent+=electrostatic_energy(a, b, params.solvationparams, params.p2, params.p2, box_dimension, params.cutoff2);
                return;
            }
            else if(b.nodetype==3){
                (*pmy_energycomponents).ccunlikedifferent+=electrostatic_energy(a, b, params.solvationparams, params.p2, params.p3, box_dimension, params.cutoff2);
                return;
            }
            else my_exit("don't have energy functions for these node types!");
        }
        else if(a.nodetype==3){
            if(b.nodetype==0){
                (*pmy_energycomponents).cndifferent+=hard_energy(a.r, b.r, params.p3.rhard, params.rhard0, box_dimension);
                return;
            }
            else if(b.nodetype==1){
                (*pmy_energycomponents).cndifferent+=hard_energy(a.r, b.r, params.p3.rhard, params.rhard1, box_dimension);
                return;
            }
            else if(b.nodetype==2){
                (*pmy_energycomponents).ccunlikedifferent+=electrostatic_energy(a, b, params.solvationparams, params.p3, params.p2, box_dimension, params.cutoff2);
                return;
            }
            else if(b.nodetype==3){
                (*pmy_energycomponents).cclikedifferent+=electrostatic_energy(a, b, params.solvationparams, params.p3, params.p3, box_dimension, params.cutoff2);
                return;
            }
            else my_exit("don't have energy functions for these node types!");
        }
        else my_exit("don't have energy functions for these node types!");
    }
}

void calc_nonbonded_energy_components_cellstruct_leafid_cluster(int n, int neighbor1, int neighbor2, int neighbor3, linkedlist link, coord *coordarray, double_triple box_dimension, nonbonded_params my_nonbonded_params, energycomponents *pmy_energycomponents, int *cluster_index_of_chain, int *size_of_cluster, int *head_of_cluster, int *list_of_chain, int *pNclusters){
	int i, j, k, i_image, j_image, k_image, index;
	coord coord_image=coordarray[n];
	for(i=-1;i<=1;i++){																	//	loop through neighboring cells
		if((link.cell[n].x==0)&&(i==-1)){
			i_image=link.cellsperside.x-1;
			coord_image.r.x=coordarray[n].r.x+box_dimension.x;
		}
		else if((link.cell[n].x==link.cellsperside.x-1)&&(i==1)){
			i_image=0;
			coord_image.r.x=coordarray[n].r.x-box_dimension.x;
		}
		else{
            coord_image.r.x=coordarray[n].r.x;
            i_image=link.cell[n].x+i;
        }
		for(j=-1;j<=1;j++){
			if((link.cell[n].y==0)&&(j==-1)){
				j_image=link.cellsperside.y-1;
				coord_image.r.y=coordarray[n].r.y+box_dimension.y;
			}
			else if((link.cell[n].y==link.cellsperside.y-1)&&(j==1)){
				j_image=0;
				coord_image.r.y=coordarray[n].r.y-box_dimension.y;
			}
			else{
                coord_image.r.y=coordarray[n].r.y;
                j_image=link.cell[n].y+j;
            }
			for(k=-1;k<=1;k++){
				if((link.cell[n].z==0)&&(k==-1)){
					k_image=link.cellsperside.z-1;
					coord_image.r.z=coordarray[n].r.z+box_dimension.z;
				}
				else if((link.cell[n].z==link.cellsperside.z-1)&&(k==1)){
					k_image=0;
					coord_image.r.z=coordarray[n].r.z-box_dimension.z;
				}
				else{
                    coord_image.r.z=coordarray[n].r.z;
                    k_image=link.cell[n].z+k;
                }
				index=link.head[i_image][j_image][k_image];
				while(index>=0){
					if((index!=n)&&(index!=neighbor1)&&(index!=neighbor2)&&(index!=neighbor3)){					//	don't calculate nonbonded energy for bonded neighbors
						calc_nonbonded_energy_pair_components_leafid_cluster(coord_image, coordarray[index], box_dimension, my_nonbonded_params, pmy_energycomponents, cluster_index_of_chain, size_of_cluster, head_of_cluster, list_of_chain, pNclusters);
					}
					index=link.list[index];
				}
			}
		}
	}
}

void calc_cellstruct_leafid_cluster(coord *coordarray, bonded_params my_bonded_params, nonbonded_params my_nonbonded_params, int Nchains, int *chainlength, monomernodes **monomerid, double_triple box_dimension, energycomponents *pmy_energycomponents, linkedlistset linkset, double_triple *pavboxdimension, int *clusterdistribution){
	int i, j, left=-1, twoleft=-1, threeleft=-1, fourleft=-1, backboneid, sidechainid;
	(*pavboxdimension).x+=box_dimension.x;
	(*pavboxdimension).y+=box_dimension.y;
	(*pavboxdimension).z+=box_dimension.z;
    
    declare_array(int, cluster_index_of_chain, Nchains);
    declare_array(int, size_of_cluster, Nchains);
    declare_array(int, head_of_cluster, Nchains);
    declare_array(int, list_of_chain, Nchains);
    for(i=0;i<Nchains;i++){
        cluster_index_of_chain[i]=-1;
        list_of_chain[i]=-1;
    }
    int Nclusters=0;
    
	for(i=0;i<Nchains;i++){
		if(chainlength[i]>2){
			for(j=0;j<chainlength[i];j++){				
                backboneid=monomerid[i][j].backbone;
				if(my_nonbonded_params.solvationparams.interface==1) calc_onebody_energy_components(coordarray[backboneid], chooseonebody(coordarray[backboneid].nodetype, my_nonbonded_params), my_nonbonded_params.solvationparams, pmy_energycomponents);
				
				//	using old sidechainid:
				//	backbone only has short-ranged interactions (zero)
                
				sidechainid=monomerid[i][j].sidechain;
                if(coordarray[sidechainid].nodetype>=0){
                    if(my_nonbonded_params.solvationparams.interface==1) calc_onebody_energy_components(coordarray[sidechainid], chooseonebody(coordarray[sidechainid].nodetype, my_nonbonded_params), my_nonbonded_params.solvationparams, pmy_energycomponents);
                    if(coordarray[sidechainid].nodetype==1){
						calc_nonbonded_energy_components_cellstruct_leafid_cluster(sidechainid, backboneid, -1, -1, linkset.phenylphenyl.core, coordarray, box_dimension, my_nonbonded_params, pmy_energycomponents, cluster_index_of_chain, size_of_cluster, head_of_cluster, list_of_chain, &Nclusters);
                    }
					if(coordarray[sidechainid].nodetype>1){
						calc_nonbonded_energy_components_cellstruct_leafid_cluster(sidechainid, backboneid, -1, -1, linkset.chargedcharged.core, coordarray, box_dimension, my_nonbonded_params, pmy_energycomponents, cluster_index_of_chain, size_of_cluster, head_of_cluster, list_of_chain, &Nclusters);
					}
                    if(coordarray[sidechainid].nodetype>=0) (*pmy_energycomponents).sidechain+=calc_sidechain_bonded_energy(my_bonded_params.sidechain[coordarray[sidechainid].nodetype], coordarray[backboneid], coordarray[sidechainid], box_dimension);
                }
				if(j==2){
                    (*pmy_energycomponents).backbone+=calc_backbone_bonded_energy(my_bonded_params, coordarray[left], coordarray[twoleft], coordarray[backboneid], box_dimension);
                }
				else if(j>=3){
                    (*pmy_energycomponents).backbone+=calc_backbone_bonded_energy_norim1i(my_bonded_params, coordarray[left], coordarray[twoleft], coordarray[backboneid], box_dimension);
                }
				fourleft=threeleft;
				threeleft=twoleft;
				twoleft=left;
				left=backboneid;
			}
            
			//	backbone only has short-ranged interactions (zero)
            
		}
		else if(chainlength[i]==2){
			for(j=0;j<chainlength[i];j++){
				backboneid=monomerid[i][j].backbone;
				if(my_nonbonded_params.solvationparams.interface==1) calc_onebody_energy_components(coordarray[backboneid], chooseonebody(coordarray[backboneid].nodetype, my_nonbonded_params), my_nonbonded_params.solvationparams, pmy_energycomponents);
				sidechainid=monomerid[i][j].sidechain;
                if(coordarray[sidechainid].nodetype>=0){
                    if(my_nonbonded_params.solvationparams.interface==1) calc_onebody_energy_components(coordarray[sidechainid], chooseonebody(coordarray[sidechainid].nodetype, my_nonbonded_params), my_nonbonded_params.solvationparams, pmy_energycomponents);
                    
                    if(coordarray[sidechainid].nodetype==1) calc_nonbonded_energy_components_cellstruct_leafid_cluster(sidechainid, backboneid, -1, -1, linkset.phenylphenyl.core, coordarray, box_dimension, my_nonbonded_params, pmy_energycomponents, cluster_index_of_chain, size_of_cluster, head_of_cluster, list_of_chain, &Nclusters);
                    if(coordarray[sidechainid].nodetype>1) calc_nonbonded_energy_components_cellstruct_leafid_cluster(sidechainid, backboneid, -1, -1, linkset.chargedcharged.core, coordarray, box_dimension, my_nonbonded_params, pmy_energycomponents, cluster_index_of_chain, size_of_cluster, head_of_cluster, list_of_chain, &Nclusters);
                    
                    if(coordarray[sidechainid].nodetype>=0) (*pmy_energycomponents).sidechain+=calc_sidechain_bonded_energy(my_bonded_params.sidechain[coordarray[sidechainid].nodetype], coordarray[backboneid], coordarray[sidechainid], box_dimension);
                }
				if(j==1){
					
					//	backbone only has short-ranged interactions (zero)
                    
					(*pmy_energycomponents).backbone+=calc_backbone_bonded_energy_onlyleft(my_bonded_params, coordarray[backboneid], coordarray[left], box_dimension);
				}
				left=backboneid;
			}
			
			//	backbone only has short-ranged interactions (zero)
            
		}
		else{
			for(j=0;j<chainlength[i];j++){
				backboneid=monomerid[i][j].backbone;
				if(my_nonbonded_params.solvationparams.interface==1) calc_onebody_energy_components(coordarray[backboneid], chooseonebody(coordarray[backboneid].nodetype, my_nonbonded_params), my_nonbonded_params.solvationparams, pmy_energycomponents);
				sidechainid=monomerid[i][j].sidechain;
                if(coordarray[sidechainid].nodetype>=0){
                    if(my_nonbonded_params.solvationparams.interface==1) calc_onebody_energy_components(coordarray[sidechainid], chooseonebody(coordarray[sidechainid].nodetype, my_nonbonded_params), my_nonbonded_params.solvationparams, pmy_energycomponents);
                    
                    if(coordarray[sidechainid].nodetype==1) calc_nonbonded_energy_components_cellstruct_leafid_cluster(sidechainid, backboneid, -1, -1, linkset.phenylphenyl.core, coordarray, box_dimension, my_nonbonded_params, pmy_energycomponents, cluster_index_of_chain, size_of_cluster, head_of_cluster, list_of_chain, &Nclusters);
                    if(coordarray[sidechainid].nodetype>1) calc_nonbonded_energy_components_cellstruct_leafid_cluster(sidechainid, backboneid, -1, -1, linkset.chargedcharged.core, coordarray, box_dimension, my_nonbonded_params, pmy_energycomponents, cluster_index_of_chain, size_of_cluster, head_of_cluster, list_of_chain, &Nclusters);
                    
                    if(coordarray[sidechainid].nodetype>=0) (*pmy_energycomponents).sidechain+=calc_sidechain_bonded_energy(my_bonded_params.sidechain[coordarray[sidechainid].nodetype], coordarray[backboneid], coordarray[sidechainid], box_dimension);
                }
			}
        }
	}
    clusterdistribution[0]+=Nchains;
    for(i=0;i<Nclusters;i++){
        if(size_of_cluster[i]==1){
            my_exit("size shouldn't be 1");
        }
        else if(size_of_cluster[i]>1){
            clusterdistribution[size_of_cluster[i]-1]++;
            clusterdistribution[0]-=size_of_cluster[i];
        }
    }
    
    free(cluster_index_of_chain);
    free(size_of_cluster);
    free(head_of_cluster);
    free(list_of_chain);
}

void calc_nonbonded_energy_components_cellstruct_leafid(int n, int neighbor1, int neighbor2, int neighbor3, linkedlist link, coord *coordarray, double_triple box_dimension, nonbonded_params my_nonbonded_params, energycomponents *pmy_energycomponents){
	int i, j, k, i_image, j_image, k_image, index;
	coord coord_image=coordarray[n];
	for(i=-1;i<=1;i++){																	//	loop through neighboring cells
		if((link.cell[n].x==0)&&(i==-1)){
			i_image=link.cellsperside.x-1;
			coord_image.r.x=coordarray[n].r.x+box_dimension.x;
		}
		else if((link.cell[n].x==link.cellsperside.x-1)&&(i==1)){
			i_image=0;
			coord_image.r.x=coordarray[n].r.x-box_dimension.x;
		}
		else{
            coord_image.r.x=coordarray[n].r.x;
            i_image=link.cell[n].x+i;
        }
		for(j=-1;j<=1;j++){
			if((link.cell[n].y==0)&&(j==-1)){
				j_image=link.cellsperside.y-1;
				coord_image.r.y=coordarray[n].r.y+box_dimension.y;
			}
			else if((link.cell[n].y==link.cellsperside.y-1)&&(j==1)){
				j_image=0;
				coord_image.r.y=coordarray[n].r.y-box_dimension.y;
			}
			else{
                coord_image.r.y=coordarray[n].r.y;
                j_image=link.cell[n].y+j;
            }
			for(k=-1;k<=1;k++){
				if((link.cell[n].z==0)&&(k==-1)){
					k_image=link.cellsperside.z-1;
					coord_image.r.z=coordarray[n].r.z+box_dimension.z;
				}
				else if((link.cell[n].z==link.cellsperside.z-1)&&(k==1)){
					k_image=0;
					coord_image.r.z=coordarray[n].r.z-box_dimension.z;
				}
				else{
                    coord_image.r.z=coordarray[n].r.z;
                    k_image=link.cell[n].z+k;
                }
				index=link.head[i_image][j_image][k_image];
				while(index>=0){
					if((index!=n)&&(index!=neighbor1)&&(index!=neighbor2)&&(index!=neighbor3)){					//	don't calculate nonbonded energy for bonded neighbors
						calc_nonbonded_energy_pair_components_leafid(coord_image, coordarray[index], box_dimension, my_nonbonded_params, pmy_energycomponents);
					}
					index=link.list[index];
				}
			}
		}
	}
}

void calc_nonbonded_energy_pair_components_leafid(coord a, coord b, double_triple box_dimension, nonbonded_params params, energycomponents *pmy_energycomponents){
    if(a.leafid==b.leafid){
		if(a.chainid==b.chainid){
			if(a.nodetype==0){
				if(b.nodetype==0){
					(*pmy_energycomponents).nnsamepoly+=hard_energy(a.r, b.r, params.rhard0, params.rhard0, box_dimension);
					return;
				}
				else if(b.nodetype==1){
					(*pmy_energycomponents).nnsamepoly+=hard_energy(a.r, b.r, params.rhard0, params.rhard1, box_dimension);
					return;
				}
				else if(b.nodetype==2){
					(*pmy_energycomponents).cnsamepoly+=hard_energy(a.r, b.r, params.rhard0, params.p2.rhard, box_dimension);
					return;
				}
				else if(b.nodetype==3){
					(*pmy_energycomponents).cnsamepoly+=hard_energy(a.r, b.r, params.rhard0, params.p3.rhard, box_dimension);
					return;
				}
				else my_exit("don't have energy functions for these node types!");
			}
			else if(a.nodetype==1){
				if(b.nodetype==0){
					(*pmy_energycomponents).nnsamepoly+=hard_energy(a.r, b.r, params.rhard1, params.rhard0, box_dimension);
					return;
				}
				else if(b.nodetype==1){
					(*pmy_energycomponents).nnsamepoly+=phenylphenyl_energy(a, b, params.p11vac, box_dimension, params.phenylcutoff2, params.solvationparams);
					return;
				}
				else if(b.nodetype==2){
					(*pmy_energycomponents).cnsamepoly+=hard_energy(a.r, b.r, params.rhard1, params.p2.rhard, box_dimension);
					return;
				}
				else if(b.nodetype==3){
					(*pmy_energycomponents).cnsamepoly+=hard_energy(a.r, b.r, params.rhard1, params.p3.rhard, box_dimension);
					return;
				}
				else my_exit("don't have energy functions for these node types!");
			}
			else if(a.nodetype==2){
				if(b.nodetype==0){
					(*pmy_energycomponents).cnsamepoly+=hard_energy(a.r, b.r, params.p2.rhard, params.rhard0, box_dimension);
					return;
				}
				else if(b.nodetype==1){
					(*pmy_energycomponents).cnsamepoly+=hard_energy(a.r, b.r, params.p2.rhard, params.rhard1, box_dimension);
					return;
				}
				else if(b.nodetype==2){
					(*pmy_energycomponents).cclikesamepoly+=electrostatic_energy(a, b, params.solvationparams, params.p2, params.p2, box_dimension, params.cutoff2);
					return;
				}
				else if(b.nodetype==3){
					(*pmy_energycomponents).ccunlikesamepoly+=electrostatic_energy(a, b, params.solvationparams, params.p2, params.p3, box_dimension, params.cutoff2);
					return;
				}
				else my_exit("don't have energy functions for these node types!");
			}
			else if(a.nodetype==3){
				if(b.nodetype==0){
					(*pmy_energycomponents).cnsamepoly+=hard_energy(a.r, b.r, params.p3.rhard, params.rhard0, box_dimension);
					return;
				}
				else if(b.nodetype==1){
					(*pmy_energycomponents).cnsamepoly+=hard_energy(a.r, b.r, params.p3.rhard, params.rhard1, box_dimension);
					return;
				}
				else if(b.nodetype==2){
					(*pmy_energycomponents).ccunlikesamepoly+=electrostatic_energy(a, b, params.solvationparams, params.p3, params.p2, box_dimension, params.cutoff2);
					return;
				}
				else if(b.nodetype==3){
					(*pmy_energycomponents).cclikesamepoly+=electrostatic_energy(a, b, params.solvationparams, params.p3, params.p3, box_dimension, params.cutoff2);
					return;
				}
				else my_exit("don't have energy functions for these node types!");
			}
			else my_exit("don't have energy functions for these node types!");
		}
		else{
			if(a.nodetype==0){
				if(b.nodetype==0){
					(*pmy_energycomponents).nn+=hard_energy(a.r, b.r, params.rhard0, params.rhard0, box_dimension);
					return;
				}
				else if(b.nodetype==1){
					(*pmy_energycomponents).nn+=hard_energy(a.r, b.r, params.rhard0, params.rhard1, box_dimension);
					return;
				}
				else if(b.nodetype==2){
					(*pmy_energycomponents).cn+=hard_energy(a.r, b.r, params.rhard0, params.p2.rhard, box_dimension);
					return;
				}
				else if(b.nodetype==3){
					(*pmy_energycomponents).cn+=hard_energy(a.r, b.r, params.rhard0, params.p3.rhard, box_dimension);
					return;
				}
				else my_exit("don't have energy functions for these node types!");
			}
			else if(a.nodetype==1){
				if(b.nodetype==0){
					(*pmy_energycomponents).nn+=hard_energy(a.r, b.r, params.rhard1, params.rhard0, box_dimension);
					return;
				}
				else if(b.nodetype==1){
					(*pmy_energycomponents).nn+=phenylphenyl_energy(a, b, params.p11vac, box_dimension, params.phenylcutoff2, params.solvationparams);
					return;
				}
				else if(b.nodetype==2){
					(*pmy_energycomponents).cn+=hard_energy(a.r, b.r, params.rhard1, params.p2.rhard, box_dimension);
					return;
				}
				else if(b.nodetype==3){
					(*pmy_energycomponents).cn+=hard_energy(a.r, b.r, params.rhard1, params.p3.rhard, box_dimension);
					return;
				}
				else my_exit("don't have energy functions for these node types!");
			}
			else if(a.nodetype==2){
				if(b.nodetype==0){
					(*pmy_energycomponents).cn+=hard_energy(a.r, b.r, params.p2.rhard, params.rhard0, box_dimension);
					return;
				}
				else if(b.nodetype==1){
					(*pmy_energycomponents).cn+=hard_energy(a.r, b.r, params.p2.rhard, params.rhard1, box_dimension);
					return;
				}
				else if(b.nodetype==2){
					(*pmy_energycomponents).cclike+=electrostatic_energy(a, b, params.solvationparams, params.p2, params.p2, box_dimension, params.cutoff2);
					return;
				}
				else if(b.nodetype==3){
					(*pmy_energycomponents).ccunlike+=electrostatic_energy(a, b, params.solvationparams, params.p2, params.p3, box_dimension, params.cutoff2);
					return;
				}
				else my_exit("don't have energy functions for these node types!");
			}
			else if(a.nodetype==3){
				if(b.nodetype==0){
					(*pmy_energycomponents).cn+=hard_energy(a.r, b.r, params.p3.rhard, params.rhard0, box_dimension);
					return;
				}
				else if(b.nodetype==1){
					(*pmy_energycomponents).cn+=hard_energy(a.r, b.r, params.p3.rhard, params.rhard1, box_dimension);
					return;
				}
				else if(b.nodetype==2){
					(*pmy_energycomponents).ccunlike+=electrostatic_energy(a, b, params.solvationparams, params.p3, params.p2, box_dimension, params.cutoff2);
					return;
				}
				else if(b.nodetype==3){
					(*pmy_energycomponents).cclike+=electrostatic_energy(a, b, params.solvationparams, params.p3, params.p3, box_dimension, params.cutoff2);
					return;
				}
				else my_exit("don't have energy functions for these node types!");
			}
			else my_exit("don't have energy functions for these node types!");
		}
    }
    else if((a.leafid)/2==(b.leafid)/2){
        if(a.nodetype==0){
            if(b.nodetype==0){
                (*pmy_energycomponents).nncross+=hard_energy(a.r, b.r, params.rhard0, params.rhard0, box_dimension);
                return;
            }
            else if(b.nodetype==1){
                (*pmy_energycomponents).nncross+=hard_energy(a.r, b.r, params.rhard0, params.rhard1, box_dimension);
                return;
            }
            else if(b.nodetype==2){
                (*pmy_energycomponents).cncross+=hard_energy(a.r, b.r, params.rhard0, params.p2.rhard, box_dimension);
                return;
            }
            else if(b.nodetype==3){
                (*pmy_energycomponents).cncross+=hard_energy(a.r, b.r, params.rhard0, params.p3.rhard, box_dimension);
                return;
            }
            else my_exit("don't have energy functions for these node types!");
        }
        else if(a.nodetype==1){
            if(b.nodetype==0){
                (*pmy_energycomponents).nncross+=hard_energy(a.r, b.r, params.rhard1, params.rhard0, box_dimension);
                return;
            }
            else if(b.nodetype==1){
                //(*pmy_energycomponents).nncross+=phenylphenyl_energy(a, b, params.p11vac, box_dimension, params.phenylcutoff2);
                
				(*pmy_energycomponents).nncross+=phenylphenyl_energy(a, b, params.p11vac, box_dimension, params.phenylcutoff2, params.solvationparams);
                
				return;
            }
            else if(b.nodetype==2){
                (*pmy_energycomponents).cncross+=hard_energy(a.r, b.r, params.rhard1, params.p2.rhard, box_dimension);
                return;
            }
            else if(b.nodetype==3){
                (*pmy_energycomponents).cncross+=hard_energy(a.r, b.r, params.rhard1, params.p3.rhard, box_dimension);
                return;
            }
            else my_exit("don't have energy functions for these node types!");
        }
        else if(a.nodetype==2){
            if(b.nodetype==0){
                (*pmy_energycomponents).cncross+=hard_energy(a.r, b.r, params.p2.rhard, params.rhard0, box_dimension);
                return;
            }
            else if(b.nodetype==1){
                (*pmy_energycomponents).cncross+=hard_energy(a.r, b.r, params.p2.rhard, params.rhard1, box_dimension);
                return;
            }
            else if(b.nodetype==2){
                (*pmy_energycomponents).cclikecross+=electrostatic_energy(a, b, params.solvationparams, params.p2, params.p2, box_dimension, params.cutoff2);
                return;
            }
            else if(b.nodetype==3){
                (*pmy_energycomponents).ccunlikecross+=electrostatic_energy(a, b, params.solvationparams, params.p2, params.p3, box_dimension, params.cutoff2);
                return;
            }
            else my_exit("don't have energy functions for these node types!");
        }
        else if(a.nodetype==3){
            if(b.nodetype==0){
                (*pmy_energycomponents).cncross+=hard_energy(a.r, b.r, params.p3.rhard, params.rhard0, box_dimension);
                return;
            }
            else if(b.nodetype==1){
                (*pmy_energycomponents).cncross+=hard_energy(a.r, b.r, params.p3.rhard, params.rhard1, box_dimension);
                return;
            }
            else if(b.nodetype==2){
                (*pmy_energycomponents).ccunlikecross+=electrostatic_energy(a, b, params.solvationparams, params.p3, params.p2, box_dimension, params.cutoff2);
                return;
            }
            else if(b.nodetype==3){
                (*pmy_energycomponents).cclikecross+=electrostatic_energy(a, b, params.solvationparams, params.p3, params.p3, box_dimension, params.cutoff2);
                return;
            }
            else my_exit("don't have energy functions for these node types!");
        }
        else my_exit("don't have energy functions for these node types!");
    }
    else{
        if(a.nodetype==0){
            if(b.nodetype==0){
                (*pmy_energycomponents).nndifferent+=hard_energy(a.r, b.r, params.rhard0, params.rhard0, box_dimension);
                return;
            }
            else if(b.nodetype==1){
                (*pmy_energycomponents).nndifferent+=hard_energy(a.r, b.r, params.rhard0, params.rhard1, box_dimension);
                return;
            }
            else if(b.nodetype==2){
                (*pmy_energycomponents).cndifferent+=hard_energy(a.r, b.r, params.rhard0, params.p2.rhard, box_dimension);
                return;
            }
            else if(b.nodetype==3){
                (*pmy_energycomponents).cndifferent+=hard_energy(a.r, b.r, params.rhard0, params.p3.rhard, box_dimension);
                return;
            }
            else my_exit("don't have energy functions for these node types!");
        }
        else if(a.nodetype==1){
            if(b.nodetype==0){
                (*pmy_energycomponents).nndifferent+=hard_energy(a.r, b.r, params.rhard1, params.rhard0, box_dimension);
                return;
            }
            else if(b.nodetype==1){
				(*pmy_energycomponents).nndifferent+=phenylphenyl_energy(a, b, params.p11vac, box_dimension, params.phenylcutoff2, params.solvationparams);
				return;
            }
            else if(b.nodetype==2){
                (*pmy_energycomponents).cndifferent+=hard_energy(a.r, b.r, params.rhard1, params.p2.rhard, box_dimension);
                return;
            }
            else if(b.nodetype==3){
                (*pmy_energycomponents).cndifferent+=hard_energy(a.r, b.r, params.rhard1, params.p3.rhard, box_dimension);
                return;
            }
            else my_exit("don't have energy functions for these node types!");
        }
        else if(a.nodetype==2){
            if(b.nodetype==0){
                (*pmy_energycomponents).cndifferent+=hard_energy(a.r, b.r, params.p2.rhard, params.rhard0, box_dimension);
                return;
            }
            else if(b.nodetype==1){
                (*pmy_energycomponents).cndifferent+=hard_energy(a.r, b.r, params.p2.rhard, params.rhard1, box_dimension);
                return;
            }
            else if(b.nodetype==2){
                (*pmy_energycomponents).cclikedifferent+=electrostatic_energy(a, b, params.solvationparams, params.p2, params.p2, box_dimension, params.cutoff2);
                return;
            }
            else if(b.nodetype==3){
                (*pmy_energycomponents).ccunlikedifferent+=electrostatic_energy(a, b, params.solvationparams, params.p2, params.p3, box_dimension, params.cutoff2);
                return;
            }
            else my_exit("don't have energy functions for these node types!");
        }
        else if(a.nodetype==3){
            if(b.nodetype==0){
                (*pmy_energycomponents).cndifferent+=hard_energy(a.r, b.r, params.p3.rhard, params.rhard0, box_dimension);
                return;
            }
            else if(b.nodetype==1){
                (*pmy_energycomponents).cndifferent+=hard_energy(a.r, b.r, params.p3.rhard, params.rhard1, box_dimension);
                return;
            }
            else if(b.nodetype==2){
                (*pmy_energycomponents).ccunlikedifferent+=electrostatic_energy(a, b, params.solvationparams, params.p3, params.p2, box_dimension, params.cutoff2);
                return;
            }
            else if(b.nodetype==3){
                (*pmy_energycomponents).cclikedifferent+=electrostatic_energy(a, b, params.solvationparams, params.p3, params.p3, box_dimension, params.cutoff2);
                return;
            }
            else my_exit("don't have energy functions for these node types!");
        }
        else my_exit("don't have energy functions for these node types!");
    }
}

void output_cluster(char *clusterfile, int *clusterdistribution, int cycle, double factor, int Nchains){
	int i;
	FILE *outp;
	outp=fopen(clusterfile, "a");
    fprintf(outp, "%i", cycle);
    for(i=0;i<Nchains;i++){
        if(clusterdistribution[i]>0){
            fprintf(outp, "\t%.4f", clusterdistribution[i]*factor);
            clusterdistribution[i]=0;
        }
        else fprintf(outp, "\t0");
    }
    fprintf(outp, "\n");
    fclose(outp);
}

void output_timeseries(char *totalfilename, energycomponents *pmy_energycomponents, int cycle, monomertypes monomercount, double factor, int Nnodetypes, double_triple *pavboxdimension, int runningtime, double *avmesharea, double *avmesharea2, double *avmeancurvature, double *avmeancurvaturesquared, double *avbondlength, double *avbondlength2, int Ntriangles, double *surfacetensionenergy, double *pboundaryenergy, double *meancurvatureenergy, double *meancurvaturesquaredenergy, double *uniformbondenergy, double *pdvolume, double *pdvolumeenergy){
	int i;
	FILE *outp;
	outp=fopen(totalfilename, "a");
	fprintf(outp, "%i\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f", cycle, (*pmy_energycomponents).backbone*factor/monomercount.monomers, (*pmy_energycomponents).sidechain*factor/monomercount.monomers, 0.5*(*pmy_energycomponents).nnsamepoly*factor/monomercount.monomers, 0.5*(*pmy_energycomponents).cnsamepoly*factor/monomercount.monomers, 0.5*(*pmy_energycomponents).cclikesamepoly*factor/monomercount.monomers, 0.5*(*pmy_energycomponents).ccunlikesamepoly*factor/monomercount.monomers, 0.5*(*pmy_energycomponents).nn*factor/monomercount.monomers, 0.5*(*pmy_energycomponents).cn*factor/monomercount.monomers, 0.5*(*pmy_energycomponents).cclike*factor/monomercount.monomers, 0.5*(*pmy_energycomponents).ccunlike*factor/monomercount.monomers);
    for(i=0;i<2;i++){
		if(monomercount.monomers>0) fprintf(outp, "\t%.8f", surfacetensionenergy[i]*factor/monomercount.monomers);
		else fprintf(outp, "\t%.8f", surfacetensionenergy[i]*factor/(0.5*Ntriangles));
		surfacetensionenergy[i]=0;
	}
    if(monomercount.monomers>0) fprintf(outp, "\t%.8f", (*pboundaryenergy)*factor/monomercount.monomers);
    else fprintf(outp, "\t%.8f", (*pboundaryenergy)*factor/(0.5*Ntriangles));
    (*pboundaryenergy)=0;
    for(i=0;i<2;i++){
		if(monomercount.monomers>0) fprintf(outp, "\t%.8f", meancurvatureenergy[i]*factor/monomercount.monomers);
		else fprintf(outp, "\t%.8f", meancurvatureenergy[i]*factor/(0.5*Ntriangles));
		meancurvatureenergy[i]=0;
	}
    for(i=0;i<2;i++){
		if(monomercount.monomers>0) fprintf(outp, "\t%.8f", meancurvaturesquaredenergy[i]*factor/monomercount.monomers);
		else fprintf(outp, "\t%.8f", meancurvaturesquaredenergy[i]*factor/(0.5*Ntriangles));
		meancurvaturesquaredenergy[i]=0;
	}
	for(i=0;i<6;i++){
		if(monomercount.monomers>0){
            fprintf(outp, "\t%.8f", uniformbondenergy[i]*factor/monomercount.monomers);          //  report per monomer, not per mesh point, so it compares to everything else
        }
		else fprintf(outp, "\t%.8f", uniformbondenergy[i]*factor/(0.5*Ntriangles));             //  if bare interface, report per mesh point
        uniformbondenergy[i]=0;
	}
	if(monomercount.monomers>0) fprintf(outp, "\t%.8f", (*pdvolumeenergy)*factor/monomercount.monomers);
	else fprintf(outp, "\t%.8f", (*pdvolumeenergy)*factor/(0.5*Ntriangles));
	for(i=0;i<Nnodetypes;i++){
        fprintf(outp, "\t%.8f\t%.8f", (*pmy_energycomponents).solvation[i]*factor/monomercount.monomers, (*pmy_energycomponents).interface[i]*factor/monomercount.monomers);
		(*pmy_energycomponents).solvation[i]=0;
		(*pmy_energycomponents).interface[i]=0;
    }
    
    for(i=0;i<2;i++){
        avmesharea[i]*=(factor/Ntriangles);     //  per triangle (Ntriangles is number of triangles per surface)
        avmesharea2[i]*=(factor/Ntriangles);
    }
    for(i=0;i<6;i++){
        avbondlength[i]*=(factor/(0.5*Ntriangles));     //  per bond (Ntriangles/2 bonds in each of three directions for each surface)
        avbondlength2[i]*=(factor/(0.5*Ntriangles));
    }
    
    //  Don't need to divide meancurvature and meancurvaturesquared because I already divided these area-integrated values by the area
    
	fprintf(outp, "\t%.8f\t%.8f\t%.8f", (*pavboxdimension).x*factor, (*pavboxdimension).y*factor, (*pavboxdimension).z*factor);
	(*pavboxdimension).x=(*pavboxdimension).y=(*pavboxdimension).z=0;
    fprintf(outp, "\t%.8f\t%.8f\t%.8f\t%.8f\t%i", avmesharea[0], avmesharea[1], sqrt(avmesharea2[0]-avmesharea[0]*avmesharea[0]), sqrt(avmesharea2[1]-avmesharea[1]*avmesharea[1]), Ntriangles);
    for(i=0;i<6;i++){
        fprintf(outp, "\t%.8f", avbondlength[i]);
    }
    for(i=0;i<6;i++){
        fprintf(outp, "\t%.8f", sqrt(avbondlength2[i]-avbondlength[i]*avbondlength[i]));
    }
    for(i=0;i<2;i++){
        fprintf(outp, "\t%.8f", avmeancurvature[i]*factor);
    }
    for(i=0;i<2;i++){
        fprintf(outp, "\t%.8f", avmeancurvaturesquared[i]*factor);
    }
    fprintf(outp, "\t%i", Ntriangles/2);
	fprintf(outp, "\t%.8f", (*pdvolume)*factor/Ntriangles);
	fprintf(outp, "\t%i\n", runningtime);
	fclose(outp);
	(*pmy_energycomponents).backbone=0;
	(*pmy_energycomponents).sidechain=0;
	(*pmy_energycomponents).nn=0;
	(*pmy_energycomponents).ccunlike=0;
	(*pmy_energycomponents).cclike=0;
	(*pmy_energycomponents).cn=0;
	(*pmy_energycomponents).nnsamepoly=0;
	(*pmy_energycomponents).ccunlikesamepoly=0;
	(*pmy_energycomponents).cclikesamepoly=0;
	(*pmy_energycomponents).cnsamepoly=0;
    for(i=0;i<2;i++){
        avmesharea[i]=avmesharea2[i]=0;
        avbondlength[i]=avbondlength2[i]=0;
        avmeancurvature[i]=avmeancurvaturesquared[i]=0;
    }
	(*pdvolume)=(*pdvolumeenergy)=0;
}

void output_timeseries_leaves(char *totalfilename, energycomponents *pmy_energycomponents, int cycle, monomertypes monomercount, double factor, int Nnodetypes, double_triple *pavboxdimension, int runningtime){
	int i;
	FILE *outp;
	outp=fopen(totalfilename, "a");
	fprintf(outp, "%i\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f", cycle, (*pmy_energycomponents).backbone*factor/monomercount.monomers, (*pmy_energycomponents).sidechain*factor/monomercount.monomers, 0.5*(*pmy_energycomponents).nnsamepoly*factor/monomercount.monomers, 0.5*(*pmy_energycomponents).cnsamepoly*factor/monomercount.monomers, 0.5*(*pmy_energycomponents).cclikesamepoly*factor/monomercount.monomers, 0.5*(*pmy_energycomponents).ccunlikesamepoly*factor/monomercount.monomers, 0.5*(*pmy_energycomponents).nn*factor/monomercount.monomers, 0.5*(*pmy_energycomponents).cn*factor/monomercount.monomers, 0.5*(*pmy_energycomponents).cclike*factor/monomercount.monomers, 0.5*(*pmy_energycomponents).ccunlike*factor/monomercount.monomers, 0.5*(*pmy_energycomponents).nncross*factor/monomercount.monomers, 0.5*(*pmy_energycomponents).cncross*factor/monomercount.monomers, 0.5*(*pmy_energycomponents).cclikecross*factor/monomercount.monomers, 0.5*(*pmy_energycomponents).ccunlikecross*factor/monomercount.monomers);
	for(i=0;i<Nnodetypes;i++){
        fprintf(outp, "\t%.8f\t%.8f", (*pmy_energycomponents).solvation[i]*factor/monomercount.monomers, (*pmy_energycomponents).interface[i]*factor/monomercount.monomers);
		(*pmy_energycomponents).solvation[i]=0;
		(*pmy_energycomponents).interface[i]=0;
    }
	fprintf(outp, "\t%.8f\t%.8f\t%.8f", (*pavboxdimension).x*factor, (*pavboxdimension).y*factor, (*pavboxdimension).z*factor);
	(*pavboxdimension).x=(*pavboxdimension).y=(*pavboxdimension).z=0;
	fprintf(outp, "\t%i\n", runningtime);
	fclose(outp);
	(*pmy_energycomponents).backbone=0;
	(*pmy_energycomponents).sidechain=0;
	(*pmy_energycomponents).nn=0;
	(*pmy_energycomponents).ccunlike=0;
	(*pmy_energycomponents).cclike=0;
	(*pmy_energycomponents).cn=0;
    
	(*pmy_energycomponents).nncross=0;
	(*pmy_energycomponents).ccunlikecross=0;
	(*pmy_energycomponents).cclikecross=0;
	(*pmy_energycomponents).cncross=0;
    
	(*pmy_energycomponents).nndifferent=0;
	(*pmy_energycomponents).ccunlikedifferent=0;
	(*pmy_energycomponents).cclikedifferent=0;
	(*pmy_energycomponents).cndifferent=0;
	
	(*pmy_energycomponents).nnsamepoly=0;
	(*pmy_energycomponents).ccunlikesamepoly=0;
	(*pmy_energycomponents).cclikesamepoly=0;
	(*pmy_energycomponents).cnsamepoly=0;
}

void output_timeseries_sheets(char *totalfilename, energycomponents *pmy_energycomponents, int cycle, monomertypes monomercount, double factor, int Nnodetypes, double_triple *pavboxdimension, int runningtime){
	int i;
	FILE *outp;
	outp=fopen(totalfilename, "a");
	fprintf(outp, "%i\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f", cycle, (*pmy_energycomponents).backbone*factor/monomercount.monomers, (*pmy_energycomponents).sidechain*factor/monomercount.monomers, 0.5*(*pmy_energycomponents).nnsamepoly*factor/monomercount.monomers, 0.5*(*pmy_energycomponents).cnsamepoly*factor/monomercount.monomers, 0.5*(*pmy_energycomponents).cclikesamepoly*factor/monomercount.monomers, 0.5*(*pmy_energycomponents).ccunlikesamepoly*factor/monomercount.monomers, 0.5*(*pmy_energycomponents).nn*factor/monomercount.monomers, 0.5*(*pmy_energycomponents).cn*factor/monomercount.monomers, 0.5*(*pmy_energycomponents).cclike*factor/monomercount.monomers, 0.5*(*pmy_energycomponents).ccunlike*factor/monomercount.monomers, 0.5*(*pmy_energycomponents).nncross*factor/monomercount.monomers, 0.5*(*pmy_energycomponents).cncross*factor/monomercount.monomers, 0.5*(*pmy_energycomponents).cclikecross*factor/monomercount.monomers, 0.5*(*pmy_energycomponents).ccunlikecross*factor/monomercount.monomers, 0.5*(*pmy_energycomponents).nndifferent*factor/monomercount.monomers, 0.5*(*pmy_energycomponents).cndifferent*factor/monomercount.monomers, 0.5*(*pmy_energycomponents).cclikedifferent*factor/monomercount.monomers, 0.5*(*pmy_energycomponents).ccunlikedifferent*factor/monomercount.monomers);
	for(i=0;i<Nnodetypes;i++){
        fprintf(outp, "\t%.8f\t%.8f", (*pmy_energycomponents).solvation[i]*factor/monomercount.monomers, (*pmy_energycomponents).interface[i]*factor/monomercount.monomers);
		(*pmy_energycomponents).solvation[i]=0;
		(*pmy_energycomponents).interface[i]=0;
    }
	fprintf(outp, "\t%.8f\t%.8f\t%.8f", (*pavboxdimension).x*factor, (*pavboxdimension).y*factor, (*pavboxdimension).z*factor);
	(*pavboxdimension).x=(*pavboxdimension).y=(*pavboxdimension).z=0;
	fprintf(outp, "\t%i\n", runningtime);
	fclose(outp);
	(*pmy_energycomponents).backbone=0;
	(*pmy_energycomponents).sidechain=0;
	(*pmy_energycomponents).nn=0;
	(*pmy_energycomponents).ccunlike=0;
	(*pmy_energycomponents).cclike=0;
	(*pmy_energycomponents).cn=0;
    
	(*pmy_energycomponents).nncross=0;
	(*pmy_energycomponents).ccunlikecross=0;
	(*pmy_energycomponents).cclikecross=0;
	(*pmy_energycomponents).cncross=0;
    
	(*pmy_energycomponents).nndifferent=0;
	(*pmy_energycomponents).ccunlikedifferent=0;
	(*pmy_energycomponents).cclikedifferent=0;
	(*pmy_energycomponents).cndifferent=0;
	
	(*pmy_energycomponents).nnsamepoly=0;
	(*pmy_energycomponents).ccunlikesamepoly=0;
	(*pmy_energycomponents).cclikesamepoly=0;
	(*pmy_energycomponents).cnsamepoly=0;	
}

void output_trajectory(int Nchains, int *chainlength, coord *coordarray, monomernodes **monomerid, double_triple box_dimension, char *filename, int *pframe, int cycle){
	int i, j;
    FILE *outp;
    if((*pframe)==0) outp=fopen(filename, "w");
    else outp=fopen(filename, "a");
    for(i=0;i<Nchains;i++){
        for(j=0;j<chainlength[i];j++){
            fprintf(outp, "%i %i %i %i %.6f %.6f %.6f %.6f %.6f %.6f\n", coordarray[monomerid[i][j].backbone].leafid, i, j, coordarray[monomerid[i][j].backbone].nodetype, coordarray[monomerid[i][j].backbone].r.x, coordarray[monomerid[i][j].backbone].r.y, coordarray[monomerid[i][j].backbone].r.z, coordarray[monomerid[i][j].backbone].n.x, coordarray[monomerid[i][j].backbone].n.y, coordarray[monomerid[i][j].backbone].n.z);
            if(coordarray[monomerid[i][j].sidechain].nodetype>=0) fprintf(outp, "%i %i %i %i %.6f %.6f %.6f %.6f %.6f %.6f\n", coordarray[monomerid[i][j].sidechain].leafid, i, j, coordarray[monomerid[i][j].sidechain].nodetype, coordarray[monomerid[i][j].sidechain].r.x, coordarray[monomerid[i][j].sidechain].r.y, coordarray[monomerid[i][j].sidechain].r.z, coordarray[monomerid[i][j].sidechain].n.x, coordarray[monomerid[i][j].sidechain].n.y, coordarray[monomerid[i][j].sidechain].n.z);
            else fprintf(outp, "%i %i %i %i 0 0 0 0 0 0\n", coordarray[monomerid[i][j].sidechain].leafid, i, j, coordarray[monomerid[i][j].sidechain].nodetype);
        }
    }
    fprintf(outp, "%i box %.6f %.6f %.6f\n\n", cycle, box_dimension.x, box_dimension.y, box_dimension.z);
    fclose(outp);
    (*pframe)++;
}

void calc_displacement(coord left, coord right, displacementparams *pdisplacement, double boxwidth){
    double sep=right.r.x-left.r.x;
	double dif=sep-(*pdisplacement).lastsep;
	recenter(dif, boxwidth);
	(*pdisplacement).lastsep+=dif;
	sep=(*pdisplacement).lastsep;
    (*pdisplacement).av+=sep;
    (*pdisplacement).count++;
}

void output_displacement_timeseries(char *filename, int cycle, displacementparams *pdisplacement){
	(*pdisplacement).av/=(1.*(*pdisplacement).count);
	FILE *outp;
	outp=fopen(filename, "a");
	fprintf(outp, "%i\t%.8f\n", cycle, (*pdisplacement).av);
	(*pdisplacement).av=(*pdisplacement).count=0;
	fclose(outp);
}

double calc_sidechain_bonded_energy(sidechain_params sidechain, coord backbonecoord, coord sidechaincoord, double_triple box_dimension){
 	double_triple sep=subtract_double_triple(sidechaincoord.r, backbonecoord.r);
	recenter_double_triple(&sep, box_dimension);
	double norm2=dot_product(sep, sep);
	double rparallel=dot_product(sep, backbonecoord.n);
	double rperp=sqrt(norm2-rparallel*rparallel);
    double npar=dot_product(sidechaincoord.n, backbonecoord.n);
	double rdotn=dot_product(sep, sidechaincoord.n);
    
	double energy=sidechain.k1*pow(sqrt(pow(rperp-sidechain.rperp0, 2)+pow(rparallel-sidechain.rpar0, 2))-sidechain.r0, 2);
	double arctan=atan(rparallel/rperp);
	energy+=(sidechain.J10+arctan*(sidechain.J11+arctan*(sidechain.J12+arctan*(sidechain.J13+arctan*sidechain.J14))));
	if(sidechaincoord.orientationtype==0){              //  parallel
		energy+=(sidechain.k2*pow(npar-(sidechain.r20+rparallel*(sidechain.r21+rparallel*sidechain.r22)), 2));
		energy+=(sidechain.k3*pow(rdotn-(sidechain.r30+npar*(sidechain.r31+npar*sidechain.r32)), 2));
	}
	if(sidechaincoord.orientationtype==1){              //  perpendicular
        energy+=(rparallel*(sidechain.J201+rparallel*(sidechain.J202+rparallel*(sidechain.J203+rparallel*(sidechain.J204+rparallel*(sidechain.J205+rparallel*sidechain.J206))))));
		double npar2=npar*npar;
		double npar4=npar2*npar2;
		energy+=(npar2*(sidechain.J220+rparallel*(sidechain.J221+rparallel*sidechain.J222)));
		energy+=(npar4*sidechain.J240);
		
		//	J30x are redundant, 300, 302, 304 added to 200, 220, and 240, respectively
		
		energy+=(rdotn*npar*(sidechain.J311+npar2*sidechain.J313));
		double rdotnpower=rdotn*rdotn;
		energy+=(rdotnpower*(sidechain.J320+npar2*sidechain.J322));
		rdotnpower*=rdotn;
		energy+=(rdotnpower*npar*sidechain.J331);
		rdotnpower*=rdotn;
		energy+=(rdotnpower*sidechain.J340);
	}
	return energy;
}

double calc_backbone_bonded_energy(bonded_params my_bonded_params, coord centercoord, coord leftcoord, coord rightcoord, double_triple box_dimension){
	double_triple riip1=subtract_double_triple(rightcoord.r, centercoord.r);
	recenter_double_triple(&riip1, box_dimension);
	double_triple rim1i=subtract_double_triple(centercoord.r, leftcoord.r);
	recenter_double_triple(&rim1i, box_dimension);
	double riip1norm=norm(riip1);
	double rim1inorm=norm(rim1i);
	double_triple rhatiip1=scalar_multiply_double_triple(riip1, 1./riip1norm);
	double_triple rhatim1i=scalar_multiply_double_triple(rim1i, 1./rim1inorm);
	
    double energy=my_bonded_params.K10+riip1norm*(my_bonded_params.K11+riip1norm*(my_bonded_params.K12+riip1norm*(my_bonded_params.K13+riip1norm*my_bonded_params.K14)));	//	onlyright
    
    energy+=(my_bonded_params.K10+rim1inorm*(my_bonded_params.K11+rim1inorm*(my_bonded_params.K12+rim1inorm*(my_bonded_params.K13+rim1inorm*my_bonded_params.K14))));	//	onlyleft
	
	double riip1dotni=dot_product(rhatiip1, centercoord.n);
	double riip1dotnip1=dot_product(rhatiip1, rightcoord.n);
	double rim1idotni=dot_product(rhatim1i, centercoord.n);
	double rim1idotnim1=dot_product(rhatim1i, leftcoord.n);
    
	energy+=(my_bonded_params.kr*pow(riip1dotni-(riip1norm-my_bonded_params.r0r)/my_bonded_params.sr, 2));			//	onlyright, right interaction
	energy+=(my_bonded_params.kl*pow(riip1dotnip1-(riip1norm-my_bonded_params.r0l)/my_bonded_params.sl, 2));		//	onlyright, left interaction
	energy+=(my_bonded_params.kl*pow(rim1idotni-(rim1inorm-my_bonded_params.r0l)/my_bonded_params.sl, 2));			//	onlyright, left interaction
	energy+=(my_bonded_params.kr*pow(rim1idotnim1-(rim1inorm-my_bonded_params.r0r)/my_bonded_params.sr, 2));		//	onlyright, right interaction
	
	double dotcrossproduct=dot_product(rhatim1i, rhatiip1)-rim1idotni*riip1dotni;
	
	energy+=(my_bonded_params.K20+dotcrossproduct*(my_bonded_params.K21+dotcrossproduct*(my_bonded_params.K22+dotcrossproduct*(my_bonded_params.K23+dotcrossproduct*my_bonded_params.K24))));
    
	return my_bonded_params.factor*energy;
}

double calc_backbone_bonded_energy_norim1i(bonded_params my_bonded_params, coord centercoord, coord leftcoord, coord rightcoord, double_triple box_dimension){
	double_triple riip1=subtract_double_triple(rightcoord.r, centercoord.r);
	recenter_double_triple(&riip1, box_dimension);
	double_triple rim1i=subtract_double_triple(centercoord.r, leftcoord.r);
	recenter_double_triple(&rim1i, box_dimension);
	double riip1norm=norm(riip1);
	double rim1inorm=norm(rim1i);
	double_triple rhatiip1=scalar_multiply_double_triple(riip1, 1./riip1norm);
	double_triple rhatim1i=scalar_multiply_double_triple(rim1i, 1./rim1inorm);
	
    double energy=my_bonded_params.K10+riip1norm*(my_bonded_params.K11+riip1norm*(my_bonded_params.K12+riip1norm*(my_bonded_params.K13+riip1norm*my_bonded_params.K14)));	//	onlyright
    	
	double riip1dotni=dot_product(rhatiip1, centercoord.n);
	double riip1dotnip1=dot_product(rhatiip1, rightcoord.n);
	double rim1idotni=dot_product(rhatim1i, centercoord.n);
	
	energy+=(my_bonded_params.kr*pow(riip1dotni-(riip1norm-my_bonded_params.r0r)/my_bonded_params.sr, 2));			//	onlyright, right interaction
    energy+=(my_bonded_params.kl*pow(riip1dotnip1-(riip1norm-my_bonded_params.r0l)/my_bonded_params.sl, 2));		//	onlyright, left interaction
	
	double dotcrossproduct=dot_product(rhatim1i, rhatiip1)-rim1idotni*riip1dotni;
	
	energy+=(my_bonded_params.K20+dotcrossproduct*(my_bonded_params.K21+dotcrossproduct*(my_bonded_params.K22+dotcrossproduct*(my_bonded_params.K23+dotcrossproduct*my_bonded_params.K24))));
	
	return my_bonded_params.factor*energy;
}

double calc_backbone_bonded_energy_onlyleft(bonded_params my_bonded_params, coord centercoord, coord leftcoord, double_triple box_dimension){
	double_triple rim1i=subtract_double_triple(centercoord.r, leftcoord.r);
	recenter_double_triple(&rim1i, box_dimension);
	double rim1inorm=norm(rim1i);
	double_triple rhatim1i=scalar_multiply_double_triple(rim1i, 1./rim1inorm);
	
    double energy=my_bonded_params.K10+rim1inorm*(my_bonded_params.K11+rim1inorm*(my_bonded_params.K12+rim1inorm*(my_bonded_params.K13+rim1inorm*my_bonded_params.K14)));	//	onlyleft
	
	double rim1idotni=dot_product(rhatim1i, centercoord.n);
	double rim1idotnim1=dot_product(rhatim1i, leftcoord.n);
    
	energy+=(my_bonded_params.kl*pow(rim1idotni-(rim1inorm-my_bonded_params.r0l)/my_bonded_params.sl, 2));			//	onlyleft, left interaction
	energy+=(my_bonded_params.kr*pow(rim1idotnim1-(rim1inorm-my_bonded_params.r0r)/my_bonded_params.sr, 2));		//	onlyleft, right interaction
    
	return my_bonded_params.factor*energy;
}

double calc_energy_difference_changer_hard(int n, coord coord_old, coord coord_new, coord *coordarray, double_triple box_dimension, bonded_params my_bonded_params, nonbonded_params my_nonbonded_params, monomernodes **monomerid, int *chainlength, linkedlist link){
	int moverchain=coord_new.chainid;
	int movermonomer=coord_new.monomerid;
	int neighbor, twoleftneighbor=-1, leftneighbor=-1, rightneighbor=-1, tworightneighbor=-1;
	double difference=0;
	if(coord_new.nodetype==0){		//	backbone
		neighbor=monomerid[moverchain][movermonomer].sidechain;
        if(coordarray[neighbor].nodetype>=0){
			difference=calc_sidechain_bonded_energy(my_bonded_params.sidechain[coordarray[neighbor].nodetype], coord_new, coordarray[neighbor], box_dimension);
			difference-=calc_sidechain_bonded_energy(my_bonded_params.sidechain[coordarray[neighbor].nodetype], coord_old, coordarray[neighbor], box_dimension);
        }
		if(movermonomer>0){
			leftneighbor=monomerid[moverchain][movermonomer-1].backbone;
			if(movermonomer>1) twoleftneighbor=monomerid[moverchain][movermonomer-2].backbone;
		}
		if(movermonomer<chainlength[moverchain]-1){
			rightneighbor=monomerid[moverchain][movermonomer+1].backbone;
			if(movermonomer<chainlength[moverchain]-2) tworightneighbor=monomerid[moverchain][movermonomer+2].backbone;
		}
		difference+=calc_backbone_bonded_difference_fourneighbors_changer(twoleftneighbor, leftneighbor, rightneighbor, tworightneighbor, my_bonded_params, coordarray, coord_old, coord_new, box_dimension);
		difference+=calc_hard_energy_cellstruct(neighbor, leftneighbor, rightneighbor, link.number_neighbors[n], link.neighbor[n], coord_new, coordarray, box_dimension, my_nonbonded_params);
	}
	else if(coord_new.nodetype>0){							//	sidechain
		neighbor=monomerid[moverchain][movermonomer].backbone;
        if(coord_new.nodetype>=0){
            difference=calc_sidechain_bonded_energy(my_bonded_params.sidechain[coord_new.nodetype], coordarray[neighbor], coord_new, box_dimension);
            difference-=calc_sidechain_bonded_energy(my_bonded_params.sidechain[coord_new.nodetype], coordarray[neighbor], coord_old, box_dimension);
        }
		difference+=calc_hard_energy_cellstruct(neighbor, -1, -1, link.number_neighbors[n], link.neighbor[n], coord_new, coordarray, box_dimension, my_nonbonded_params);
	}
    else my_exit("trying to move a site with nodetype<0! (calc_energy_difference_changer_hard)");
	if(my_nonbonded_params.solvationparams.interface==1) difference+=calc_onebody_difference(coord_old, coord_new, chooseonebody(coord_new.nodetype, my_nonbonded_params), my_nonbonded_params.solvationparams);
	return my_bonded_params.factor*difference;
}

double calc_energy_difference_changer_phenyl(int n, coord coord_old, coord coord_new, coord *coordarray, double_triple box_dimension, bonded_params my_bonded_params, nonbonded_params my_nonbonded_params, monomernodes **monomerid, int *chainlength, linkedlist link1, linkedlist link2){
	int i;
	int moverchain=coord_new.chainid;
	int movermonomer=coord_new.monomerid;
	int neighbor, twoleftneighbor=-1, leftneighbor=-1, rightneighbor=-1, tworightneighbor=-1;
	double difference=0;
	if(coord_new.nodetype==0){		//	backbone
		neighbor=monomerid[moverchain][movermonomer].sidechain;
        if(coordarray[neighbor].nodetype>=0){
            difference=calc_sidechain_bonded_energy(my_bonded_params.sidechain[coordarray[neighbor].nodetype], coord_new, coordarray[neighbor], box_dimension);
            difference-=calc_sidechain_bonded_energy(my_bonded_params.sidechain[coordarray[neighbor].nodetype], coord_old, coordarray[neighbor], box_dimension);
        }
		if(movermonomer>0){
			leftneighbor=monomerid[moverchain][movermonomer-1].backbone;
			if(movermonomer>1) twoleftneighbor=monomerid[moverchain][movermonomer-2].backbone;
		}
		if(movermonomer<chainlength[moverchain]-1){
			rightneighbor=monomerid[moverchain][movermonomer+1].backbone;
			if(movermonomer<chainlength[moverchain]-2) tworightneighbor=monomerid[moverchain][movermonomer+2].backbone;
		}
		difference+=calc_backbone_bonded_difference_fourneighbors_changer(twoleftneighbor, leftneighbor, rightneighbor, tworightneighbor, my_bonded_params, coordarray, coord_old, coord_new, box_dimension);
		difference+=calc_hard_energy_cellstruct(neighbor, leftneighbor, rightneighbor, link1.number_neighbors[n], link1.neighbor[n], coord_new, coordarray, box_dimension, my_nonbonded_params);
		difference+=calc_phenyl_energy_cellstruct(neighbor, leftneighbor, rightneighbor, link2.number_neighbors[n], link2.neighbor[n], coord_new, coordarray, box_dimension, my_nonbonded_params);
		for(i=0;i<link2.number_neighbors[n];i++) difference-=link2.pair_energy[n][i];
	}
	else if(coord_new.nodetype>0){							//	sidechain
		neighbor=monomerid[moverchain][movermonomer].backbone;
        if(coord_new.nodetype>=0){
            difference=calc_sidechain_bonded_energy(my_bonded_params.sidechain[coord_new.nodetype], coordarray[neighbor], coord_new, box_dimension);
            difference-=calc_sidechain_bonded_energy(my_bonded_params.sidechain[coord_new.nodetype], coordarray[neighbor], coord_old, box_dimension);
        }
		difference+=calc_hard_energy_cellstruct(neighbor, -1, -1, link1.number_neighbors[n], link1.neighbor[n], coord_new, coordarray, box_dimension, my_nonbonded_params);
		difference+=calc_phenyl_energy_cellstruct(neighbor, -1, -1, link2.number_neighbors[n], link2.neighbor[n], coord_new, coordarray, box_dimension, my_nonbonded_params);
		for(i=0;i<link2.number_neighbors[n];i++) difference-=link2.pair_energy[n][i];
	}
    else my_exit("trying to move a site with nodetype<0! (calc_energy_difference_changer_phenyl)");
	if(my_nonbonded_params.solvationparams.interface==1) difference+=calc_onebody_difference(coord_old, coord_new, chooseonebody(coord_new.nodetype, my_nonbonded_params), my_nonbonded_params.solvationparams);
	return difference;
}

double calc_energy_difference_changer_charged(int n, coord coord_old, coord coord_new, coord *coordarray, double_triple box_dimension, bonded_params my_bonded_params, nonbonded_params my_nonbonded_params, monomernodes **monomerid, int *chainlength, linkedlist link1, linkedlist link2){
	int i;
	int moverchain=coord_new.chainid;
	int movermonomer=coord_new.monomerid;
	int neighbor, twoleftneighbor=-1, leftneighbor=-1, rightneighbor=-1, tworightneighbor=-1;
	double difference=0;
	if(coord_new.nodetype==0){		//	backbone
		neighbor=monomerid[moverchain][movermonomer].sidechain;
        if(coordarray[neighbor].nodetype>=0){
            difference=calc_sidechain_bonded_energy(my_bonded_params.sidechain[coordarray[neighbor].nodetype], coord_new, coordarray[neighbor], box_dimension);
            difference-=calc_sidechain_bonded_energy(my_bonded_params.sidechain[coordarray[neighbor].nodetype], coord_old, coordarray[neighbor], box_dimension);
        }
		if(movermonomer>0){
			leftneighbor=monomerid[moverchain][movermonomer-1].backbone;
			if(movermonomer>1) twoleftneighbor=monomerid[moverchain][movermonomer-2].backbone;
		}
		if(movermonomer<chainlength[moverchain]-1){
			rightneighbor=monomerid[moverchain][movermonomer+1].backbone;
			if(movermonomer<chainlength[moverchain]-2) tworightneighbor=monomerid[moverchain][movermonomer+2].backbone;
		}
		difference+=calc_backbone_bonded_difference_fourneighbors_changer(twoleftneighbor, leftneighbor, rightneighbor, tworightneighbor, my_bonded_params, coordarray, coord_old, coord_new, box_dimension);
		difference+=calc_hard_energy_cellstruct(neighbor, leftneighbor, rightneighbor, link1.number_neighbors[n], link1.neighbor[n], coord_new, coordarray, box_dimension, my_nonbonded_params);
		difference+=calc_charged_energy_cellstruct(neighbor, leftneighbor, rightneighbor, link2.number_neighbors[n], link2.neighbor[n], coord_new, coordarray, box_dimension, my_nonbonded_params);
		for(i=0;i<link2.number_neighbors[n];i++) difference-=link2.pair_energy[n][i];
	}
	else if(coord_new.nodetype>0){							//	sidechain
		neighbor=monomerid[moverchain][movermonomer].backbone;
        if(coord_new.nodetype>=0){
            difference=calc_sidechain_bonded_energy(my_bonded_params.sidechain[coord_new.nodetype], coordarray[neighbor], coord_new, box_dimension);
            difference-=calc_sidechain_bonded_energy(my_bonded_params.sidechain[coord_new.nodetype], coordarray[neighbor], coord_old, box_dimension);
        }
		difference+=calc_hard_energy_cellstruct(neighbor, -1, -1, link1.number_neighbors[n], link1.neighbor[n], coord_new, coordarray, box_dimension, my_nonbonded_params);
		difference+=calc_charged_energy_cellstruct(neighbor, -1, -1, link2.number_neighbors[n], link2.neighbor[n], coord_new, coordarray, box_dimension, my_nonbonded_params);
		for(i=0;i<link2.number_neighbors[n];i++) difference-=link2.pair_energy[n][i];
	}
    else my_exit("trying to move a site with nodetype<0! (calc_energy_difference_changer_charged)");
	if(my_nonbonded_params.solvationparams.interface==1) difference+=calc_onebody_difference(coord_old, coord_new, chooseonebody(coord_new.nodetype, my_nonbonded_params), my_nonbonded_params.solvationparams);
	return difference;
}

double calc_energy_difference_changen_hard(int n, coord coord_old, coord coord_new, coord *coordarray, double_triple box_dimension, bonded_params my_bonded_params, nonbonded_params my_nonbonded_params, monomernodes **monomerid, int *chainlength, linkedlist link){
	int moverchain=coord_new.chainid;
	int movermonomer=coord_new.monomerid;
	int neighbor, leftneighbor=-1, rightneighbor=-1;
	double difference=0;
	if(coord_new.nodetype==0){		//	backbone
		neighbor=monomerid[moverchain][movermonomer].sidechain;
        if(coordarray[neighbor].nodetype>=0){
            difference=calc_sidechain_bonded_energy(my_bonded_params.sidechain[coordarray[neighbor].nodetype], coord_new, coordarray[neighbor], box_dimension);
            difference-=calc_sidechain_bonded_energy(my_bonded_params.sidechain[coordarray[neighbor].nodetype], coord_old, coordarray[neighbor], box_dimension);
        }
		if(movermonomer>0) leftneighbor=monomerid[moverchain][movermonomer-1].backbone;
		if(movermonomer<chainlength[moverchain]-1) rightneighbor=monomerid[moverchain][movermonomer+1].backbone;
		difference+=calc_backbone_bonded_difference_twoneighbors_changen(leftneighbor, rightneighbor, my_bonded_params, coordarray, coord_old, coord_new, box_dimension);
	}
	else if(coord_new.nodetype>0){							//	sidechain
		neighbor=monomerid[moverchain][movermonomer].backbone;
        if(coord_new.nodetype>=0){
            difference=calc_sidechain_bonded_energy(my_bonded_params.sidechain[coord_new.nodetype], coordarray[neighbor], coord_new, box_dimension);
            difference-=calc_sidechain_bonded_energy(my_bonded_params.sidechain[coord_new.nodetype], coordarray[neighbor], coord_old, box_dimension);
        }
	}
    else my_exit("trying to move a site with nodetype<0! (calc_energy_difference_changen_hard)");
	if(my_nonbonded_params.solvationparams.interface==1) difference+=calc_onebody_difference(coord_old, coord_new, chooseonebody(coord_new.nodetype, my_nonbonded_params), my_nonbonded_params.solvationparams);
	return difference;
}

double calc_energy_difference_changen_phenyl(int n, coord coord_old, coord coord_new, coord *coordarray, double_triple box_dimension, bonded_params my_bonded_params, nonbonded_params my_nonbonded_params, monomernodes **monomerid, int *chainlength, linkedlist link){
	int i;
	int moverchain=coord_new.chainid;
	int movermonomer=coord_new.monomerid;
	int neighbor, leftneighbor=-1, rightneighbor=-1;
	double difference=0;
	if(coord_new.nodetype==0){		//	backbone
		neighbor=monomerid[moverchain][movermonomer].sidechain;
        if(coordarray[neighbor].nodetype>=0){
            difference=calc_sidechain_bonded_energy(my_bonded_params.sidechain[coordarray[neighbor].nodetype], coord_new, coordarray[neighbor], box_dimension);
            difference-=calc_sidechain_bonded_energy(my_bonded_params.sidechain[coordarray[neighbor].nodetype], coord_old, coordarray[neighbor], box_dimension);
        }
		if(movermonomer>0) leftneighbor=monomerid[moverchain][movermonomer-1].backbone;
		if(movermonomer<chainlength[moverchain]-1) rightneighbor=monomerid[moverchain][movermonomer+1].backbone;
		difference+=calc_backbone_bonded_difference_twoneighbors_changen(leftneighbor, rightneighbor, my_bonded_params, coordarray, coord_old, coord_new, box_dimension);
		difference+=calc_phenyl_energy_cellstruct(neighbor, leftneighbor, rightneighbor, link.number_neighbors[n], link.neighbor[n], coord_new, coordarray, box_dimension, my_nonbonded_params);
		for(i=0;i<link.number_neighbors[n];i++) difference-=link.pair_energy[n][i];
	}
	else if(coord_new.nodetype>0){							//	sidechain
		neighbor=monomerid[moverchain][movermonomer].backbone;
		difference=calc_sidechain_bonded_energy(my_bonded_params.sidechain[coord_new.nodetype], coordarray[neighbor], coord_new, box_dimension);
		difference-=calc_sidechain_bonded_energy(my_bonded_params.sidechain[coord_new.nodetype], coordarray[neighbor], coord_old, box_dimension);
		difference+=calc_phenyl_energy_cellstruct(neighbor, -1, -1, link.number_neighbors[n], link.neighbor[n], coord_new, coordarray, box_dimension, my_nonbonded_params);
		for(i=0;i<link.number_neighbors[n];i++) difference-=link.pair_energy[n][i];
	}
    else my_exit("trying to move a site with nodetype<0! (calc_energy_difference_changen_phenyl)");
	if(my_nonbonded_params.solvationparams.interface==1) difference+=calc_onebody_difference(coord_old, coord_new, chooseonebody(coord_new.nodetype, my_nonbonded_params), my_nonbonded_params.solvationparams);
	return difference;
}

double calc_energy_difference_changen_charged(int n, coord coord_old, coord coord_new, coord *coordarray, double_triple box_dimension, bonded_params my_bonded_params, nonbonded_params my_nonbonded_params, monomernodes **monomerid, int *chainlength, linkedlist link){
	int i;
	int moverchain=coord_new.chainid;
	int movermonomer=coord_new.monomerid;
	int neighbor, leftneighbor=-1, rightneighbor=-1;
	double difference=0;
	if(coord_new.nodetype==0){		//	backbone
		neighbor=monomerid[moverchain][movermonomer].sidechain;
        if(coordarray[neighbor].nodetype>=0){
            difference=calc_sidechain_bonded_energy(my_bonded_params.sidechain[coordarray[neighbor].nodetype], coord_new, coordarray[neighbor], box_dimension);
            difference-=calc_sidechain_bonded_energy(my_bonded_params.sidechain[coordarray[neighbor].nodetype], coord_old, coordarray[neighbor], box_dimension);
        }
		if(movermonomer>0) leftneighbor=monomerid[moverchain][movermonomer-1].backbone;
		if(movermonomer<chainlength[moverchain]-1) rightneighbor=monomerid[moverchain][movermonomer+1].backbone;
		difference+=calc_backbone_bonded_difference_twoneighbors_changen(leftneighbor, rightneighbor, my_bonded_params, coordarray, coord_old, coord_new, box_dimension);
		difference+=calc_charged_energy_cellstruct(neighbor, leftneighbor, rightneighbor, link.number_neighbors[n], link.neighbor[n], coord_new, coordarray, box_dimension, my_nonbonded_params);
		for(i=0;i<link.number_neighbors[n];i++) difference-=link.pair_energy[n][i];
	}
	else if(coord_new.nodetype>0){							//	sidechain
		neighbor=monomerid[moverchain][movermonomer].backbone;
        if(coord_new.nodetype>=0){
            difference=calc_sidechain_bonded_energy(my_bonded_params.sidechain[coord_new.nodetype], coordarray[neighbor], coord_new, box_dimension);
            difference-=calc_sidechain_bonded_energy(my_bonded_params.sidechain[coord_new.nodetype], coordarray[neighbor], coord_old, box_dimension);
        }
		difference+=calc_charged_energy_cellstruct(neighbor, -1, -1, link.number_neighbors[n], link.neighbor[n], coord_new, coordarray, box_dimension, my_nonbonded_params);
		for(i=0;i<link.number_neighbors[n];i++) difference-=link.pair_energy[n][i];
	}
    else my_exit("trying to move a site with nodetype<0! (calc_energy_difference_changen_charged)");
	if(my_nonbonded_params.solvationparams.interface==1) difference+=calc_onebody_difference(coord_old, coord_new, chooseonebody(coord_new.nodetype, my_nonbonded_params), my_nonbonded_params.solvationparams);
	return difference;
}

double calc_nonbonded_energy_pair(coord a, coord b, double_triple box_dimension, nonbonded_params params){
	if(a.nodetype==0){
		if(b.nodetype==0){
			return hard_energy(a.r, b.r, params.rhard0, params.rhard0, box_dimension);
		}
		else if(b.nodetype==1){
			return hard_energy(a.r, b.r, params.rhard0, params.rhard1, box_dimension);
		}
		else if(b.nodetype==2){
			return hard_energy(a.r, b.r, params.rhard0, params.p2.rhard, box_dimension);
		}
		else if(b.nodetype==3){
			return hard_energy(a.r, b.r, params.rhard0, params.p3.rhard, box_dimension);
		}
		else my_exit("don't have energy functions for these node types!");
	}
	else if(a.nodetype==1){
		if(b.nodetype==0){
			return hard_energy(a.r, b.r, params.rhard1, params.rhard0, box_dimension);
		}
		else if(b.nodetype==1){
			return phenylphenyl_energy(a, b, params.p11vac, box_dimension, params.phenylcutoff2, params.solvationparams);
            
		}
		else if(b.nodetype==2){
			return hard_energy(a.r, b.r, params.rhard1, params.p2.rhard, box_dimension);
		}
		else if(b.nodetype==3){
			return hard_energy(a.r, b.r, params.rhard1, params.p3.rhard, box_dimension);
		}
		else my_exit("don't have energy functions for these node types!");
	}
	else if(a.nodetype==2){
		if(b.nodetype==0){
			return hard_energy(a.r, b.r, params.p2.rhard, params.rhard0, box_dimension);
		}
		else if(b.nodetype==1){
			return hard_energy(a.r, b.r, params.p2.rhard, params.rhard1, box_dimension);
		}
		else if(b.nodetype==2){
			return electrostatic_energy(a, b, params.solvationparams, params.p2, params.p2, box_dimension, params.cutoff2);
		}
		else if(b.nodetype==3){
			return electrostatic_energy(a, b, params.solvationparams, params.p2, params.p3, box_dimension, params.cutoff2);
		}
		else my_exit("don't have energy functions for these node types!");
	}
	else if(a.nodetype==3){
		if(b.nodetype==0){
			return hard_energy(a.r, b.r, params.p3.rhard, params.rhard0, box_dimension);
		}
		else if(b.nodetype==1){
			return hard_energy(a.r, b.r, params.p3.rhard, params.rhard1, box_dimension);
		}
		else if(b.nodetype==2){
			return electrostatic_energy(a, b, params.solvationparams, params.p3, params.p2, box_dimension, params.cutoff2);
		}
		else if(b.nodetype==3){
			return electrostatic_energy(a, b, params.solvationparams, params.p3, params.p3, box_dimension, params.cutoff2);
		}
		else my_exit("don't have energy functions for these node types!");
	}
	else my_exit("don't have energy functions for these node types!");
	return 0;
}

double hard_energy(double_triple ar, double_triple br, double arhard, double brhard, double_triple box_dimension){
	double_triple sep=subtract_double_triple(br, ar);
	recenter_double_triple(&sep, box_dimension);
	double rmag=norm(sep);
	if(rmag>(arhard+brhard)){
		return 0;
	}
	else{
		return 1.0/0.0;
	}
}

double hard_energy_sumradii(double_triple ar, double_triple br, double sum, double_triple box_dimension){
	double_triple sep=subtract_double_triple(br, ar);
	recenter_double_triple(&sep, box_dimension);
	double rmag=norm(sep);
	if(rmag>(sum)){
		return 0;
	}
	else{
		return 1.0/0.0;
	}
}

double calc_hard_energy_cellstruct(int neighbor1, int neighbor2, int neighbor3, int numberneighbors, int *neighborlist, coord coord_new, coord *coordarray, double_triple box_dimension, nonbonded_params my_nonbonded_params){
    int i, index;
    double result=0;
    for(i=0;i<numberneighbors;i++){
        index=neighborlist[i];
        if((index!=neighbor1)&&(neighborlist[i]!=neighbor2)&&(neighborlist[i]!=neighbor3)){
            if(coord_new.nodetype==0){
                if(coordarray[index].nodetype==0){
                    result+=hard_energy_sumradii(coord_new.r, coordarray[index].r, 2.*my_nonbonded_params.rhard0, box_dimension);
                }
                else if(coordarray[index].nodetype==1){
                    result+=hard_energy_sumradii(coord_new.r, coordarray[index].r, my_nonbonded_params.rhard0+my_nonbonded_params.rhard1, box_dimension);
                }
                else if(coordarray[index].nodetype==2){
                    result+=hard_energy_sumradii(coord_new.r, coordarray[index].r, my_nonbonded_params.rhard0+my_nonbonded_params.p2.rhard, box_dimension);
                }
                else if(coordarray[index].nodetype==3){
                    result+=hard_energy_sumradii(coord_new.r, coordarray[index].r, my_nonbonded_params.rhard0+my_nonbonded_params.p3.rhard, box_dimension);
                }
            }
            else if(coord_new.nodetype==1){
                if(coordarray[index].nodetype==0){
                    result+=hard_energy_sumradii(coord_new.r, coordarray[index].r, my_nonbonded_params.rhard1+my_nonbonded_params.rhard0, box_dimension);
                }
                else if(coordarray[index].nodetype==2){
                    result+=hard_energy_sumradii(coord_new.r, coordarray[index].r, my_nonbonded_params.rhard1+my_nonbonded_params.p2.rhard, box_dimension);
                }
                else if(coordarray[index].nodetype==3){
                    result+=hard_energy_sumradii(coord_new.r, coordarray[index].r, my_nonbonded_params.rhard1+my_nonbonded_params.p3.rhard, box_dimension);
                }
            }
            else if(coord_new.nodetype==2){
                if(coordarray[index].nodetype==0){
                    result+=hard_energy_sumradii(coord_new.r, coordarray[index].r, my_nonbonded_params.p2.rhard+my_nonbonded_params.rhard0, box_dimension);
                }
                else if(coordarray[index].nodetype==1){
                    result+=hard_energy_sumradii(coord_new.r, coordarray[index].r, my_nonbonded_params.p2.rhard+my_nonbonded_params.rhard1, box_dimension);
                }
            }
            else if(coord_new.nodetype==3){
                if(coordarray[index].nodetype==0){
                    result+=hard_energy_sumradii(coord_new.r, coordarray[index].r, my_nonbonded_params.p3.rhard+my_nonbonded_params.rhard0, box_dimension);
                }
                else if(coordarray[index].nodetype==1){
                    result+=hard_energy_sumradii(coord_new.r, coordarray[index].r, my_nonbonded_params.p3.rhard+my_nonbonded_params.rhard1, box_dimension);
                }
            }
        }
    }
    return result;
}

double calc_hard_energy_cellstruct_otherpolymer(int numberneighbors, int *neighborlist, coord coord_new, coord *coordarray, double_triple box_dimension, nonbonded_params my_nonbonded_params){
	int i, index;
	double result=0;
	for(i=0;i<numberneighbors;i++){
		index=neighborlist[i];
        if(coordarray[index].chainid!=coord_new.chainid){
			if(coord_new.nodetype==0){
				if(coordarray[index].nodetype==0){
					result+=hard_energy_sumradii(coord_new.r, coordarray[index].r, 2.*my_nonbonded_params.rhard0, box_dimension);
				}
				else if(coordarray[index].nodetype==1){
					result+=hard_energy_sumradii(coord_new.r, coordarray[index].r, my_nonbonded_params.rhard0+my_nonbonded_params.rhard1, box_dimension);
				}
				else if(coordarray[index].nodetype==2){
					result+=hard_energy_sumradii(coord_new.r, coordarray[index].r, my_nonbonded_params.rhard0+my_nonbonded_params.p2.rhard, box_dimension);
				}
				else if(coordarray[index].nodetype==3){
					result+=hard_energy_sumradii(coord_new.r, coordarray[index].r, my_nonbonded_params.rhard0+my_nonbonded_params.p3.rhard, box_dimension);
				}
			}
			else if(coord_new.nodetype==1){
				if(coordarray[index].nodetype==0){
					result+=hard_energy_sumradii(coord_new.r, coordarray[index].r, my_nonbonded_params.rhard1+my_nonbonded_params.rhard0, box_dimension);
				}
				else if(coordarray[index].nodetype==2){
					result+=hard_energy_sumradii(coord_new.r, coordarray[index].r, my_nonbonded_params.rhard1+my_nonbonded_params.p2.rhard, box_dimension);
				}
				else if(coordarray[index].nodetype==3){
					result+=hard_energy_sumradii(coord_new.r, coordarray[index].r, my_nonbonded_params.rhard1+my_nonbonded_params.p3.rhard, box_dimension);
				}
			}
			else if(coord_new.nodetype==2){
				if(coordarray[index].nodetype==0){
					result+=hard_energy_sumradii(coord_new.r, coordarray[index].r, my_nonbonded_params.p2.rhard+my_nonbonded_params.rhard0, box_dimension);
				}
				else if(coordarray[index].nodetype==1){
					result+=hard_energy_sumradii(coord_new.r, coordarray[index].r, my_nonbonded_params.p2.rhard+my_nonbonded_params.rhard1, box_dimension);
				}
			}
			else if(coord_new.nodetype==3){
				if(coordarray[index].nodetype==0){
					result+=hard_energy_sumradii(coord_new.r, coordarray[index].r, my_nonbonded_params.p3.rhard+my_nonbonded_params.rhard0, box_dimension);
				}
				else if(coordarray[index].nodetype==1){
					result+=hard_energy_sumradii(coord_new.r, coordarray[index].r, my_nonbonded_params.p3.rhard+my_nonbonded_params.rhard1, box_dimension);
				}
			}
		}
	}
	return result;
}

double calc_hard_energy_cellstruct_givenpolymer(int numberneighbors, int *neighborlist, coord coord_new, coord *coordarray, double_triple box_dimension, nonbonded_params my_nonbonded_params, int polymer){
	int i, index;
	double result=0;
	for(i=0;i<numberneighbors;i++){
		index=neighborlist[i];
        if(coordarray[index].chainid==polymer){
			if(coord_new.nodetype==0){
				if(coordarray[index].nodetype==0){
					result+=hard_energy_sumradii(coord_new.r, coordarray[index].r, 2.*my_nonbonded_params.rhard0, box_dimension);
				}
				else if(coordarray[index].nodetype==1){
					result+=hard_energy_sumradii(coord_new.r, coordarray[index].r, my_nonbonded_params.rhard0+my_nonbonded_params.rhard1, box_dimension);
				}
				else if(coordarray[index].nodetype==2){
					result+=hard_energy_sumradii(coord_new.r, coordarray[index].r, my_nonbonded_params.rhard0+my_nonbonded_params.p2.rhard, box_dimension);
				}
				else if(coordarray[index].nodetype==3){
					result+=hard_energy_sumradii(coord_new.r, coordarray[index].r, my_nonbonded_params.rhard0+my_nonbonded_params.p3.rhard, box_dimension);
				}
			}
			else if(coord_new.nodetype==1){
				if(coordarray[index].nodetype==0){
					result+=hard_energy_sumradii(coord_new.r, coordarray[index].r, my_nonbonded_params.rhard1+my_nonbonded_params.rhard0, box_dimension);
				}
				else if(coordarray[index].nodetype==2){
					result+=hard_energy_sumradii(coord_new.r, coordarray[index].r, my_nonbonded_params.rhard1+my_nonbonded_params.p2.rhard, box_dimension);
				}
				else if(coordarray[index].nodetype==3){
					result+=hard_energy_sumradii(coord_new.r, coordarray[index].r, my_nonbonded_params.rhard1+my_nonbonded_params.p3.rhard, box_dimension);
				}
			}
			else if(coord_new.nodetype==2){
				if(coordarray[index].nodetype==0){
					result+=hard_energy_sumradii(coord_new.r, coordarray[index].r, my_nonbonded_params.p2.rhard+my_nonbonded_params.rhard0, box_dimension);
				}
				else if(coordarray[index].nodetype==1){
					result+=hard_energy_sumradii(coord_new.r, coordarray[index].r, my_nonbonded_params.p2.rhard+my_nonbonded_params.rhard1, box_dimension);
				}
			}
			else if(coord_new.nodetype==3){
				if(coordarray[index].nodetype==0){
					result+=hard_energy_sumradii(coord_new.r, coordarray[index].r, my_nonbonded_params.p3.rhard+my_nonbonded_params.rhard0, box_dimension);
				}
				else if(coordarray[index].nodetype==1){
					result+=hard_energy_sumradii(coord_new.r, coordarray[index].r, my_nonbonded_params.p3.rhard+my_nonbonded_params.rhard1, box_dimension);
				}
			}
		}
	}
	return result;
}

double calc_hard_energy_cellstruct_otherbilayer(int numberneighbors, int *neighborlist, coord coord_new, coord *coordarray, double_triple box_dimension, nonbonded_params my_nonbonded_params){
	int i, index, sheet=coord_new.leafid/2;
	double result=0;
	for(i=0;i<numberneighbors;i++){
		index=neighborlist[i];
        if(coordarray[index].leafid/2!=sheet){
			if(coord_new.nodetype==0){
				if(coordarray[index].nodetype==0){
					result+=hard_energy_sumradii(coord_new.r, coordarray[index].r, 2.*my_nonbonded_params.rhard0, box_dimension);
				}
				else if(coordarray[index].nodetype==1){
					result+=hard_energy_sumradii(coord_new.r, coordarray[index].r, my_nonbonded_params.rhard0+my_nonbonded_params.rhard1, box_dimension);
				}
				else if(coordarray[index].nodetype==2){
					result+=hard_energy_sumradii(coord_new.r, coordarray[index].r, my_nonbonded_params.rhard0+my_nonbonded_params.p2.rhard, box_dimension);
				}
				else if(coordarray[index].nodetype==3){
					result+=hard_energy_sumradii(coord_new.r, coordarray[index].r, my_nonbonded_params.rhard0+my_nonbonded_params.p3.rhard, box_dimension);
				}
			}
			else if(coord_new.nodetype==1){
				if(coordarray[index].nodetype==0){
					result+=hard_energy_sumradii(coord_new.r, coordarray[index].r, my_nonbonded_params.rhard1+my_nonbonded_params.rhard0, box_dimension);
				}
				else if(coordarray[index].nodetype==2){
					result+=hard_energy_sumradii(coord_new.r, coordarray[index].r, my_nonbonded_params.rhard1+my_nonbonded_params.p2.rhard, box_dimension);
				}
				else if(coordarray[index].nodetype==3){
					result+=hard_energy_sumradii(coord_new.r, coordarray[index].r, my_nonbonded_params.rhard1+my_nonbonded_params.p3.rhard, box_dimension);
				}
			}
			else if(coord_new.nodetype==2){
				if(coordarray[index].nodetype==0){
					result+=hard_energy_sumradii(coord_new.r, coordarray[index].r, my_nonbonded_params.p2.rhard+my_nonbonded_params.rhard0, box_dimension);
				}
				else if(coordarray[index].nodetype==1){
					result+=hard_energy_sumradii(coord_new.r, coordarray[index].r, my_nonbonded_params.p2.rhard+my_nonbonded_params.rhard1, box_dimension);
				}
			}
			else if(coord_new.nodetype==3){
				if(coordarray[index].nodetype==0){
					result+=hard_energy_sumradii(coord_new.r, coordarray[index].r, my_nonbonded_params.p3.rhard+my_nonbonded_params.rhard0, box_dimension);
				}
				else if(coordarray[index].nodetype==1){
					result+=hard_energy_sumradii(coord_new.r, coordarray[index].r, my_nonbonded_params.p3.rhard+my_nonbonded_params.rhard1, box_dimension);
				}
			}
		}
	}
	return result;
}

double calc_hard_energy_cellstruct_specificbilayer(int numberneighbors, int *neighborlist, coord coord_new, coord *coordarray, double_triple box_dimension, nonbonded_params my_nonbonded_params, int othersheet){
	int i, index;
	double result=0;
	for(i=0;i<numberneighbors;i++){
		index=neighborlist[i];
        if(coordarray[index].leafid/2==othersheet){
			if(coord_new.nodetype==0){
				if(coordarray[index].nodetype==0){
					result+=hard_energy_sumradii(coord_new.r, coordarray[index].r, 2.*my_nonbonded_params.rhard0, box_dimension);
				}
				else if(coordarray[index].nodetype==1){
					result+=hard_energy_sumradii(coord_new.r, coordarray[index].r, my_nonbonded_params.rhard0+my_nonbonded_params.rhard1, box_dimension);
				}
				else if(coordarray[index].nodetype==2){
					result+=hard_energy_sumradii(coord_new.r, coordarray[index].r, my_nonbonded_params.rhard0+my_nonbonded_params.p2.rhard, box_dimension);
				}
				else if(coordarray[index].nodetype==3){
					result+=hard_energy_sumradii(coord_new.r, coordarray[index].r, my_nonbonded_params.rhard0+my_nonbonded_params.p3.rhard, box_dimension);
				}
			}
			else if(coord_new.nodetype==1){
				if(coordarray[index].nodetype==0){
					result+=hard_energy_sumradii(coord_new.r, coordarray[index].r, my_nonbonded_params.rhard1+my_nonbonded_params.rhard0, box_dimension);
				}
				else if(coordarray[index].nodetype==2){
					result+=hard_energy_sumradii(coord_new.r, coordarray[index].r, my_nonbonded_params.rhard1+my_nonbonded_params.p2.rhard, box_dimension);
				}
				else if(coordarray[index].nodetype==3){
					result+=hard_energy_sumradii(coord_new.r, coordarray[index].r, my_nonbonded_params.rhard1+my_nonbonded_params.p3.rhard, box_dimension);
				}
			}
			else if(coord_new.nodetype==2){
				if(coordarray[index].nodetype==0){
					result+=hard_energy_sumradii(coord_new.r, coordarray[index].r, my_nonbonded_params.p2.rhard+my_nonbonded_params.rhard0, box_dimension);
				}
				else if(coordarray[index].nodetype==1){
					result+=hard_energy_sumradii(coord_new.r, coordarray[index].r, my_nonbonded_params.p2.rhard+my_nonbonded_params.rhard1, box_dimension);
				}
			}
			else if(coord_new.nodetype==3){
				if(coordarray[index].nodetype==0){
					result+=hard_energy_sumradii(coord_new.r, coordarray[index].r, my_nonbonded_params.p3.rhard+my_nonbonded_params.rhard0, box_dimension);
				}
				else if(coordarray[index].nodetype==1){
					result+=hard_energy_sumradii(coord_new.r, coordarray[index].r, my_nonbonded_params.p3.rhard+my_nonbonded_params.rhard1, box_dimension);
				}
			}
		}
	}
	return result;
}

double electrostatic_energy(coord a, coord b, solvation_parameters solvationparams, electrostaticparam pa, electrostaticparam pb, double_triple box_dimension, double cutoff2){
    
	// assuming a is amino
	
    double_triple r=subtract_double_triple(b.r, a.r);
	recenter_double_triple(&r, box_dimension);
	double rsq=dot_product(r, r);
	if(rsq>cutoff2) return 0;
	double rmag=sqrt(rsq);
    double rminusrhard=rmag-(pa.rhard+pb.rhard);
	if(rminusrhard<0) return 1.0/0.0;
	double eps;
	r=subtract_double_triple(add_double_triple(b.r, scalar_multiply_double_triple(b.n, pb.shift)), add_double_triple(a.r, scalar_multiply_double_triple(a.n, pa.shift)));
    
	recenter_double_triple(&r, box_dimension);
	
    rsq=dot_product(r, r);
	rmag=sqrt(rsq);
    double rmagminushard=rmag-(pa.rhard+pb.rhard-pa.shift-pb.shift);
	if(rmagminushard>solvationparams.saturationrange) eps=waterpermittivity;
	else eps=solvationparams.shortrangepermittivity+(waterpermittivity-solvationparams.shortrangepermittivity)*rmagminushard/solvationparams.saturationrange;
	double energy=pa.charge*pb.charge/rmag*(coulombconstant/eps);
    if(rminusrhard<solvationparams.waterradius){
        energy+=(pa.solvationenergy+pb.solvationenergy)*solvationparams.firstshellfraction*(1-4*pow((rminusrhard-0.5*solvationparams.waterradius)/solvationparams.waterradius, 2));
    }
	return pa.factor*pb.factor*energy;
}

int count_cg_atoms(int Nchains, int *chainlength, monomernodes **monomerid, coord *coordarray){
	int i, j, total=0, sidechainindex;
	for(i=0;i<Nchains;i++){
		for(j=0;j<chainlength[i];j++){
			total+=1;		//	backbone
			sidechainindex=monomerid[i][j].sidechain;
			if(coordarray[sidechainindex].nodetype>=0) total++;
		}
	}
	return total;
}

double calc_onebody_difference(coord coord_old, coord coord_new, onebodyparam p, solvation_parameters solvp){
	if(coord_old.r.z==coord_new.r.z) return 0;
	double affinityheight, interfaceheight;
    
    affinityheight=coord_new.height-p.z0;
    interfaceheight=coord_new.height-p.zinterface;
    
	double newaffinity=-p.solvationenergy-roomTinkcals*log(1+(exp(-p.solvationenergy/roomTinkcals)-1)/(1+exp(-4*affinityheight/solvp.interfacethickness)));
	double newinterface=-p.uinterface*exp(-0.5*pow(interfaceheight/p.sigmainterface, 2));
    
	affinityheight=coord_old.height-p.z0;
    interfaceheight=coord_old.height-p.zinterface;
    
	double oldaffinity=-p.solvationenergy-roomTinkcals*log(1+(exp(-p.solvationenergy/roomTinkcals)-1)/(1+exp(-4*affinityheight/solvp.interfacethickness)));
	double oldinterface=-p.uinterface*exp(-0.5*pow(interfaceheight/p.sigmainterface, 2));
	return (newaffinity+newinterface-oldaffinity-oldinterface);
}

void calc_onebody_energy_components(coord mycoord, onebodyparam p, solvation_parameters solvp, energycomponents *pmy_energycomponents){
	double affinityheight, interfaceheight;
    
    affinityheight=mycoord.height-p.z0;
    interfaceheight=mycoord.height-p.zinterface;
	double affinity=-p.solvationenergy-roomTinkcals*log(1+(exp(-p.solvationenergy/roomTinkcals)-1)/(1+exp(-4*affinityheight/solvp.interfacethickness)));
	double interface=-p.uinterface*exp(-0.5*pow(interfaceheight/p.sigmainterface, 2));
	(*pmy_energycomponents).solvation[mycoord.nodetype]+=affinity;
	(*pmy_energycomponents).interface[mycoord.nodetype]+=interface;
}

onebodyparam chooseonebody(int type, nonbonded_params p){
	if(type==0) return p.one0;
	if(type==1) return p.one1;
	if(type==2) return p.one2;
	if(type==3) return p.one3;
	printf("don't have onebody param for type %i!", type);
    exit(1);
    return p.one0;
}

double onebody_energy(coord mycoord, onebodyparam p, solvation_parameters solvp){
 	double affinityheight, interfaceheight;
    
    affinityheight=mycoord.height-p.z0;
    interfaceheight=mycoord.height-p.zinterface;

	double affinity=-p.solvationenergy-roomTinkcals*log(1+(exp(-p.solvationenergy/roomTinkcals)-1)/(1+exp(-4*affinityheight/solvp.interfacethickness)));
	double interface=-p.uinterface*exp(-0.5*pow(interfaceheight/p.sigmainterface, 2));
	return affinity+interface;
}

void outputcoords_centered(FILE *outp, char *name, double_triple pos, double_triple box_dimension){
	fprintf(outp, "%s\t%.6f\t%.6f\t%.6f\n", name, pos.x-0.5*box_dimension.x, pos.y-0.5*box_dimension.y, pos.z-0.5*box_dimension.z);
}

double_triple nearest_image(double_triple a, double_triple b, double_triple box){
	double_triple result=subtract_double_triple(a, b);
	recenter_double_triple(&result, box);
	return add_double_triple(b, result);
}

void updatecell(double_triple pos, linkedlistfull *plink, int n){
	int_triple new;
	new.x=(int) floor(pos.x/(*plink).cellwidth.x);
	new.y=(int) floor(pos.y/(*plink).cellwidth.y);
	new.z=(int) floor(pos.z/(*plink).cellwidth.z);
    
    //  must mod cell integer due to numerical imprecision
    
    new.x=mod(new.x, (*plink).core.cellsperside.x);
    new.y=mod(new.y, (*plink).core.cellsperside.y);
    new.z=mod(new.z, (*plink).core.cellsperside.z);
	if((new.x!=(*plink).core.cell[n].x)||(new.y!=(*plink).core.cell[n].y)||(new.z!=(*plink).core.cell[n].z)){
		if((*plink).reverselist[n]==-1){					//	n is head
			(*plink).core.head[(*plink).core.cell[n].x][(*plink).core.cell[n].y][(*plink).core.cell[n].z]=(*plink).core.list[n];
			if((*plink).core.list[n]!=-1) (*plink).reverselist[(*plink).core.list[n]]=-1;
		}
		else{
			(*plink).core.list[(*plink).reverselist[n]]=(*plink).core.list[n];
			if((*plink).core.list[n]!=-1) (*plink).reverselist[(*plink).core.list[n]]=(*plink).reverselist[n];
		}
		(*plink).core.list[n]=(*plink).core.head[new.x][new.y][new.z];
		if((*plink).core.list[n]!=-1) (*plink).reverselist[(*plink).core.list[n]]=n;
		(*plink).core.head[new.x][new.y][new.z]=n;
		(*plink).reverselist[n]=-1;
		(*plink).core.cell[n]=new;
	}
}

void updatecell_asymmetric(double_triple pos, linkedlistasymmetric *plink, int n){        //  for moving site
	int_triple new;
	new.x=(int) floor(pos.x/(*plink).cellwidth.x);
	new.y=(int) floor(pos.y/(*plink).cellwidth.y);
	new.z=(int) floor(pos.z/(*plink).cellwidth.z);
    
    //  must mod cell integer due to numerical imprecision
    
    new.x=mod(new.x, (*plink).core.cellsperside.x);
    new.y=mod(new.y, (*plink).core.cellsperside.y);
    new.z=mod(new.z, (*plink).core.cellsperside.z);
	if((new.x!=(*plink).cellreverse[n].x)||(new.y!=(*plink).cellreverse[n].y)||(new.z!=(*plink).cellreverse[n].z)){
		if((*plink).reverselist_reverse[n]==-1){					//	n is head
			(*plink).headreverse[(*plink).cellreverse[n].x][(*plink).cellreverse[n].y][(*plink).cellreverse[n].z]=(*plink).listreverse[n];
			if((*plink).listreverse[n]!=-1) (*plink).reverselist_reverse[(*plink).listreverse[n]]=-1;   // fixed
		}
		else{
			(*plink).listreverse[(*plink).reverselist_reverse[n]]=(*plink).listreverse[n];
			if((*plink).listreverse[n]!=-1) (*plink).reverselist_reverse[(*plink).listreverse[n]]=(*plink).reverselist_reverse[n];
		}
		(*plink).listreverse[n]=(*plink).headreverse[new.x][new.y][new.z];
		if((*plink).listreverse[n]!=-1) (*plink).reverselist_reverse[(*plink).listreverse[n]]=n;
		(*plink).headreverse[new.x][new.y][new.z]=n;
		(*plink).reverselist_reverse[n]=-1;
		(*plink).cellreverse[n]=new;
	}
}

void updatecell_asymmetric_reverse(double_triple pos, linkedlistasymmetric *plink, int n){        //  for moving mesh; looks just like regular updatecell
	int_triple new;
	new.x=(int) floor(pos.x/(*plink).cellwidth.x);
	new.y=(int) floor(pos.y/(*plink).cellwidth.y);
	new.z=(int) floor(pos.z/(*plink).cellwidth.z);
    
    //  must mod cell integer due to numerical imprecision
    
    new.x=mod(new.x, (*plink).core.cellsperside.x);
    new.y=mod(new.y, (*plink).core.cellsperside.y);
    new.z=mod(new.z, (*plink).core.cellsperside.z);
	if((new.x!=(*plink).core.cell[n].x)||(new.y!=(*plink).core.cell[n].y)||(new.z!=(*plink).core.cell[n].z)){
        
		if((*plink).reverselist[n]==-1){					//	n is head
			(*plink).core.head[(*plink).core.cell[n].x][(*plink).core.cell[n].y][(*plink).core.cell[n].z]=(*plink).core.list[n];
			if((*plink).core.list[n]!=-1) (*plink).reverselist[(*plink).core.list[n]]=-1;
		}
		else{
			(*plink).core.list[(*plink).reverselist[n]]=(*plink).core.list[n];
			if((*plink).core.list[n]!=-1) (*plink).reverselist[(*plink).core.list[n]]=(*plink).reverselist[n];
		}
        
		(*plink).core.list[n]=(*plink).core.head[new.x][new.y][new.z];
		if((*plink).core.list[n]!=-1) (*plink).reverselist[(*plink).core.list[n]]=n;
		(*plink).core.head[new.x][new.y][new.z]=n;
		(*plink).reverselist[n]=-1;
		(*plink).core.cell[n]=new;
	}
}

double total_energy_cellstruct(coord *coordarray, bonded_params my_bonded_params, nonbonded_params my_nonbonded_params, int Nchains, int *chainlength, monomernodes **monomerid, double_triple box_dimension, linkedlistset linkset){
	int i, j, left=-1, twoleft=-1, threeleft=-1, fourleft=-1, backboneid, sidechainid;
	double energy=0;
	for(i=0;i<Nchains;i++){
		if(chainlength[i]>2){
			for(j=0;j<chainlength[i];j++){
				backboneid=monomerid[i][j].backbone;
				if(my_nonbonded_params.solvationparams.interface==1) energy+=onebody_energy(coordarray[backboneid], chooseonebody(coordarray[backboneid].nodetype, my_nonbonded_params), my_nonbonded_params.solvationparams);
				if(j>0){
                    
					//	backbone only has short-ranged interactions

                    energy+=0.5*calc_hard_energy_cellstruct(sidechainid, twoleft, backboneid, linkset.shortrange.core.number_neighbors[left], linkset.shortrange.core.neighbor[left], coordarray[left], coordarray, box_dimension, my_nonbonded_params);
                }
				sidechainid=monomerid[i][j].sidechain;
				if(coordarray[sidechainid].nodetype>=0){
                    if(my_nonbonded_params.solvationparams.interface==1) energy+=onebody_energy(coordarray[sidechainid], chooseonebody(coordarray[sidechainid].nodetype, my_nonbonded_params), my_nonbonded_params.solvationparams);
                    energy+=0.5*calc_hard_energy_cellstruct(backboneid, -1, -1, linkset.shortrange.core.number_neighbors[sidechainid], linkset.shortrange.core.neighbor[sidechainid], coordarray[sidechainid], coordarray, box_dimension, my_nonbonded_params);
                    if(coordarray[sidechainid].nodetype==1) energy+=0.5*calc_phenyl_energy_cellstruct(backboneid, -1, -1, linkset.phenylphenyl.core.number_neighbors[sidechainid], linkset.phenylphenyl.core.neighbor[sidechainid], coordarray[sidechainid], coordarray, box_dimension, my_nonbonded_params);
                    else if(coordarray[sidechainid].nodetype>1) energy+=0.5*calc_charged_energy_cellstruct(backboneid, -1, -1, linkset.chargedcharged.core.number_neighbors[sidechainid], linkset.chargedcharged.core.neighbor[sidechainid], coordarray[sidechainid], coordarray, box_dimension, my_nonbonded_params);
                    
                    //  factor of 0.5 to account for double counting
                    
                    if(coordarray[sidechainid].nodetype>=0) energy+=calc_sidechain_bonded_energy(my_bonded_params.sidechain[coordarray[sidechainid].nodetype], coordarray[backboneid], coordarray[sidechainid], box_dimension);
                }
				if(j==2) energy+=calc_backbone_bonded_energy(my_bonded_params, coordarray[left], coordarray[twoleft], coordarray[backboneid], box_dimension);
				else if(j>=3) energy+=calc_backbone_bonded_energy_norim1i(my_bonded_params, coordarray[left], coordarray[twoleft], coordarray[backboneid], box_dimension);
				fourleft=threeleft;
				threeleft=twoleft;
				twoleft=left;
				left=backboneid;
			}
            
			//	backbone only has short-ranged interactions
            
            energy+=calc_hard_energy_cellstruct(sidechainid, twoleft, -1, linkset.shortrange.core.number_neighbors[left], linkset.shortrange.core.neighbor[left], coordarray[left], coordarray, box_dimension, my_nonbonded_params);
            
            //  single-counting unlike interactions, double counting like interactins, but it doesn't matter since all 0 or inf
		}
		else if(chainlength[i]==2){
			for(j=0;j<chainlength[i];j++){
				backboneid=monomerid[i][j].backbone;
				if(my_nonbonded_params.solvationparams.interface==1) energy+=onebody_energy(coordarray[backboneid], chooseonebody(coordarray[backboneid].nodetype, my_nonbonded_params), my_nonbonded_params.solvationparams);
				if(j==1){
					//	backbone only has short-ranged interactions
                    
                    energy+=calc_hard_energy_cellstruct(sidechainid, twoleft, backboneid, linkset.shortrange.core.number_neighbors[left], linkset.shortrange.core.neighbor[left], coordarray[left], coordarray, box_dimension, my_nonbonded_params);
				}
				sidechainid=monomerid[i][j].sidechain;
                if(coordarray[sidechainid].nodetype>=0){
                    if(my_nonbonded_params.solvationparams.interface==1) energy+=onebody_energy(coordarray[sidechainid], chooseonebody(coordarray[sidechainid].nodetype, my_nonbonded_params), my_nonbonded_params.solvationparams);
                    
                    energy+=0.5*calc_hard_energy_cellstruct(backboneid, -1, -1, linkset.shortrange.core.number_neighbors[sidechainid], linkset.shortrange.core.neighbor[sidechainid], coordarray[sidechainid], coordarray, box_dimension, my_nonbonded_params);
                    if(coordarray[sidechainid].nodetype==1) energy+=0.5*calc_phenyl_energy_cellstruct(backboneid, -1, -1, linkset.phenylphenyl.core.number_neighbors[sidechainid], linkset.phenylphenyl.core.neighbor[sidechainid], coordarray[sidechainid], coordarray, box_dimension, my_nonbonded_params);
                    else if(coordarray[sidechainid].nodetype>1) energy+=0.5*calc_charged_energy_cellstruct(backboneid, -1, -1, linkset.chargedcharged.core.number_neighbors[sidechainid], linkset.chargedcharged.core.neighbor[sidechainid], coordarray[sidechainid], coordarray, box_dimension, my_nonbonded_params);
                    
                    if(coordarray[sidechainid].nodetype>=0) energy+=calc_sidechain_bonded_energy(my_bonded_params.sidechain[coordarray[sidechainid].nodetype], coordarray[backboneid], coordarray[sidechainid], box_dimension);
                }
				if(j==1){
                    
                    //  single-counting unlike interactions, double counting like interactins, but it doesn't matter since all 0 or inf
                    
					energy+=calc_backbone_bonded_energy_onlyleft(my_bonded_params, coordarray[backboneid], coordarray[left], box_dimension);
				}
				twoleft=left;
				left=backboneid;
			}
            
			//	backbone only has short-ranged interactions
            
            energy+=calc_hard_energy_cellstruct(sidechainid, twoleft, -1, linkset.shortrange.core.number_neighbors[left], linkset.shortrange.core.neighbor[left], coordarray[left], coordarray, box_dimension, my_nonbonded_params);
            
            //  single-counting unlike interactions, double counting like interactins, but it doesn't matter since all 0 or inf
		}
        else{
			for(j=0;j<chainlength[i];j++){
				backboneid=monomerid[i][j].backbone;
				if(my_nonbonded_params.solvationparams.interface==1) energy+=onebody_energy(coordarray[backboneid], chooseonebody(coordarray[backboneid].nodetype, my_nonbonded_params), my_nonbonded_params.solvationparams);
				sidechainid=monomerid[i][j].sidechain;
                if(coordarray[sidechainid].nodetype>=0){
                    if(my_nonbonded_params.solvationparams.interface==1) energy+=onebody_energy(coordarray[sidechainid], chooseonebody(coordarray[sidechainid].nodetype, my_nonbonded_params), my_nonbonded_params.solvationparams);
                    
                    energy+=0.5*calc_hard_energy_cellstruct(backboneid, -1, -1, linkset.shortrange.core.number_neighbors[sidechainid], linkset.shortrange.core.neighbor[sidechainid], coordarray[sidechainid], coordarray, box_dimension, my_nonbonded_params);
                    if(coordarray[sidechainid].nodetype==1) energy+=0.5*calc_phenyl_energy_cellstruct(backboneid, -1, -1, linkset.phenylphenyl.core.number_neighbors[sidechainid], linkset.phenylphenyl.core.neighbor[sidechainid], coordarray[sidechainid], coordarray, box_dimension, my_nonbonded_params);
                    else if(coordarray[sidechainid].nodetype>1) energy+=0.5*calc_charged_energy_cellstruct(backboneid, -1, -1, linkset.chargedcharged.core.number_neighbors[sidechainid], linkset.chargedcharged.core.neighbor[sidechainid], coordarray[sidechainid], coordarray, box_dimension, my_nonbonded_params);
                }
				if(coordarray[sidechainid].nodetype>=0) energy+=calc_sidechain_bonded_energy(my_bonded_params.sidechain[coordarray[sidechainid].nodetype], coordarray[backboneid], coordarray[sidechainid], box_dimension);
                energy+=calc_hard_energy_cellstruct(sidechainid, -1, -1, linkset.shortrange.core.number_neighbors[left], linkset.shortrange.core.neighbor[left], coordarray[left], coordarray, box_dimension, my_nonbonded_params);
            }
        }
	}
	return energy;
}

double phenylphenyl_energy(coord a, coord b, GBQvac pvac, double_triple box_dimension, double cutoff2, solvation_parameters solvp){
	double epsprime, strength, GBeps, sigmaGB, LJGBarg, sixthpower, twelfthpower, energy=0, normr5, rhatd1sq, rhatd2sq;    
	double_triple r=subtract_double_triple(b.r, a.r);
	recenter_double_triple(&r, box_dimension);
	double normr2=dot_product(r, r);
	if(normr2>cutoff2) return 0;
	double normr=sqrt(normr2);
	double_triple rhat=scalar_multiply_double_triple(r, 1./normr);
	double rhatd1=dot_product(rhat, a.n);
	double rhatd2=dot_product(rhat, b.n);
	double d1d2=dot_product(a.n, b.n);
	double rhatd1plusrhatd2sq=(rhatd1+rhatd2)*(rhatd1+rhatd2);
	double rhatd1minusrhatd2sq=(rhatd1-rhatd2)*(rhatd1-rhatd2);
    double solvated, heighta, heightb, sigma, epsprimesolvent;
	if((solvp.interface==0)||(pvac.code==1)) solvated=1;
	else{
        if(pvac.code==0) solvated=0;
        else{            
            heighta=a.height;
            heightb=b.height;
            heighta-=pvac.p11solv.interpolatemiddle;
            heightb-=pvac.p11solv.interpolatemiddle;
            solvated=0.5*(1./(1+exp(4*heighta/pvac.p11solv.interpolatewidth))+1./(1+exp(4*heightb/pvac.p11solv.interpolatewidth)));
        }
	}
	if(solvated<1){
		epsprime=1-0.5*pvac.chiprime*(rhatd1plusrhatd2sq/(1+pvac.chiprime*d1d2)+rhatd1minusrhatd2sq/(1-pvac.chiprime*d1d2));
		strength=1./sqrt(1-pvac.chi*pvac.chi*d1d2*d1d2);
		GBeps=epsprime/(strength*strength);			//	strength^nu*epsprime^mu, fixed mu=1, nu=-2
		sigmaGB=1./sqrt(1.-0.5*pvac.chi*(rhatd1plusrhatd2sq/(1+pvac.chi*d1d2)+rhatd1minusrhatd2sq/(1-pvac.chi*d1d2)));
		LJGBarg=pvac.sigma0*pvac.xi/(normr+pvac.sigma0*(pvac.xi-sigmaGB));
		sixthpower=LJGBarg*LJGBarg*LJGBarg*LJGBarg*LJGBarg*LJGBarg;
		twelfthpower=sixthpower*sixthpower;
		energy+=(1.-solvated)*4*pvac.eps0*GBeps*(twelfthpower-sixthpower);
	}
	if(solvated>0){
		epsprime=1-0.5*pvac.p11solv.chiprime*(rhatd1plusrhatd2sq/(1+pvac.p11solv.chiprime*d1d2)+rhatd1minusrhatd2sq/(1-pvac.p11solv.chiprime*d1d2));
		strength=1./sqrt(1-pvac.p11solv.chi*pvac.p11solv.chi*d1d2*d1d2);
		GBeps=epsprime/(strength*strength);			//	strength^nu*epsprime^mu, fixed mu=1, nu=-2
		sigmaGB=1./sqrt(1.-0.5*pvac.p11solv.chi*(rhatd1plusrhatd2sq/(1+pvac.p11solv.chi*d1d2)+rhatd1minusrhatd2sq/(1-pvac.p11solv.chi*d1d2)));
		LJGBarg=pvac.p11solv.sigma0*pvac.p11solv.xi/(normr+pvac.p11solv.sigma0*(pvac.p11solv.xi-sigmaGB));
		sixthpower=LJGBarg*LJGBarg*LJGBarg*LJGBarg*LJGBarg*LJGBarg;
		twelfthpower=sixthpower*sixthpower;
		energy+=solvated*4*pvac.p11solv.eps0*GBeps*(twelfthpower-sixthpower);
	}
	if(solvated>0){
		sigma=pvac.p11solv.sigma0*sigmaGB;
		if(normr<sigma+2*pvac.p11solv.w){
			if(normr<0.5*sigma) energy+=1./0.;
			else if(normr<sigma+pvac.p11solv.w){
				epsprimesolvent=1-0.5*pvac.p11solv.chirep*(rhatd1plusrhatd2sq/(1+pvac.p11solv.chirep*d1d2)+rhatd1minusrhatd2sq/(1-pvac.p11solv.chirep*d1d2));
				energy+=solvated*pvac.p11solv.amprep*epsprimesolvent*(1-4/(pvac.p11solv.w*pvac.p11solv.w)*pow(normr-(sigma+0.5*pvac.p11solv.w), 2));
			}
			else{
				epsprimesolvent=1-0.5*pvac.p11solv.chiattr*(rhatd1plusrhatd2sq/(1+pvac.p11solv.chiattr*d1d2)+rhatd1minusrhatd2sq/(1-pvac.p11solv.chiattr*d1d2));
				energy-=solvated*pvac.p11solv.ampattr*epsprimesolvent*(1-4/(pvac.p11solv.w*pvac.p11solv.w)*pow(normr-(sigma+1.5*pvac.p11solv.w), 2));
			}
		}
	}
	if((solvated<1)&&(normr>0.5)){
		normr5=normr2*normr2*normr;
		rhatd1sq=rhatd1*rhatd1;
		rhatd2sq=rhatd2*rhatd2;
		energy+=((1.-solvated)*coulombconstant*0.75*pvac.Q*pvac.Q/normr5*(1+2*d1d2*d1d2-5*(rhatd1sq+rhatd2sq)-20*rhatd1*rhatd2*d1d2+35*rhatd1sq*rhatd2sq));
	}
    return pvac.factor*energy;
}

double phenylphenyl_energy_heightchanged(coord a, coord b, GBQvac pvac, double_triple box_dimension, double cutoff2, solvation_parameters solvp, int heightchangeda, double heightnewa, int heightchangedb, double heightnewb){
	double epsprime, strength, GBeps, sigmaGB, LJGBarg, sixthpower, twelfthpower, energy=0, normr5, rhatd1sq, rhatd2sq;
	double_triple r=subtract_double_triple(b.r, a.r);
	recenter_double_triple(&r, box_dimension);
	double normr2=dot_product(r, r);
	if(normr2>cutoff2) return 0;
	double normr=sqrt(normr2);
	double_triple rhat=scalar_multiply_double_triple(r, 1./normr);
	double rhatd1=dot_product(rhat, a.n);
	double rhatd2=dot_product(rhat, b.n);
	double d1d2=dot_product(a.n, b.n);
	double rhatd1plusrhatd2sq=(rhatd1+rhatd2)*(rhatd1+rhatd2);
	double rhatd1minusrhatd2sq=(rhatd1-rhatd2)*(rhatd1-rhatd2);
    double solvated, heighta, heightb, sigma, epsprimesolvent;
	if((solvp.interface==0)||(pvac.code==1)) solvated=1;
	else{
        if(pvac.code==0) solvated=0;
        else{
            heighta=a.height;
            heightb=b.height;
			if(heightchangeda==1) heighta=heightnewa;
			else my_exit("changed height doesn't match!");
			if(heightchangedb==1) heightb=heightnewb;
            heighta-=pvac.p11solv.interpolatemiddle;
            heightb-=pvac.p11solv.interpolatemiddle;
            solvated=0.5*(1./(1+exp(4*heighta/pvac.p11solv.interpolatewidth))+1./(1+exp(4*heightb/pvac.p11solv.interpolatewidth)));
        }
	}
	if(solvated<1){
		epsprime=1-0.5*pvac.chiprime*(rhatd1plusrhatd2sq/(1+pvac.chiprime*d1d2)+rhatd1minusrhatd2sq/(1-pvac.chiprime*d1d2));
		strength=1./sqrt(1-pvac.chi*pvac.chi*d1d2*d1d2);
		GBeps=epsprime/(strength*strength);			//	strength^nu*epsprime^mu, fixed mu=1, nu=-2
		sigmaGB=1./sqrt(1.-0.5*pvac.chi*(rhatd1plusrhatd2sq/(1+pvac.chi*d1d2)+rhatd1minusrhatd2sq/(1-pvac.chi*d1d2)));
		LJGBarg=pvac.sigma0*pvac.xi/(normr+pvac.sigma0*(pvac.xi-sigmaGB));
		sixthpower=LJGBarg*LJGBarg*LJGBarg*LJGBarg*LJGBarg*LJGBarg;
		twelfthpower=sixthpower*sixthpower;
		energy+=(1.-solvated)*4*pvac.eps0*GBeps*(twelfthpower-sixthpower);
	}
	if(solvated>0){
		epsprime=1-0.5*pvac.p11solv.chiprime*(rhatd1plusrhatd2sq/(1+pvac.p11solv.chiprime*d1d2)+rhatd1minusrhatd2sq/(1-pvac.p11solv.chiprime*d1d2));
		strength=1./sqrt(1-pvac.p11solv.chi*pvac.p11solv.chi*d1d2*d1d2);
		GBeps=epsprime/(strength*strength);			//	strength^nu*epsprime^mu, fixed mu=1, nu=-2
		sigmaGB=1./sqrt(1.-0.5*pvac.p11solv.chi*(rhatd1plusrhatd2sq/(1+pvac.p11solv.chi*d1d2)+rhatd1minusrhatd2sq/(1-pvac.p11solv.chi*d1d2)));
		LJGBarg=pvac.p11solv.sigma0*pvac.p11solv.xi/(normr+pvac.p11solv.sigma0*(pvac.p11solv.xi-sigmaGB));
		sixthpower=LJGBarg*LJGBarg*LJGBarg*LJGBarg*LJGBarg*LJGBarg;
		twelfthpower=sixthpower*sixthpower;
		energy+=solvated*4*pvac.p11solv.eps0*GBeps*(twelfthpower-sixthpower);
	}
	if(solvated>0){
		sigma=pvac.p11solv.sigma0*sigmaGB;
		if(normr<sigma+2*pvac.p11solv.w){
			if(normr<0.5*sigma) energy+=1./0.;
			else if(normr<sigma+pvac.p11solv.w){
				epsprimesolvent=1-0.5*pvac.p11solv.chirep*(rhatd1plusrhatd2sq/(1+pvac.p11solv.chirep*d1d2)+rhatd1minusrhatd2sq/(1-pvac.p11solv.chirep*d1d2));
				energy+=solvated*pvac.p11solv.amprep*epsprimesolvent*(1-4/(pvac.p11solv.w*pvac.p11solv.w)*pow(normr-(sigma+0.5*pvac.p11solv.w), 2));
			}
			else{
				epsprimesolvent=1-0.5*pvac.p11solv.chiattr*(rhatd1plusrhatd2sq/(1+pvac.p11solv.chiattr*d1d2)+rhatd1minusrhatd2sq/(1-pvac.p11solv.chiattr*d1d2));
				energy-=solvated*pvac.p11solv.ampattr*epsprimesolvent*(1-4/(pvac.p11solv.w*pvac.p11solv.w)*pow(normr-(sigma+1.5*pvac.p11solv.w), 2));
			}
		}
	}
	if((solvated<1)&&(normr>0.5)){
		normr5=normr2*normr2*normr;
		rhatd1sq=rhatd1*rhatd1;
		rhatd2sq=rhatd2*rhatd2;
		energy+=((1.-solvated)*coulombconstant*0.75*pvac.Q*pvac.Q/normr5*(1+2*d1d2*d1d2-5*(rhatd1sq+rhatd2sq)-20*rhatd1*rhatd2*d1d2+35*rhatd1sq*rhatd2sq));
	}
    return pvac.factor*energy;
}

void updateneighborlist(linkedlistfull *plink, int self, coord *coordarray, double_triple box_dimension){
	int i, j, k, i_image, j_image, k_image, neighbor;
	for(i=0;i<(*plink).core.number_neighbors[self];i++){									//	wipe out n from neighbors' neighbor lists
		neighbor=(*plink).core.neighbor[self][i];
		j=0;
		while((*plink).core.neighbor[neighbor][j]!=self) j++;
		j++;
		while(j<(*plink).core.number_neighbors[neighbor]){
			(*plink).core.neighbor[neighbor][j-1]=(*plink).core.neighbor[neighbor][j];
			j++;
		}
		(*plink).core.number_neighbors[neighbor]--;
	}
	(*plink).core.number_neighbors[self]=0;
	int_triple selfcell=(*plink).core.cell[self];
	double_triple selfpos=coordarray[self].r;
	for(i=-1;i<=1;i++){																	//	add n to new neighbors' neighbor list, looping through neighboring cells
		if((selfcell.x==0)&&(i==-1)){
			i_image=(*plink).core.cellsperside.x-1;
			selfpos.x+=box_dimension.x;
		}
		else if((selfcell.x==(*plink).core.cellsperside.x-1)&&(i==1)){
			i_image=0;
			selfpos.x-=box_dimension.x;
		}
		else i_image=selfcell.x+i;
		for(j=-1;j<=1;j++){
			if((selfcell.y==0)&&(j==-1)){
				j_image=(*plink).core.cellsperside.y-1;
				selfpos.y+=box_dimension.y;
			}
			else if((selfcell.y==(*plink).core.cellsperside.y-1)&&(j==1)){
				j_image=0;
				selfpos.y-=box_dimension.y;
			}
			else j_image=selfcell.y+j;
			for(k=-1;k<=1;k++){
				if((selfcell.z==0)&&(k==-1)){
					k_image=(*plink).core.cellsperside.z-1;
					selfpos.z+=box_dimension.z;
				}
				else if((selfcell.z==(*plink).core.cellsperside.z-1)&&(k==1)){
					k_image=0;
					selfpos.z-=box_dimension.z;
				}
				else k_image=selfcell.z+k;
				neighbor=(*plink).core.head[i_image][j_image][k_image];
				while(neighbor>=0){
					if(neighbor!=self){
						if(norm(subtract_double_triple(coordarray[neighbor].r, selfpos))<(*plink).mincellwidth){
							(*plink).core.neighbor[self][(*plink).core.number_neighbors[self]]=neighbor;
							(*plink).core.neighbor[neighbor][(*plink).core.number_neighbors[neighbor]]=self;
							(*plink).core.number_neighbors[self]++;
							(*plink).core.number_neighbors[neighbor]++;
						}
					}
					neighbor=(*plink).core.list[neighbor];
				}
				if((selfcell.z==0)&&(k==-1)){
					selfpos.z-=box_dimension.z;
				}
				else if((selfcell.z==(*plink).core.cellsperside.z-1)&&(k==1)){
					selfpos.z+=box_dimension.z;
				}
			}
			if((selfcell.y==0)&&(j==-1)){
				selfpos.y-=box_dimension.y;
			}
			else if((selfcell.y==(*plink).core.cellsperside.y-1)&&(j==1)){
				selfpos.y+=box_dimension.y;
			}
		}
		if((selfcell.x==0)&&(i==-1)){
			selfpos.x-=box_dimension.x;
		}
		else if((selfcell.x==(*plink).core.cellsperside.x-1)&&(i==1)){
			selfpos.x+=box_dimension.x;
		}
		
	}
}

void updateneighborlist_positions(linkedlistfull *plink, int self, double_triple *positions, double_triple box_dimension){
	int i, j, k, i_image, j_image, k_image, neighbor;
	for(i=0;i<(*plink).core.number_neighbors[self];i++){									//	wipe out n from neighbors' neighbor lists
		neighbor=(*plink).core.neighbor[self][i];
		j=0;
		while((*plink).core.neighbor[neighbor][j]!=self) j++;
		j++;
		while(j<(*plink).core.number_neighbors[neighbor]){
			(*plink).core.neighbor[neighbor][j-1]=(*plink).core.neighbor[neighbor][j];
			j++;
		}
		(*plink).core.number_neighbors[neighbor]--;
	}
	(*plink).core.number_neighbors[self]=0;
	int_triple selfcell=(*plink).core.cell[self];
	double_triple selfpos=positions[self];
	for(i=-1;i<=1;i++){																	//	add n to new neighbors' neighbor list, looping through neighboring cells
		if((selfcell.x==0)&&(i==-1)){
			i_image=(*plink).core.cellsperside.x-1;
			selfpos.x+=box_dimension.x;
		}
		else if((selfcell.x==(*plink).core.cellsperside.x-1)&&(i==1)){
			i_image=0;
			selfpos.x-=box_dimension.x;
		}
		else i_image=selfcell.x+i;
		for(j=-1;j<=1;j++){
			if((selfcell.y==0)&&(j==-1)){
				j_image=(*plink).core.cellsperside.y-1;
				selfpos.y+=box_dimension.y;
			}
			else if((selfcell.y==(*plink).core.cellsperside.y-1)&&(j==1)){
				j_image=0;
				selfpos.y-=box_dimension.y;
			}
			else j_image=selfcell.y+j;
			for(k=-1;k<=1;k++){
				if((selfcell.z==0)&&(k==-1)){
					k_image=(*plink).core.cellsperside.z-1;
					selfpos.z+=box_dimension.z;
				}
				else if((selfcell.z==(*plink).core.cellsperside.z-1)&&(k==1)){
					k_image=0;
					selfpos.z-=box_dimension.z;
				}
				else k_image=selfcell.z+k;
				neighbor=(*plink).core.head[i_image][j_image][k_image];
				while(neighbor>=0){
					if(neighbor!=self){
						if(norm(subtract_double_triple(positions[neighbor], selfpos))<(*plink).mincellwidth){
							(*plink).core.neighbor[self][(*plink).core.number_neighbors[self]]=neighbor;
							(*plink).core.neighbor[neighbor][(*plink).core.number_neighbors[neighbor]]=self;
							(*plink).core.number_neighbors[self]++;
							(*plink).core.number_neighbors[neighbor]++;
						}
					}
					neighbor=(*plink).core.list[neighbor];
				}
				if((selfcell.z==0)&&(k==-1)){
					selfpos.z-=box_dimension.z;
				}
				else if((selfcell.z==(*plink).core.cellsperside.z-1)&&(k==1)){
					selfpos.z+=box_dimension.z;
				}
			}
			if((selfcell.y==0)&&(j==-1)){
				selfpos.y-=box_dimension.y;
			}
			else if((selfcell.y==(*plink).core.cellsperside.y-1)&&(j==1)){
				selfpos.y+=box_dimension.y;
			}
		}
		if((selfcell.x==0)&&(i==-1)){
			selfpos.x-=box_dimension.x;
		}
		else if((selfcell.x==(*plink).core.cellsperside.x-1)&&(i==1)){
			selfpos.x+=box_dimension.x;
		}
		
	}
}

void updateneighborlist_asymmetric(linkedlistasymmetric *plink, int self, double_triple selfpos, double_triple *positions, double_triple box_dimension){        //  moving site
	int i, j, k, i_image, j_image, k_image, neighbor;
	for(i=0;i<(*plink).core.number_neighbors[self];i++){									//	wipe out self from neighbors' REVERSE neighbor lists
		neighbor=(*plink).core.neighbor[self][i];
		j=0;
		while((*plink).neighbor_reverse[neighbor][j]!=self) j++;
		j++;
		while(j<(*plink).number_neighbors_reverse[neighbor]){
			(*plink).neighbor_reverse[neighbor][j-1]=(*plink).neighbor_reverse[neighbor][j];
			j++;
		}
		(*plink).number_neighbors_reverse[neighbor]--;
	}
	(*plink).core.number_neighbors[self]=0;
	int_triple selfcell=(*plink).cellreverse[self];
	for(i=-1;i<=1;i++){																	//	add n to new neighbors' neighbor list, looping through neighboring cells
		if((selfcell.x==0)&&(i==-1)){
			i_image=(*plink).core.cellsperside.x-1;
			selfpos.x+=box_dimension.x;
		}
		else if((selfcell.x==(*plink).core.cellsperside.x-1)&&(i==1)){
			i_image=0;
			selfpos.x-=box_dimension.x;
		}
		else i_image=selfcell.x+i;
		for(j=-1;j<=1;j++){
			if((selfcell.y==0)&&(j==-1)){
				j_image=(*plink).core.cellsperside.y-1;
				selfpos.y+=box_dimension.y;
			}
			else if((selfcell.y==(*plink).core.cellsperside.y-1)&&(j==1)){
				j_image=0;
				selfpos.y-=box_dimension.y;
			}
			else j_image=selfcell.y+j;
			for(k=-1;k<=1;k++){
				if((selfcell.z==0)&&(k==-1)){
					k_image=(*plink).core.cellsperside.z-1;
					selfpos.z+=box_dimension.z;
				}
				else if((selfcell.z==(*plink).core.cellsperside.z-1)&&(k==1)){
					k_image=0;
					selfpos.z-=box_dimension.z;
				}
				else k_image=selfcell.z+k;
				neighbor=(*plink).core.head[i_image][j_image][k_image];
				while(neighbor>=0){
					if(neighbor!=self){
						if(norm(subtract_double_triple(positions[neighbor], selfpos))<(*plink).mincellwidth){         //  POSITIONS(neighbor)
							(*plink).core.neighbor[self][(*plink).core.number_neighbors[self]]=neighbor;
							(*plink).neighbor_reverse[neighbor][(*plink).number_neighbors_reverse[neighbor]]=self;     //  REVERSE NEIGHBOR LIST
							(*plink).core.number_neighbors[self]++;
							(*plink).number_neighbors_reverse[neighbor]++;    //  REVERSE NEIGHBOR LIST
						}
					}
					neighbor=(*plink).core.list[neighbor];
				}
				if((selfcell.z==0)&&(k==-1)){
					selfpos.z-=box_dimension.z;
				}
				else if((selfcell.z==(*plink).core.cellsperside.z-1)&&(k==1)){
					selfpos.z+=box_dimension.z;
				}
			}
			if((selfcell.y==0)&&(j==-1)){
				selfpos.y-=box_dimension.y;
			}
			else if((selfcell.y==(*plink).core.cellsperside.y-1)&&(j==1)){
				selfpos.y+=box_dimension.y;
			}
		}
		if((selfcell.x==0)&&(i==-1)){
			selfpos.x-=box_dimension.x;
		}
		else if((selfcell.x==(*plink).core.cellsperside.x-1)&&(i==1)){
			selfpos.x+=box_dimension.x;
		}		
	}
}

void updateneighborlist_asymmetric_reverse(linkedlistasymmetric *plink, int self, double_triple selfpos, coord *coordarray, double_triple box_dimension){        //  moving site
	int i, j, k, i_image, j_image, k_image, neighbor;
	for(i=0;i<(*plink).number_neighbors_reverse[self];i++){									//	wipe out self from neighbors' FORWARD neighbor lists
		neighbor=(*plink).neighbor_reverse[self][i];
		j=0;
		while((*plink).core.neighbor[neighbor][j]!=self) j++;
		j++;
		while(j<(*plink).core.number_neighbors[neighbor]){
			(*plink).core.neighbor[neighbor][j-1]=(*plink).core.neighbor[neighbor][j];
			j++;
		}
		(*plink).core.number_neighbors[neighbor]--;
	}
	(*plink).number_neighbors_reverse[self]=0;
	int_triple selfcell=(*plink).core.cell[self];
	for(i=-1;i<=1;i++){																	//	add n to new neighbors' neighbor list, looping through neighboring cells
		if((selfcell.x==0)&&(i==-1)){
			i_image=(*plink).core.cellsperside.x-1;
			selfpos.x+=box_dimension.x;
		}
		else if((selfcell.x==(*plink).core.cellsperside.x-1)&&(i==1)){
			i_image=0;
			selfpos.x-=box_dimension.x;
		}
		else i_image=selfcell.x+i;
		for(j=-1;j<=1;j++){
			if((selfcell.y==0)&&(j==-1)){
				j_image=(*plink).core.cellsperside.y-1;
				selfpos.y+=box_dimension.y;
			}
			else if((selfcell.y==(*plink).core.cellsperside.y-1)&&(j==1)){
				j_image=0;
				selfpos.y-=box_dimension.y;
			}
			else j_image=selfcell.y+j;
			for(k=-1;k<=1;k++){
				if((selfcell.z==0)&&(k==-1)){
					k_image=(*plink).core.cellsperside.z-1;
					selfpos.z+=box_dimension.z;
				}
				else if((selfcell.z==(*plink).core.cellsperside.z-1)&&(k==1)){
					k_image=0;
					selfpos.z-=box_dimension.z;
				}
				else k_image=selfcell.z+k;
				neighbor=(*plink).headreverse[i_image][j_image][k_image];
				while(neighbor>=0){
					if(neighbor!=self){
						if(norm(subtract_double_triple(coordarray[neighbor].r, selfpos))<(*plink).mincellwidth){         //  COORD(neighbor)
							(*plink).neighbor_reverse[self][(*plink).number_neighbors_reverse[self]]=neighbor;
							(*plink).core.neighbor[neighbor][(*plink).core.number_neighbors[neighbor]]=self;     //  FORWARD NEIGHBOR LIST
							(*plink).number_neighbors_reverse[self]++;
							(*plink).core.number_neighbors[neighbor]++;    //  FORWARD NEIGHBOR LIST
						}
					}
					neighbor=(*plink).listreverse[neighbor];
				}
				if((selfcell.z==0)&&(k==-1)){
					selfpos.z-=box_dimension.z;
				}
				else if((selfcell.z==(*plink).core.cellsperside.z-1)&&(k==1)){
					selfpos.z+=box_dimension.z;
				}
			}
			if((selfcell.y==0)&&(j==-1)){
				selfpos.y-=box_dimension.y;
			}
			else if((selfcell.y==(*plink).core.cellsperside.y-1)&&(j==1)){
				selfpos.y+=box_dimension.y;
			}
		}
		if((selfcell.x==0)&&(i==-1)){
			selfpos.x-=box_dimension.x;
		}
		else if((selfcell.x==(*plink).core.cellsperside.x-1)&&(i==1)){
			selfpos.x+=box_dimension.x;
		}
		
	}
}

void updateneighborlist_pairenergies_phenyl(linkedlistfull *plink, int self, coord *coordarray, double_triple box_dimension, nonbonded_params my_nonbonded_params){
	int i, j, k, i_image, j_image, k_image, neighbor;
	double energy;
	for(i=0;i<(*plink).core.number_neighbors[self];i++){									//	wipe out n from neighbors' neighbor lists
		neighbor=(*plink).core.neighbor[self][i];
		j=0;
		while((*plink).core.neighbor[neighbor][j]!=self) j++;
		j++;
		while(j<(*plink).core.number_neighbors[neighbor]){
			(*plink).core.neighbor[neighbor][j-1]=(*plink).core.neighbor[neighbor][j];
			(*plink).core.pair_energy[neighbor][j-1]=(*plink).core.pair_energy[neighbor][j];
			j++;
		}
		(*plink).core.number_neighbors[neighbor]--;
	}
	(*plink).core.number_neighbors[self]=0;
	int_triple selfcell=(*plink).core.cell[self];
	double_triple selfpos=coordarray[self].r;
	for(i=-1;i<=1;i++){																	//	add n to new neighbors' neighbor list, looping through neighboring cells
		if((selfcell.x==0)&&(i==-1)){
			i_image=(*plink).core.cellsperside.x-1;
			selfpos.x+=box_dimension.x;
		}
		else if((selfcell.x==(*plink).core.cellsperside.x-1)&&(i==1)){
			i_image=0;
			selfpos.x-=box_dimension.x;
		}
		else i_image=selfcell.x+i;
		for(j=-1;j<=1;j++){
			if((selfcell.y==0)&&(j==-1)){
				j_image=(*plink).core.cellsperside.y-1;
				selfpos.y+=box_dimension.y;
			}
			else if((selfcell.y==(*plink).core.cellsperside.y-1)&&(j==1)){
				j_image=0;
				selfpos.y-=box_dimension.y;
			}
			else j_image=selfcell.y+j;
			for(k=-1;k<=1;k++){
				if((selfcell.z==0)&&(k==-1)){
					k_image=(*plink).core.cellsperside.z-1;
					selfpos.z+=box_dimension.z;
				}
				else if((selfcell.z==(*plink).core.cellsperside.z-1)&&(k==1)){
					k_image=0;
					selfpos.z-=box_dimension.z;
				}
				else k_image=selfcell.z+k;
				neighbor=(*plink).core.head[i_image][j_image][k_image];
				while(neighbor>=0){
					if(neighbor!=self){
						if(norm(subtract_double_triple(coordarray[neighbor].r, selfpos))<(*plink).mincellwidth){
							(*plink).core.neighbor[self][(*plink).core.number_neighbors[self]]=neighbor;
							(*plink).core.neighbor[neighbor][(*plink).core.number_neighbors[neighbor]]=self;

							energy=phenylphenyl_energy(coordarray[self], coordarray[neighbor], my_nonbonded_params.p11vac, box_dimension, my_nonbonded_params.phenylcutoff2, my_nonbonded_params.solvationparams);
                            
							(*plink).core.pair_energy[self][(*plink).core.number_neighbors[self]]=energy;
							(*plink).core.pair_energy[neighbor][(*plink).core.number_neighbors[neighbor]]=energy;
							
							(*plink).core.number_neighbors[self]++;
							(*plink).core.number_neighbors[neighbor]++;
						}
					}
					neighbor=(*plink).core.list[neighbor];
				}
				if((selfcell.z==0)&&(k==-1)){
					selfpos.z-=box_dimension.z;
				}
				else if((selfcell.z==(*plink).core.cellsperside.z-1)&&(k==1)){
					selfpos.z+=box_dimension.z;
				}
			}
			if((selfcell.y==0)&&(j==-1)){
				selfpos.y-=box_dimension.y;
			}
			else if((selfcell.y==(*plink).core.cellsperside.y-1)&&(j==1)){
				selfpos.y+=box_dimension.y;
			}
		}
		if((selfcell.x==0)&&(i==-1)){
			selfpos.x-=box_dimension.x;
		}
		else if((selfcell.x==(*plink).core.cellsperside.x-1)&&(i==1)){
			selfpos.x+=box_dimension.x;
		}
		
	}
}

void updateneighborlist_pairenergies_charged(linkedlistfull *plink, int self, coord *coordarray, double_triple box_dimension, nonbonded_params my_nonbonded_params){
	int i, j, k, i_image, j_image, k_image, neighbor;
	double energy;
	for(i=0;i<(*plink).core.number_neighbors[self];i++){									//	wipe out n from neighbors' neighbor lists
		neighbor=(*plink).core.neighbor[self][i];
		j=0;
		while((*plink).core.neighbor[neighbor][j]!=self) j++;
		j++;
		while(j<(*plink).core.number_neighbors[neighbor]){
			(*plink).core.neighbor[neighbor][j-1]=(*plink).core.neighbor[neighbor][j];
			(*plink).core.pair_energy[neighbor][j-1]=(*plink).core.pair_energy[neighbor][j];
			j++;
		}
		(*plink).core.number_neighbors[neighbor]--;
	}
	(*plink).core.number_neighbors[self]=0;
	int_triple selfcell=(*plink).core.cell[self];
	double_triple selfpos=coordarray[self].r;
	for(i=-1;i<=1;i++){																	//	add n to new neighbors' neighbor list, looping through neighboring cells
		if((selfcell.x==0)&&(i==-1)){
			i_image=(*plink).core.cellsperside.x-1;
			selfpos.x+=box_dimension.x;
		}
		else if((selfcell.x==(*plink).core.cellsperside.x-1)&&(i==1)){
			i_image=0;
			selfpos.x-=box_dimension.x;
		}
		else i_image=selfcell.x+i;
		for(j=-1;j<=1;j++){
			if((selfcell.y==0)&&(j==-1)){
				j_image=(*plink).core.cellsperside.y-1;
				selfpos.y+=box_dimension.y;
			}
			else if((selfcell.y==(*plink).core.cellsperside.y-1)&&(j==1)){
				j_image=0;
				selfpos.y-=box_dimension.y;
			}
			else j_image=selfcell.y+j;
			for(k=-1;k<=1;k++){
				if((selfcell.z==0)&&(k==-1)){
					k_image=(*plink).core.cellsperside.z-1;
					selfpos.z+=box_dimension.z;
				}
				else if((selfcell.z==(*plink).core.cellsperside.z-1)&&(k==1)){
					k_image=0;
					selfpos.z-=box_dimension.z;
				}
				else k_image=selfcell.z+k;
				neighbor=(*plink).core.head[i_image][j_image][k_image];
				while(neighbor>=0){
					if(neighbor!=self){
						if(norm(subtract_double_triple(coordarray[neighbor].r, selfpos))<(*plink).mincellwidth){
							(*plink).core.neighbor[self][(*plink).core.number_neighbors[self]]=neighbor;
							(*plink).core.neighbor[neighbor][(*plink).core.number_neighbors[neighbor]]=self;
							
							if(coordarray[self].nodetype==2){
								if(coordarray[neighbor].nodetype==2){
									energy=electrostatic_energy(coordarray[self], coordarray[neighbor], my_nonbonded_params.solvationparams, my_nonbonded_params.p2, my_nonbonded_params.p2, box_dimension, my_nonbonded_params.cutoff2);
								}
								else if(coordarray[neighbor].nodetype==3){
									energy=electrostatic_energy(coordarray[self], coordarray[neighbor], my_nonbonded_params.solvationparams, my_nonbonded_params.p2, my_nonbonded_params.p3, box_dimension, my_nonbonded_params.cutoff2);
								}
							}
							else if(coordarray[self].nodetype==3){
								if(coordarray[neighbor].nodetype==2){
									energy=electrostatic_energy(coordarray[self], coordarray[neighbor], my_nonbonded_params.solvationparams, my_nonbonded_params.p3, my_nonbonded_params.p2, box_dimension, my_nonbonded_params.cutoff2);
								}
								else if(coordarray[neighbor].nodetype==3){
									energy=electrostatic_energy(coordarray[self], coordarray[neighbor], my_nonbonded_params.solvationparams, my_nonbonded_params.p3, my_nonbonded_params.p3, box_dimension, my_nonbonded_params.cutoff2);
								}
							}
                            
                            if(energy==1./0.){
                                printf("charged pair energy %f between %i and %i (polymers %i and %i) in updateneighborlist_pairenergies_charged!\n", energy, self, neighbor, coordarray[self].chainid, coordarray[neighbor].chainid);
                                exit(1);
                            }
                            
							(*plink).core.pair_energy[self][(*plink).core.number_neighbors[self]]=energy;
							(*plink).core.pair_energy[neighbor][(*plink).core.number_neighbors[neighbor]]=energy;
							
							(*plink).core.number_neighbors[self]++;
							(*plink).core.number_neighbors[neighbor]++;
						}
					}
					neighbor=(*plink).core.list[neighbor];
				}
				if((selfcell.z==0)&&(k==-1)){
					selfpos.z-=box_dimension.z;
				}
				else if((selfcell.z==(*plink).core.cellsperside.z-1)&&(k==1)){
					selfpos.z+=box_dimension.z;
				}
			}
			if((selfcell.y==0)&&(j==-1)){
				selfpos.y-=box_dimension.y;
			}
			else if((selfcell.y==(*plink).core.cellsperside.y-1)&&(j==1)){
				selfpos.y+=box_dimension.y;
			}
		}
		if((selfcell.x==0)&&(i==-1)){
			selfpos.x-=box_dimension.x;
		}
		else if((selfcell.x==(*plink).core.cellsperside.x-1)&&(i==1)){
			selfpos.x+=box_dimension.x;
		}		
		
	}
}

double calc_phenyl_energy_cellstruct(int neighbor1, int neighbor2, int neighbor3, int numberneighbors, int *neighborlist, coord coord_new, coord *coordarray, double_triple box_dimension, nonbonded_params my_nonbonded_params){
	int i, index;
	double result=0;
	for(i=0;i<numberneighbors;i++){
		index=neighborlist[i];
		if((index!=neighbor1)&&(neighborlist[i]!=neighbor2)&&(neighborlist[i]!=neighbor3)){
			
			// should only be in linked list if both types are phenyl
			
			result+=phenylphenyl_energy(coord_new, coordarray[index], my_nonbonded_params.p11vac, box_dimension, my_nonbonded_params.phenylcutoff2, my_nonbonded_params.solvationparams);
		}
	}
	return result;
}

double calc_charged_energy_cellstruct(int neighbor1, int neighbor2, int neighbor3, int numberneighbors, int *neighborlist, coord coord_new, coord *coordarray, double_triple box_dimension, nonbonded_params my_nonbonded_params){
	int i, index;
	double result=0;
	for(i=0;i<numberneighbors;i++){
		index=neighborlist[i];
		if((index!=neighbor1)&&(neighborlist[i]!=neighbor2)&&(neighborlist[i]!=neighbor3)){
			if(coord_new.nodetype==2){
				if(coordarray[index].nodetype==2){
					result+=electrostatic_energy(coord_new, coordarray[index], my_nonbonded_params.solvationparams, my_nonbonded_params.p2, my_nonbonded_params.p2, box_dimension, my_nonbonded_params.cutoff2);
				}
				else if(coordarray[index].nodetype==3){
					result+=electrostatic_energy(coord_new, coordarray[index], my_nonbonded_params.solvationparams, my_nonbonded_params.p2, my_nonbonded_params.p3, box_dimension, my_nonbonded_params.cutoff2);
				}
			}
			else if(coord_new.nodetype==3){
				if(coordarray[index].nodetype==2){
					result+=electrostatic_energy(coord_new, coordarray[index], my_nonbonded_params.solvationparams, my_nonbonded_params.p3, my_nonbonded_params.p2, box_dimension, my_nonbonded_params.cutoff2);
				}
				else if(coordarray[index].nodetype==3){
					result+=electrostatic_energy(coord_new, coordarray[index], my_nonbonded_params.solvationparams, my_nonbonded_params.p3, my_nonbonded_params.p3, box_dimension, my_nonbonded_params.cutoff2);
				}
			}
		}
	}
	return result;
}

double calc_phenyl_energy_cellstruct_otherpolymer(int numberneighbors, int *neighborlist, coord coord_new, coord *coordarray, double_triple box_dimension, nonbonded_params my_nonbonded_params){
	int i, index;
	double result=0;
	for(i=0;i<numberneighbors;i++){
		index=neighborlist[i];
        if(coordarray[index].chainid!=coord_new.chainid){
			
			// should only be in linked list if both types are phenyl
			
			result+=phenylphenyl_energy(coord_new, coordarray[index], my_nonbonded_params.p11vac, box_dimension, my_nonbonded_params.phenylcutoff2, my_nonbonded_params.solvationparams);
		}
	}
	return result;
}

double calc_phenyl_energy_cellstruct_givenpolymer(int numberneighbors, int *neighborlist, coord coord_new, coord *coordarray, double_triple box_dimension, nonbonded_params my_nonbonded_params, int polymer){
	int i, index;
	double result=0;
	for(i=0;i<numberneighbors;i++){
		index=neighborlist[i];
        if(coordarray[index].chainid==polymer){
			
			// should only be in linked list if both types are phenyl
			
			result+=phenylphenyl_energy(coord_new, coordarray[index], my_nonbonded_params.p11vac, box_dimension, my_nonbonded_params.phenylcutoff2, my_nonbonded_params.solvationparams);
		}
	}
	return result;
}

double calc_charged_energy_cellstruct_otherpolymer(int numberneighbors, int *neighborlist, coord coord_new, coord *coordarray, double_triple box_dimension, nonbonded_params my_nonbonded_params){
	int i, index;
	double result=0;
	for(i=0;i<numberneighbors;i++){
		index=neighborlist[i];
        if(coordarray[index].chainid!=coord_new.chainid){
			if(coord_new.nodetype==2){
				if(coordarray[index].nodetype==2){
					result+=electrostatic_energy(coord_new, coordarray[index], my_nonbonded_params.solvationparams, my_nonbonded_params.p2, my_nonbonded_params.p2, box_dimension, my_nonbonded_params.cutoff2);
				}
				else if(coordarray[index].nodetype==3){
					result+=electrostatic_energy(coord_new, coordarray[index], my_nonbonded_params.solvationparams, my_nonbonded_params.p2, my_nonbonded_params.p3, box_dimension, my_nonbonded_params.cutoff2);
				}
			}
			else if(coord_new.nodetype==3){
				if(coordarray[index].nodetype==2){
					result+=electrostatic_energy(coord_new, coordarray[index], my_nonbonded_params.solvationparams, my_nonbonded_params.p3, my_nonbonded_params.p2, box_dimension, my_nonbonded_params.cutoff2);
				}
				else if(coordarray[index].nodetype==3){
					result+=electrostatic_energy(coord_new, coordarray[index], my_nonbonded_params.solvationparams, my_nonbonded_params.p3, my_nonbonded_params.p3, box_dimension, my_nonbonded_params.cutoff2);
				}
			}
		}
	}
	return result;
}

double calc_charged_energy_cellstruct_givenpolymer(int numberneighbors, int *neighborlist, coord coord_new, coord *coordarray, double_triple box_dimension, nonbonded_params my_nonbonded_params, int polymer){
	int i, index;
	double result=0;
	for(i=0;i<numberneighbors;i++){
		index=neighborlist[i];
        if(coordarray[index].chainid==polymer){
			if(coord_new.nodetype==2){
				if(coordarray[index].nodetype==2){
					result+=electrostatic_energy(coord_new, coordarray[index], my_nonbonded_params.solvationparams, my_nonbonded_params.p2, my_nonbonded_params.p2, box_dimension, my_nonbonded_params.cutoff2);
				}
				else if(coordarray[index].nodetype==3){
					result+=electrostatic_energy(coord_new, coordarray[index], my_nonbonded_params.solvationparams, my_nonbonded_params.p2, my_nonbonded_params.p3, box_dimension, my_nonbonded_params.cutoff2);
				}
			}
			else if(coord_new.nodetype==3){
				if(coordarray[index].nodetype==2){
					result+=electrostatic_energy(coord_new, coordarray[index], my_nonbonded_params.solvationparams, my_nonbonded_params.p3, my_nonbonded_params.p2, box_dimension, my_nonbonded_params.cutoff2);
				}
				else if(coordarray[index].nodetype==3){
					result+=electrostatic_energy(coord_new, coordarray[index], my_nonbonded_params.solvationparams, my_nonbonded_params.p3, my_nonbonded_params.p3, box_dimension, my_nonbonded_params.cutoff2);
				}
			}
		}
	}
	return result;
}

double calc_phenyl_energy_cellstruct_otherbilayer(int numberneighbors, int *neighborlist, coord coord_new, coord *coordarray, double_triple box_dimension, nonbonded_params my_nonbonded_params){
	int i, index, sheet=coord_new.leafid/2;
	double result=0;
	for(i=0;i<numberneighbors;i++){
		index=neighborlist[i];
        if(coordarray[index].leafid/2!=sheet){
			
			// should only be in linked list if both types are phenyl
			
			result+=phenylphenyl_energy(coord_new, coordarray[index], my_nonbonded_params.p11vac, box_dimension, my_nonbonded_params.phenylcutoff2, my_nonbonded_params.solvationparams);
		}
	}
	return result;
}

double calc_phenyl_energy_cellstruct_specificbilayer(int numberneighbors, int *neighborlist, coord coord_new, coord *coordarray, double_triple box_dimension, nonbonded_params my_nonbonded_params, int othersheet){
	int i, index;
	double result=0;
	for(i=0;i<numberneighbors;i++){
		index=neighborlist[i];
        if(coordarray[index].leafid/2==othersheet){
			
			// should only be in linked list if both types are phenyl
			
			result+=phenylphenyl_energy(coord_new, coordarray[index], my_nonbonded_params.p11vac, box_dimension, my_nonbonded_params.phenylcutoff2, my_nonbonded_params.solvationparams);
		}
	}
	return result;
}

double calc_charged_energy_cellstruct_otherbilayer(int numberneighbors, int *neighborlist, coord coord_new, coord *coordarray, double_triple box_dimension, nonbonded_params my_nonbonded_params){
	int i, index, sheet=coord_new.leafid/2;
	double result=0;
	for(i=0;i<numberneighbors;i++){
		index=neighborlist[i];
        if(coordarray[index].leafid/2!=sheet){
			if(coord_new.nodetype==2){
				if(coordarray[index].nodetype==2){
					result+=electrostatic_energy(coord_new, coordarray[index], my_nonbonded_params.solvationparams, my_nonbonded_params.p2, my_nonbonded_params.p2, box_dimension, my_nonbonded_params.cutoff2);
				}
				else if(coordarray[index].nodetype==3){
					result+=electrostatic_energy(coord_new, coordarray[index], my_nonbonded_params.solvationparams, my_nonbonded_params.p2, my_nonbonded_params.p3, box_dimension, my_nonbonded_params.cutoff2);
				}
			}
			else if(coord_new.nodetype==3){
				if(coordarray[index].nodetype==2){
					result+=electrostatic_energy(coord_new, coordarray[index], my_nonbonded_params.solvationparams, my_nonbonded_params.p3, my_nonbonded_params.p2, box_dimension, my_nonbonded_params.cutoff2);
				}
				else if(coordarray[index].nodetype==3){
					result+=electrostatic_energy(coord_new, coordarray[index], my_nonbonded_params.solvationparams, my_nonbonded_params.p3, my_nonbonded_params.p3, box_dimension, my_nonbonded_params.cutoff2);
				}
			}
		}
	}
	return result;
}

double calc_charged_energy_cellstruct_specificbilayer(int numberneighbors, int *neighborlist, coord coord_new, coord *coordarray, double_triple box_dimension, nonbonded_params my_nonbonded_params, int othersheet){
	int i, index;
	double result=0;
	for(i=0;i<numberneighbors;i++){
		index=neighborlist[i];
        if(coordarray[index].leafid/2==othersheet){
			if(coord_new.nodetype==2){
				if(coordarray[index].nodetype==2){
					result+=electrostatic_energy(coord_new, coordarray[index], my_nonbonded_params.solvationparams, my_nonbonded_params.p2, my_nonbonded_params.p2, box_dimension, my_nonbonded_params.cutoff2);
				}
				else if(coordarray[index].nodetype==3){
					result+=electrostatic_energy(coord_new, coordarray[index], my_nonbonded_params.solvationparams, my_nonbonded_params.p2, my_nonbonded_params.p3, box_dimension, my_nonbonded_params.cutoff2);
				}
			}
			else if(coord_new.nodetype==3){
				if(coordarray[index].nodetype==2){
					result+=electrostatic_energy(coord_new, coordarray[index], my_nonbonded_params.solvationparams, my_nonbonded_params.p3, my_nonbonded_params.p2, box_dimension, my_nonbonded_params.cutoff2);
				}
				else if(coordarray[index].nodetype==3){
					result+=electrostatic_energy(coord_new, coordarray[index], my_nonbonded_params.solvationparams, my_nonbonded_params.p3, my_nonbonded_params.p3, box_dimension, my_nonbonded_params.cutoff2);
				}
			}
		}
	}
	return result;
}

double calc_energy_pair(coord a, coord b, nonbonded_params my_nonbonded_params, double_triple box_dimension){
	double result=0;
	if(a.nodetype==0){
		if(b.nodetype==0){
			result+=hard_energy_sumradii(a.r, b.r, 2.*my_nonbonded_params.rhard0, box_dimension);
		}
		else if(b.nodetype==1){
			result+=hard_energy_sumradii(a.r, b.r, my_nonbonded_params.rhard0+my_nonbonded_params.rhard1, box_dimension);
		}
		else if(b.nodetype==2){
			result+=hard_energy_sumradii(a.r, b.r, my_nonbonded_params.rhard0+my_nonbonded_params.p2.rhard, box_dimension);
		}
		else if(b.nodetype==3){
			result+=hard_energy_sumradii(a.r, b.r, my_nonbonded_params.rhard0+my_nonbonded_params.p3.rhard, box_dimension);
		}
	}
	else if(a.nodetype==1){
		if(b.nodetype==0){
			result+=hard_energy_sumradii(a.r, b.r, my_nonbonded_params.rhard1+my_nonbonded_params.rhard0, box_dimension);
		}
		else if(b.nodetype==1){
			result+=phenylphenyl_energy(a, b, my_nonbonded_params.p11vac, box_dimension, my_nonbonded_params.phenylcutoff2, my_nonbonded_params.solvationparams);
		}
		else if(b.nodetype==2){
			result+=hard_energy_sumradii(a.r, b.r, my_nonbonded_params.rhard1+my_nonbonded_params.p2.rhard, box_dimension);
		}
		else if(b.nodetype==3){
			result+=hard_energy_sumradii(a.r, b.r, my_nonbonded_params.rhard1+my_nonbonded_params.p3.rhard, box_dimension);
		}
	}
	else if(a.nodetype==2){
		if(b.nodetype==0){
			result+=hard_energy_sumradii(a.r, b.r, my_nonbonded_params.p2.rhard+my_nonbonded_params.rhard0, box_dimension);
		}
		else if(b.nodetype==1){
			result+=hard_energy_sumradii(a.r, b.r, my_nonbonded_params.p2.rhard+my_nonbonded_params.rhard1, box_dimension);
		}
		else if(b.nodetype==2){
			result+=electrostatic_energy(a, b, my_nonbonded_params.solvationparams, my_nonbonded_params.p2, my_nonbonded_params.p2, box_dimension, my_nonbonded_params.cutoff2);
		}
		else if(b.nodetype==3){
			result+=electrostatic_energy(a, b, my_nonbonded_params.solvationparams, my_nonbonded_params.p2, my_nonbonded_params.p3, box_dimension, my_nonbonded_params.cutoff2);
		}
	}
	else if(a.nodetype==3){
		if(b.nodetype==0){
			result+=hard_energy_sumradii(a.r, b.r, my_nonbonded_params.p3.rhard+my_nonbonded_params.rhard0, box_dimension);
		}
		else if(b.nodetype==1){
			result+=hard_energy_sumradii(a.r, b.r, my_nonbonded_params.p3.rhard+my_nonbonded_params.rhard1, box_dimension);
		}
		else if(b.nodetype==2){
			result+=electrostatic_energy(a, b, my_nonbonded_params.solvationparams, my_nonbonded_params.p3, my_nonbonded_params.p2, box_dimension, my_nonbonded_params.cutoff2);
		}
		else if(b.nodetype==3){
			result+=electrostatic_energy(a, b, my_nonbonded_params.solvationparams, my_nonbonded_params.p3, my_nonbonded_params.p3, box_dimension, my_nonbonded_params.cutoff2);
		}
	}
	return result;
}

double calc_backbone_bonded_difference_fourneighbors_changer(int twoleft, int left, int right, int tworight, bonded_params my_bonded_params, coord *coordarray, coord oldcoord, coord newcoord, double_triple box_dimension){
	double difference=0;
	double riip1normnew, riip1normold, riip1dotninew, riip1dotniold, riip1dotnip1new, riip1dotnip1old, rim1inormnew, rim1inormold, rim1idotninew, rim1idotniold, rim1idotnim1new, rim1idotnim1old, dotcrossproduct;
	double_triple riip1new, riip1old, rhatiip1new, rhatiip1old, rim1inew, rim1iold, rhatim1inew, rhatim1iold;
	double rip1ip2dotnip1, rip1ip2dotnip2, rim2im1dotnim1, rim2im1dotnim2;
	double_triple rip1ip2, rhatip1ip2, rim2im1, rhatim2im1;
	if(right>=0){
		riip1new=subtract_double_triple(coordarray[right].r, newcoord.r);
		riip1old=subtract_double_triple(coordarray[right].r, oldcoord.r);
		recenter_double_triple(&riip1new, box_dimension);
		recenter_double_triple(&riip1old, box_dimension);
		riip1normnew=norm(riip1new);
		riip1normold=norm(riip1old);
		rhatiip1new=scalar_multiply_double_triple(riip1new, 1./riip1normnew);
		rhatiip1old=scalar_multiply_double_triple(riip1old, 1./riip1normold);
		
        difference+=my_bonded_params.factor*(my_bonded_params.K10+riip1normnew*(my_bonded_params.K11+riip1normnew*(my_bonded_params.K12+riip1normnew*(my_bonded_params.K13+riip1normnew*my_bonded_params.K14))));	//	onlyright
        difference-=my_bonded_params.factor*(my_bonded_params.K10+riip1normold*(my_bonded_params.K11+riip1normold*(my_bonded_params.K12+riip1normold*(my_bonded_params.K13+riip1normold*my_bonded_params.K14))));	//	onlyright
		
		riip1dotninew=dot_product(rhatiip1new, newcoord.n);
		riip1dotniold=dot_product(rhatiip1old, newcoord.n);
		riip1dotnip1new=dot_product(rhatiip1new, coordarray[right].n);
		riip1dotnip1old=dot_product(rhatiip1old, coordarray[right].n);
        
        difference+=(my_bonded_params.kr*pow(riip1dotninew-(riip1normnew-my_bonded_params.r0r)/my_bonded_params.sr, 2));			
        difference-=(my_bonded_params.kr*pow(riip1dotniold-(riip1normold-my_bonded_params.r0r)/my_bonded_params.sr, 2));			
        difference+=(my_bonded_params.kl*pow(riip1dotnip1new-(riip1normnew-my_bonded_params.r0l)/my_bonded_params.sl, 2));
        difference-=(my_bonded_params.kl*pow(riip1dotnip1old-(riip1normold-my_bonded_params.r0l)/my_bonded_params.sl, 2));
        
		if(tworight>=0){
			rip1ip2=subtract_double_triple(coordarray[tworight].r, coordarray[right].r);
			recenter_double_triple(&rip1ip2, box_dimension);
			rhatip1ip2=scalar_multiply_double_triple(rip1ip2, 1./norm(rip1ip2));            
			rip1ip2dotnip1=dot_product(rhatip1ip2, coordarray[right].n);
			rip1ip2dotnip2=dot_product(rhatip1ip2, coordarray[tworight].n);
            
			dotcrossproduct=dot_product(rhatiip1new, rhatip1ip2)-riip1dotnip1new*rip1ip2dotnip1;
            difference+=(my_bonded_params.K20+dotcrossproduct*(my_bonded_params.K21+dotcrossproduct*(my_bonded_params.K22+dotcrossproduct*(my_bonded_params.K23+dotcrossproduct*my_bonded_params.K24))));
            
			dotcrossproduct=dot_product(rhatiip1old, rhatip1ip2)-riip1dotnip1old*rip1ip2dotnip1;
            difference-=(my_bonded_params.K20+dotcrossproduct*(my_bonded_params.K21+dotcrossproduct*(my_bonded_params.K22+dotcrossproduct*(my_bonded_params.K23+dotcrossproduct*my_bonded_params.K24))));
		}
	}
	if(left>=0){
		rim1inew=subtract_double_triple(newcoord.r, coordarray[left].r);
		rim1iold=subtract_double_triple(oldcoord.r, coordarray[left].r);
		recenter_double_triple(&rim1inew, box_dimension);
		recenter_double_triple(&rim1iold, box_dimension);
		rim1inormnew=norm(rim1inew);
		rim1inormold=norm(rim1iold);
		rhatim1inew=scalar_multiply_double_triple(rim1inew, 1./rim1inormnew);
		rhatim1iold=scalar_multiply_double_triple(rim1iold, 1./rim1inormold);
		
        difference+=my_bonded_params.factor*(my_bonded_params.K10+rim1inormnew*(my_bonded_params.K11+rim1inormnew*(my_bonded_params.K12+rim1inormnew*(my_bonded_params.K13+rim1inormnew*my_bonded_params.K14))));
        difference-=my_bonded_params.factor*(my_bonded_params.K10+rim1inormold*(my_bonded_params.K11+rim1inormold*(my_bonded_params.K12+rim1inormold*(my_bonded_params.K13+rim1inormold*my_bonded_params.K14))));
        
		rim1idotninew=dot_product(rhatim1inew, newcoord.n);
		rim1idotniold=dot_product(rhatim1iold, oldcoord.n);
		rim1idotnim1new=dot_product(rhatim1inew, coordarray[left].n);
		rim1idotnim1old=dot_product(rhatim1iold, coordarray[left].n);
        
        difference+=(my_bonded_params.kl*pow(rim1idotninew-(rim1inormnew-my_bonded_params.r0l)/my_bonded_params.sl, 2));
        difference-=(my_bonded_params.kl*pow(rim1idotniold-(rim1inormold-my_bonded_params.r0l)/my_bonded_params.sl, 2));
        difference+=(my_bonded_params.kr*pow(rim1idotnim1new-(rim1inormnew-my_bonded_params.r0r)/my_bonded_params.sr, 2));
        difference-=(my_bonded_params.kr*pow(rim1idotnim1old-(rim1inormold-my_bonded_params.r0r)/my_bonded_params.sr, 2));
		
		if(right>=0){
			dotcrossproduct=dot_product(rhatim1inew, rhatiip1new)-rim1idotninew*riip1dotninew;
            difference+=(my_bonded_params.K20+dotcrossproduct*(my_bonded_params.K21+dotcrossproduct*(my_bonded_params.K22+dotcrossproduct*(my_bonded_params.K23+dotcrossproduct*my_bonded_params.K24))));
            
			dotcrossproduct=dot_product(rhatim1iold, rhatiip1old)-rim1idotniold*riip1dotniold;
            difference-=(my_bonded_params.K20+dotcrossproduct*(my_bonded_params.K21+dotcrossproduct*(my_bonded_params.K22+dotcrossproduct*(my_bonded_params.K23+dotcrossproduct*my_bonded_params.K24))));
		}
		if(twoleft>=0){
			rim2im1=subtract_double_triple(coordarray[left].r, coordarray[twoleft].r);
			recenter_double_triple(&rim2im1, box_dimension);
			rhatim2im1=scalar_multiply_double_triple(rim2im1, 1./norm(rim2im1));
			
			rim2im1dotnim1=dot_product(rhatim2im1, coordarray[left].n);
			rim2im1dotnim2=dot_product(rhatim2im1, coordarray[twoleft].n);
            
			dotcrossproduct=dot_product(rhatim2im1, rhatim1inew)-rim2im1dotnim1*rim1idotnim1new;
            difference+=(my_bonded_params.K20+dotcrossproduct*(my_bonded_params.K21+dotcrossproduct*(my_bonded_params.K22+dotcrossproduct*(my_bonded_params.K23+dotcrossproduct*my_bonded_params.K24))));
            
			dotcrossproduct=dot_product(rhatim2im1, rhatim1iold)-rim2im1dotnim1*rim1idotnim1old;
            difference-=(my_bonded_params.K20+dotcrossproduct*(my_bonded_params.K21+dotcrossproduct*(my_bonded_params.K22+dotcrossproduct*(my_bonded_params.K23+dotcrossproduct*my_bonded_params.K24))));
		}
	}
	return difference;
}

double calc_backbone_bonded_difference_twoneighbors_changen(int left, int right, bonded_params my_bonded_params, coord *coordarray, coord oldcoord, coord newcoord, double_triple box_dimension){
	double difference=0;
	double riip1dotninew, riip1dotniold, rim1idotninew, rim1idotniold, dotcrossproduct, mydot, riip1norm, rim1inorm;
	double_triple riip1, rhatiip1, rim1i, rhatim1i;
	if(right>=0){
		riip1=subtract_double_triple(coordarray[right].r, newcoord.r);
		recenter_double_triple(&riip1, box_dimension);
		riip1norm=norm(riip1);
		rhatiip1=scalar_multiply_double_triple(riip1, 1./norm(riip1));
        
		riip1dotninew=dot_product(rhatiip1, newcoord.n);
		riip1dotniold=dot_product(rhatiip1, oldcoord.n);
        
        difference+=(my_bonded_params.kr*pow(riip1dotninew-(riip1norm-my_bonded_params.r0r)/my_bonded_params.sr, 2));
        difference-=(my_bonded_params.kr*pow(riip1dotniold-(riip1norm-my_bonded_params.r0r)/my_bonded_params.sr, 2));
	}
	if(left>=0){
		rim1i=subtract_double_triple(newcoord.r, coordarray[left].r);
		recenter_double_triple(&rim1i, box_dimension);
		rim1inorm=norm(rim1i);
		rhatim1i=scalar_multiply_double_triple(rim1i, 1./norm(rim1i));
        
		rim1idotninew=dot_product(rhatim1i, newcoord.n);
		rim1idotniold=dot_product(rhatim1i, oldcoord.n);
        difference+=(my_bonded_params.kl*pow(rim1idotninew-(rim1inorm-my_bonded_params.r0l)/my_bonded_params.sl, 2));
        difference-=(my_bonded_params.kl*pow(rim1idotniold-(rim1inorm-my_bonded_params.r0l)/my_bonded_params.sl, 2));
		
		if(right>=0){
			mydot=dot_product(rhatim1i, rhatiip1);
			dotcrossproduct=mydot-rim1idotninew*riip1dotninew;
            difference+=(my_bonded_params.K20+dotcrossproduct*(my_bonded_params.K21+dotcrossproduct*(my_bonded_params.K22+dotcrossproduct*(my_bonded_params.K23+dotcrossproduct*my_bonded_params.K24))));
			
            dotcrossproduct=mydot-rim1idotniold*riip1dotniold;
            difference-=(my_bonded_params.K20+dotcrossproduct*(my_bonded_params.K21+dotcrossproduct*(my_bonded_params.K22+dotcrossproduct*(my_bonded_params.K23+dotcrossproduct*my_bonded_params.K24))));
        }
	}
	return difference;
}

int count_allatom_atoms(int Nchains, int *chainlength, monomernodes **monomerid, coord *coordarray){
	int i, j, total=0, sidechainindex;
	for(i=0;i<Nchains;i++){
		for(j=0;j<chainlength[i];j++){
            if(j==0) total+=6;
            else if(j==chainlength[i]-1){
                total+=5;
            }
			else total+=7;		//	backbone
			sidechainindex=monomerid[i][j].sidechain;
            if(coordarray[sidechainindex].nodetype>=0) total+=4;           //  ethyl H's
			if(coordarray[sidechainindex].orientationtype==1){				//	perpendicular
				if(coordarray[sidechainindex].nodetype==1) total+=12;				//	phenyl
			}
			if(coordarray[sidechainindex].orientationtype==0){				//	parallel
				if(coordarray[sidechainindex].nodetype==2) total+=5;				//	amino
				if(coordarray[sidechainindex].nodetype==3) total+=4;				//	carboxyl
				if(coordarray[sidechainindex].nodetype==4) total+=1;				//	methyl
			}
		}
	}
	return total;
}

void output_xyz_allatom(int Nchains, int *chainlength, monomernodes **monomerid, int Natoms, coord *coordarray, double_triple box_dimension, char *xyzname, char *sourcename, allatomparams params, int *pframe){
	int i, j, k, backboneindex, sidechainindex, atomcount=0, lastbackbonecount, Cgammacount;
	double_triple sidevector, director, thirdvector, loc, backbonevector, com, backbonebonded, sidechainbonded, lastbackbonebonded, firstbackbonebonded, sep, lastbackbone, Cbetapos, Cgammapos, ethylvec, ethylplane, ethyloutofplane;
    
	FILE *outp, *sourceoutp;
    if((*pframe)==0){
		sourceoutp=fopen(sourcename, "w");
		fprintf(sourceoutp, "pbc set {%.6f %.6f %.6f} -first %i -last %i\n", box_dimension.x, box_dimension.y, box_dimension.z, *pframe, *pframe);
		fprintf(sourceoutp, "topo clearbonds\n");
		fprintf(sourceoutp, "set sel [atomselect top \" name B\"]\n$sel set radius %.6f\n", params.Crad);
		fprintf(sourceoutp, "set sel [atomselect top \" name G\"]\n$sel set radius %.6f\n", params.Crad);
		fprintf(sourceoutp, "set sel [atomselect top \" name P\"]\n$sel set radius %.6f\n", params.Crad);
		fprintf(sourceoutp, "set sel [atomselect top \" name A\"]\n$sel set radius %.6f\n", params.Crad);
		fprintf(sourceoutp, "set sel [atomselect top \" name L\"]\n$sel set radius %.6f\n", params.Crad);
		fprintf(sourceoutp, "set sel [atomselect top \" name D\"]\n$sel set radius %.6f\n", params.Crad);
		fprintf(sourceoutp, "set sel [atomselect top \" name C\"]\n$sel set radius %.6f\n", params.Crad);
		fprintf(sourceoutp, "set sel [atomselect top \" name O\"]\n$sel set radius %.6f\n", params.Orad);
		fprintf(sourceoutp, "set sel [atomselect top \" name R\"]\n$sel set radius %.6f\n", params.Orad);
		fprintf(sourceoutp, "set sel [atomselect top \" name N\"]\n$sel set radius %.6f\n", params.Nrad);
		fprintf(sourceoutp, "set sel [atomselect top \" name M\"]\n$sel set radius %.6f\n", params.Nrad);
        fprintf(sourceoutp, "set sel [atomselect top \" name H\"]\n$sel set radius %.6f\n", params.Hrad);
        fprintf(sourceoutp, "set sel [atomselect top \" name Q\"]\n$sel set radius %.6f\n", params.Hrad);
        fprintf(sourceoutp, "set sel [atomselect top \" name I\"]\n$sel set radius %.6f\n", params.Hrad);
        fprintf(sourceoutp, "set sel [atomselect top \" name J\"]\n$sel set radius %.6f\n", params.Hrad);
        fprintf(sourceoutp, "set sel [atomselect top \" name K\"]\n$sel set radius %.6f\n", params.Hrad);
		fprintf(sourceoutp, "color Name N white\n");      // backbone N
		fprintf(sourceoutp, "color Name B white\n");     //  C beta
		fprintf(sourceoutp, "color Name G white\n");     //  phenyl C gamma
		fprintf(sourceoutp, "color Name H white\n");    //  backbone H
		fprintf(sourceoutp, "color Name C white\n");     //  backbone C
		fprintf(sourceoutp, "color Name O white\n");      //  backbone (carbonyl) O
		fprintf(sourceoutp, "color Name P yellow\n");     //  phenyl C
		fprintf(sourceoutp, "color Name Q yellow\n");     //  phenyl H
		fprintf(sourceoutp, "color Name A white\n");     //  amino C gamma
		fprintf(sourceoutp, "color Name M blue\n");     //  amino N
		fprintf(sourceoutp, "color Name I blue\n");     //  amino H
		fprintf(sourceoutp, "color Name L white\n");     //  carboxyl C gamma
		fprintf(sourceoutp, "color Name D red\n");     //  carboxyl C
		fprintf(sourceoutp, "color Name R red\n");     //  carboxyl O
		fprintf(sourceoutp, "color Name J white\n");     //  beta H
		fprintf(sourceoutp, "color Name K white\n");     //  gamma H

		outp=fopen(xyzname, "w");
	}
	else{
		sourceoutp=fopen(sourcename, "a");
		fprintf(sourceoutp, "pbc set {%.6f %.6f %.6f} -first %i -last %i\n", box_dimension.x, box_dimension.y, box_dimension.z, *pframe, *pframe);
		fclose(sourceoutp);
		outp=fopen(xyzname, "a");
	}
	fprintf(outp, "%i\n%i\t%.6f\t%.6f\t%.6f\n", Natoms, Nchains, box_dimension.x, box_dimension.y, box_dimension.z);
	for(i=0;i<Nchains;i++){
		com.x=com.y=com.z=0;
		for(j=0;j<chainlength[i];j++){
			if(j==0){
				com=coordarray[monomerid[i][j].backbone].r;
				lastbackbone=coordarray[monomerid[i][j].backbone].r;
			}
			else{
				sep=subtract_double_triple(coordarray[monomerid[i][j].backbone].r, lastbackbone);
				recenter_double_triple(&sep, box_dimension);
				lastbackbone=add_double_triple(lastbackbone, sep);
				com=add_double_triple(com, lastbackbone);
			}
		}
		com=scalar_multiply_double_triple(com, (1./chainlength[i]));
		fmod_double_triple(&com, box_dimension);
        
        //  work left from middle to keep entire polymer in same box as center of mass
        
        for(j=chainlength[i]/2;j>=0;j--){
			backboneindex=monomerid[i][j].backbone;
			if(j==chainlength[i]/2) firstbackbonebonded=backbonebonded=nearest_image(coordarray[backboneindex].r, com, box_dimension);
			else backbonebonded=nearest_image(coordarray[backboneindex].r, lastbackbonebonded, box_dimension);
			lastbackbonebonded=backbonebonded;
        }
        
        //  now go forward, keeping polymer in same unit cell
        
        for(j=0;j<chainlength[i];j++){
			backboneindex=monomerid[i][j].backbone;
			sidechainindex=monomerid[i][j].sidechain;
            backbonebonded=nearest_image(coordarray[backboneindex].r, lastbackbonebonded, box_dimension);
            
            director=coordarray[backboneindex].n;
            if(j>0) backbonevector=subtract_double_triple(coordarray[backboneindex].r, coordarray[monomerid[i][j-1].backbone].r);
            else backbonevector=subtract_double_triple(coordarray[monomerid[i][j+1].backbone].r, coordarray[backboneindex].r);
            
            recenter_double_triple(&backbonevector, box_dimension);
            normalize(&backbonevector);
            thirdvector=cross_product(backbonevector, director);
            normalize(&thirdvector);                                                            //  this is normal to the plane where I'll put all the atoms
            sidevector=cross_product(director, thirdvector);
            normalize(&sidevector);
            
			lastbackbonebonded=backbonebonded;
            
			outputcoords_centered(outp, "N", backbonebonded, box_dimension);				//	backbone N
            if(coordarray[sidechainindex].nodetype>=0){
                loc=scalar_multiply_double_triple(director, params.Cbetaspacing);
                Cbetapos=loc=add_double_triple(backbonebonded, loc);
                outputcoords_centered(outp, "B", loc, box_dimension);                       //  Cbeta
            }
            else{
                loc=scalar_multiply_double_triple(director, params.CterminuslongHdistance);
                loc=add_double_triple(backbonebonded, loc);
                outputcoords_centered(outp, "H", loc, box_dimension);                       //  one of two terminal H's for short R terminus
            }
            if(j==0){
                loc=add_double_triple(scalar_multiply_double_triple(director, params.NterminusHdepth), scalar_multiply_double_triple(sidevector, params.NterminusHwidth));
                loc=add_double_triple(backbonebonded, loc);
                outputcoords_centered(outp, "H", loc, box_dimension);                       //  left terminal H
            }
            else{
                loc=add_double_triple(scalar_multiply_double_triple(director, params.Cnyldepth), scalar_multiply_double_triple(sidevector, params.Cnylwidth));
                loc=add_double_triple(backbonebonded, loc);
                outputcoords_centered(outp, "C", loc, box_dimension);
                loc=add_double_triple(scalar_multiply_double_triple(director, params.Onyldepth), scalar_multiply_double_triple(sidevector, params.Onylwidth));
                loc=add_double_triple(backbonebonded, loc);
                outputcoords_centered(outp, "O", loc, box_dimension);                       //  carbonyl
            }
            if(j==chainlength[i]-1){
                if(coordarray[sidechainindex].nodetype>=0){
                    loc=add_double_triple(scalar_multiply_double_triple(director, params.CterminusHdepth), scalar_multiply_double_triple(sidevector, params.CterminusHwidth));
                    loc=add_double_triple(backbonebonded, loc);
                    outputcoords_centered(outp, "H", loc, box_dimension);                       //  terminal H
                }
                else{
                    loc=add_double_triple(scalar_multiply_double_triple(director, params.CterminuslongHdepth), scalar_multiply_double_triple(sidevector, params.CterminuslongHwidth));
                    loc=add_double_triple(backbonebonded, loc);
                    outputcoords_centered(outp, "H", loc, box_dimension);                       //  terminal H
                }
            }
            else{                                                                               //  C
                loc=add_double_triple(scalar_multiply_double_triple(director, params.Cdepth), scalar_multiply_double_triple(sidevector, params.Cwidth));
                loc=add_double_triple(backbonebonded, loc);
                outputcoords_centered(outp, "C", loc, box_dimension);
                loc=add_double_triple(add_double_triple(scalar_multiply_double_triple(director, params.CHdepth), scalar_multiply_double_triple(sidevector, params.CHwidth)), scalar_multiply_double_triple(thirdvector, params.CHoutofplane));
                loc=add_double_triple(backbonebonded, loc);
                outputcoords_centered(outp, "H", loc, box_dimension);
                loc=add_double_triple(add_double_triple(scalar_multiply_double_triple(director, params.CHdepth), scalar_multiply_double_triple(sidevector, params.CHwidth)), scalar_multiply_double_triple(thirdvector, -params.CHoutofplane));
                loc=add_double_triple(backbonebonded, loc);
                outputcoords_centered(outp, "H", loc, box_dimension);
            }
            
            if((*pframe)==0){
                if(j==0){
                    fprintf(sourceoutp, "topo addbond %i %i\n", atomcount, atomcount+1);                //  N - Cbeta
                    fprintf(sourceoutp, "topo addbond %i %i\n", atomcount, atomcount+2);                //  N - terminal H
                    fprintf(sourceoutp, "topo addbond %i %i\n", atomcount, atomcount+3);                //  N - C
                    fprintf(sourceoutp, "topo addbond %i %i\n", atomcount+3, atomcount+4);              //  C - H
                    fprintf(sourceoutp, "topo addbond %i %i\n", atomcount+3, atomcount+5);              //  C - H
                }
                else if(j==1){
                    fprintf(sourceoutp, "topo addbond %i %i\n", lastbackbonecount+3, atomcount+2);  //  previous C - current Cnyl
                    fprintf(sourceoutp, "topo addbond %i %i\n", atomcount, atomcount+1);                //  N - Cbeta
                    fprintf(sourceoutp, "topo addbond %i %i\n", atomcount, atomcount+2);                //  N - Cnyl
                    fprintf(sourceoutp, "topo addbond %i %i\n", atomcount+2, atomcount+3);              //  Cnyl - O
                    fprintf(sourceoutp, "topo addbond %i %i\n", atomcount, atomcount+4);                //  N - C
                    fprintf(sourceoutp, "topo addbond %i %i\n", atomcount+4, atomcount+5);              //  C - H
                    fprintf(sourceoutp, "topo addbond %i %i\n", atomcount+4, atomcount+6);              //  C - H
                }
                else if(j==chainlength[i]-1){
                    fprintf(sourceoutp, "topo addbond %i %i\n", lastbackbonecount+4, atomcount+2);  //  previous C - current Cnyl
                    fprintf(sourceoutp, "topo addbond %i %i\n", atomcount, atomcount+1);                //  N - Cbeta
                    fprintf(sourceoutp, "topo addbond %i %i\n", atomcount, atomcount+2);                //  N - Cnyl
                    fprintf(sourceoutp, "topo addbond %i %i\n", atomcount+2, atomcount+3);              //  Cnyl - O
                    fprintf(sourceoutp, "topo addbond %i %i\n", atomcount, atomcount+4);                //  N - terminal H
                }
                else{
                    fprintf(sourceoutp, "topo addbond %i %i\n", lastbackbonecount+4, atomcount+2);  //  previous C - current Cnyl
                    fprintf(sourceoutp, "topo addbond %i %i\n", atomcount, atomcount+1);                //  N - Cbeta
                    fprintf(sourceoutp, "topo addbond %i %i\n", atomcount, atomcount+2);                //  N - Cnyl
                    fprintf(sourceoutp, "topo addbond %i %i\n", atomcount+2, atomcount+3);              //  Cnyl - O
                    fprintf(sourceoutp, "topo addbond %i %i\n", atomcount, atomcount+4);                //  N - C
                    fprintf(sourceoutp, "topo addbond %i %i\n", atomcount+4, atomcount+5);              //  C - H
                    fprintf(sourceoutp, "topo addbond %i %i\n", atomcount+4, atomcount+6);              //  C - H
                }
            }
            
            lastbackbonecount=atomcount;
            if(j==0) atomcount+=6;
            else if(j==chainlength[i]-1) atomcount+=5;
			else atomcount+=7;		//	backbone
			if(coordarray[sidechainindex].orientationtype==1){				//	perpendicular
				sidevector=subtract_double_triple(coordarray[sidechainindex].r, coordarray[backboneindex].r);
				recenter_double_triple(&sidevector, box_dimension);
				sidechainbonded=add_double_triple(backbonebonded, sidevector);
				director=coordarray[sidechainindex].n;
				thirdvector=cross_product(sidevector, director);
				normalize(&thirdvector);
				sidevector=cross_product(director, thirdvector);
				normalize(&sidevector);
				if(coordarray[sidechainindex].nodetype==1){					//	phenyl
                    Cgammacount=atomcount+6;
                    if((*pframe)==0){
                        fprintf(sourceoutp, "topo addbond %i %i\n", lastbackbonecount+1, atomcount+6);          //  beta C to gamma C
                        fprintf(sourceoutp, "topo addbond %i %i\n", atomcount+6, atomcount);            //  gamma C to aromatic ring
                        for(k=0;k<5;k++){
                            fprintf(sourceoutp, "topo addbond %i %i\n", atomcount+k, atomcount+k+1);    //  ring
                        }
                        fprintf(sourceoutp, "topo addbond %i %i\n", atomcount+5, atomcount);            //  complete ring
                        for(k=1;k<6;k++){
                            fprintf(sourceoutp, "topo addbond %i %i\n", atomcount+k, atomcount+k+6);    //  hydrogens
                        }
                    }
                    atomcount+=12;
					loc=scalar_multiply_double_triple(sidevector, -params.phenylspacing);
					loc=add_double_triple(sidechainbonded, loc);
					outputcoords_centered(outp, "P", loc, box_dimension);
					loc=add_double_triple(scalar_multiply_double_triple(sidevector, -0.5*params.phenylspacing), scalar_multiply_double_triple(thirdvector, -0.5*sqrt(3.)*params.phenylspacing));
					loc=add_double_triple(sidechainbonded, loc);
					outputcoords_centered(outp, "P", loc, box_dimension);
					loc=add_double_triple(scalar_multiply_double_triple(sidevector, 0.5*params.phenylspacing), scalar_multiply_double_triple(thirdvector, -0.5*sqrt(3.)*params.phenylspacing));
					loc=add_double_triple(sidechainbonded, loc);
					outputcoords_centered(outp, "P", loc, box_dimension);
					loc=scalar_multiply_double_triple(sidevector, params.phenylspacing);
					loc=add_double_triple(sidechainbonded, loc);
					outputcoords_centered(outp, "P", loc, box_dimension);
					loc=add_double_triple(scalar_multiply_double_triple(sidevector, 0.5*params.phenylspacing), scalar_multiply_double_triple(thirdvector, 0.5*sqrt(3.)*params.phenylspacing));
					loc=add_double_triple(sidechainbonded, loc);
					outputcoords_centered(outp, "P", loc, box_dimension);
					loc=add_double_triple(scalar_multiply_double_triple(sidevector, -0.5*params.phenylspacing), scalar_multiply_double_triple(thirdvector, 0.5*sqrt(3.)*params.phenylspacing));
					loc=add_double_triple(sidechainbonded, loc);
					outputcoords_centered(outp, "P", loc, box_dimension);
                    
                    loc=scalar_multiply_double_triple(sidevector, -params.phenylCgammaspacing);
					Cgammapos=loc=add_double_triple(sidechainbonded, loc);
					outputcoords_centered(outp, "G", loc, box_dimension);
                    
					loc=add_double_triple(scalar_multiply_double_triple(sidevector, -0.5*params.phenylHspacing), scalar_multiply_double_triple(thirdvector, -0.5*sqrt(3.)*params.phenylHspacing));
					loc=add_double_triple(sidechainbonded, loc);
					outputcoords_centered(outp, "Q", loc, box_dimension);
					loc=add_double_triple(scalar_multiply_double_triple(sidevector, 0.5*params.phenylHspacing), scalar_multiply_double_triple(thirdvector, -0.5*sqrt(3.)*params.phenylHspacing));
					loc=add_double_triple(sidechainbonded, loc);
					outputcoords_centered(outp, "Q", loc, box_dimension);
					loc=scalar_multiply_double_triple(sidevector, params.phenylHspacing);
					loc=add_double_triple(sidechainbonded, loc);
					outputcoords_centered(outp, "Q", loc, box_dimension);
					loc=add_double_triple(scalar_multiply_double_triple(sidevector, 0.5*params.phenylHspacing), scalar_multiply_double_triple(thirdvector, 0.5*sqrt(3.)*params.phenylHspacing));
					loc=add_double_triple(sidechainbonded, loc);
					outputcoords_centered(outp, "Q", loc, box_dimension);
					loc=add_double_triple(scalar_multiply_double_triple(sidevector, -0.5*params.phenylHspacing), scalar_multiply_double_triple(thirdvector, 0.5*sqrt(3.)*params.phenylHspacing));
					loc=add_double_triple(sidechainbonded, loc);
					outputcoords_centered(outp, "Q", loc, box_dimension);
                    
                }
			}
			if(coordarray[sidechainindex].orientationtype==0){		//	parallel
				sidevector=subtract_double_triple(coordarray[sidechainindex].r, coordarray[backboneindex].r);
				recenter_double_triple(&sidevector, box_dimension);
				sidechainbonded=add_double_triple(backbonebonded, sidevector);
				if(j>0){
					if(j<chainlength[i]-1) backbonevector=subtract_double_triple(coordarray[monomerid[i][j+1].backbone].r, coordarray[monomerid[i][j-1].backbone].r);
					else backbonevector=subtract_double_triple(coordarray[monomerid[i][j].backbone].r, coordarray[monomerid[i][j-1].backbone].r);
				}
				else{
					if(j<chainlength[i]-1) backbonevector=subtract_double_triple(coordarray[monomerid[i][j+1].backbone].r, coordarray[monomerid[i][j].backbone].r);
					else{
						backbonevector.x=1;	//arbitrary
						backbonevector.y=backbonevector.z=0;
					}
				}
				recenter_double_triple(&backbonevector, box_dimension);
				normalize(&backbonevector);
				director=coordarray[sidechainindex].n;
				thirdvector=cross_product(backbonevector, director);
				normalize(&thirdvector);
				sidevector=cross_product(director, thirdvector);
				normalize(&sidevector);
				if(coordarray[sidechainindex].nodetype==2){					//	amino
                    Cgammacount=atomcount+1;
                    if((*pframe)==0){
                        fprintf(sourceoutp, "topo addbond %i %i\n", lastbackbonecount+1, atomcount+1);
                        fprintf(sourceoutp, "topo addbond %i %i\n", atomcount+1, atomcount);
                        fprintf(sourceoutp, "topo addbond %i %i\n", atomcount, atomcount+2);
                        fprintf(sourceoutp, "topo addbond %i %i\n", atomcount, atomcount+3);
                        fprintf(sourceoutp, "topo addbond %i %i\n", atomcount, atomcount+4);
                    }
                    atomcount+=5;
					loc=scalar_multiply_double_triple(director, params.aminoNdepth);
					loc=add_double_triple(sidechainbonded, loc);
					outputcoords_centered(outp, "M", loc, box_dimension);
					loc=scalar_multiply_double_triple(director, params.aminoCdepth);
					Cgammapos=loc=add_double_triple(sidechainbonded, loc);
					outputcoords_centered(outp, "A", loc, box_dimension);
                    
					loc=add_double_triple(scalar_multiply_double_triple(director, params.aminoHdepth), scalar_multiply_double_triple(sidevector, params.aminoHwidth));
					loc=add_double_triple(sidechainbonded, loc);
					outputcoords_centered(outp, "I", loc, box_dimension);
					loc=add_double_triple(add_double_triple(scalar_multiply_double_triple(director, params.aminoHdepth), scalar_multiply_double_triple(sidevector, -0.5*params.aminoHwidth)), scalar_multiply_double_triple(thirdvector, 0.5*sqrt(3.)*params.aminoHwidth));
					loc=add_double_triple(sidechainbonded, loc);
					outputcoords_centered(outp, "I", loc, box_dimension);
					loc=add_double_triple(add_double_triple(scalar_multiply_double_triple(director, params.aminoHdepth), scalar_multiply_double_triple(sidevector, -0.5*params.aminoHwidth)), scalar_multiply_double_triple(thirdvector, -0.5*sqrt(3.)*params.aminoHwidth));
					loc=add_double_triple(sidechainbonded, loc);
					outputcoords_centered(outp, "I", loc, box_dimension);
				}
				if(coordarray[sidechainindex].nodetype==3){					//	carboxyl
                    Cgammacount=atomcount;
                    if((*pframe)==0){
                        fprintf(sourceoutp, "topo addbond %i %i\n", lastbackbonecount+1, atomcount);
                        fprintf(sourceoutp, "topo addbond %i %i\n", atomcount, atomcount+1);
                        fprintf(sourceoutp, "topo addbond %i %i\n", atomcount+1, atomcount+2);
                        fprintf(sourceoutp, "topo addbond %i %i\n", atomcount+1, atomcount+3);
                    }
                    atomcount+=4;
					loc=scalar_multiply_double_triple(director, params.carboxylbackCdepth);
					Cgammapos=loc=add_double_triple(sidechainbonded, loc);
					outputcoords_centered(outp, "L", loc, box_dimension);
					loc=scalar_multiply_double_triple(director, params.carboxylforwardCdepth);
					loc=add_double_triple(sidechainbonded, loc);
					outputcoords_centered(outp, "D", loc, box_dimension);
					loc=add_double_triple(scalar_multiply_double_triple(director, params.carboxylOdepth), scalar_multiply_double_triple(sidevector, -params.carboxylOwidth));
					loc=add_double_triple(sidechainbonded, loc);
					outputcoords_centered(outp, "R", loc, box_dimension);
					loc=add_double_triple(scalar_multiply_double_triple(director, params.carboxylOdepth), scalar_multiply_double_triple(sidevector, params.carboxylOwidth));
					loc=add_double_triple(sidechainbonded, loc);
					outputcoords_centered(outp, "R", loc, box_dimension);
				}
			}
            if(coordarray[sidechainindex].nodetype>=0){
                
                //  C beta ethyl H's
                
                ethylvec=normed(add_double_triple(normed(subtract_double_triple(Cbetapos, backbonebonded)), normed(subtract_double_triple(Cbetapos, Cgammapos))));
                ethylplane=subtract_double_triple(Cgammapos, backbonebonded);
                ethyloutofplane=normed(cross_product(ethylvec, ethylplane));
                loc=add_double_triple(Cbetapos, add_double_triple(scalar_multiply_double_triple(ethylvec, params.ethylHdepth), scalar_multiply_double_triple(ethyloutofplane, params.ethylHoutofplane)));
                outputcoords_centered(outp, "J", loc, box_dimension);
                loc=add_double_triple(Cbetapos, add_double_triple(scalar_multiply_double_triple(ethylvec, params.ethylHdepth), scalar_multiply_double_triple(ethyloutofplane, -params.ethylHoutofplane)));
                outputcoords_centered(outp, "J", loc, box_dimension);
                if((*pframe)==0){
                    fprintf(sourceoutp, "topo addbond %i %i\n", lastbackbonecount+1, atomcount);
                    fprintf(sourceoutp, "topo addbond %i %i\n", lastbackbonecount+1, atomcount+1);
                }
                
                //  C gamma ethyl H's
                
                ethylvec=normed(add_double_triple(normed(subtract_double_triple(Cgammapos, Cbetapos)), normed(subtract_double_triple(Cgammapos, sidechainbonded))));
                ethylplane=subtract_double_triple(sidechainbonded, Cbetapos);
                ethyloutofplane=normed(cross_product(ethylvec, ethylplane));
                loc=add_double_triple(Cgammapos, add_double_triple(scalar_multiply_double_triple(ethylvec, params.ethylHdepth), scalar_multiply_double_triple(ethyloutofplane, params.ethylHoutofplane)));
                outputcoords_centered(outp, "K", loc, box_dimension);
                loc=add_double_triple(Cgammapos, add_double_triple(scalar_multiply_double_triple(ethylvec, params.ethylHdepth), scalar_multiply_double_triple(ethyloutofplane, -params.ethylHoutofplane)));
                outputcoords_centered(outp, "K", loc, box_dimension);
                if((*pframe)==0){
                    fprintf(sourceoutp, "topo addbond %i %i\n", Cgammacount, atomcount+2);
                    fprintf(sourceoutp, "topo addbond %i %i\n", Cgammacount, atomcount+3);
                }
                atomcount+=4;
            }
		}
	}
    if((*pframe)==0) fclose(sourceoutp);
	fclose(outp);
    (*pframe)++;
}

void calc_com_trajectory(int Nnodes, coord *coordarray, double_triple box_dimension, double_triple *pavcom, int *pcomcount){
    int i, realnodes=0;
    double_triple com, sep, pos;
    for(i=0;i<Nnodes;i++){
        if(coordarray[i].nodetype>=0){
            realnodes++;
            if(i==0){
                pos=coordarray[i].r;
                com=pos;
            }
            else{
                sep=subtract_double_triple(coordarray[i].r, coordarray[i-1].r);
                recenter_double_triple(&sep, box_dimension);
                pos=add_double_triple(pos, sep);
                com=add_double_triple(com, pos);
            }
        }
	}
	com=scalar_multiply_double_triple(com, (1./realnodes));
    fmod_double_triple(&com, box_dimension);
    (*pavcom)=add_double_triple((*pavcom), com);
    (*pcomcount)++;
}

void output_com(char *filename, double_triple *pavcom, int *pcomcount, int cycle, int interface, double *avmeshheight, int meshnumber){
    FILE *outp;
    outp=fopen(filename, "a");
    int i;
    fprintf(outp, "%i %f %f %f", cycle, (*pavcom).x/(*pcomcount), (*pavcom).y/(*pcomcount), (*pavcom).z/(*pcomcount));
    if(interface==1){
        for(i=0;i<2;i++){
            fprintf(outp, " %f", avmeshheight[i]/(1.*meshnumber/2*(*pcomcount)));
            avmeshheight[i]=0;
        }
    }
    fprintf(outp, "\n");
    (*pavcom).x=(*pavcom).y=(*pavcom).z=0;
    (*pcomcount)=0;
    fclose(outp);
}

void expose_input_polymer_to_interface(char *filename, int *pNnodes, int *pNchains, int **chainlength, monomernodes ***monomerid, coord **coordarray, double_triple *pbox_dimension, int *pinterface, double *pvacuumthickness, int *pNnodetypes, monomertypes *pmonomercount, long long int *pt, int *pframe, int *prunningtime, double newheight, double newwidth, double newdepth, double newvacuumthickness){
	int i, j, inputinterface;
	FILE *inp;
	inp=fopen(filename, "r");
	fscanf(inp, "%lf %lf %lf", &((*pbox_dimension).x), &((*pbox_dimension).y), &((*pbox_dimension).z));
	fscanf(inp, "%i", &(inputinterface));
	fscanf(inp, "%i %i %i %i %i %i", &(*pNnodes), &(*pNchains), &((*pmonomercount).charged), &((*pmonomercount).nonpolar), &((*pmonomercount).monomers), &(*pNnodetypes));
	if((*pNchains)!=1){
		printf("Nchains=%i in expose_input_polymer_to_interface!\n", *pNchains);
		exit(1);
	}
	for(i=0;i<(*pNnodetypes);i++) fscanf(inp, " %i", &((*pmonomercount).type[i]));
	(*coordarray)=xcalloc((*pNnodes), sizeof(coord));
	(*chainlength)=xcalloc((*pNchains), sizeof(int));
	(*monomerid)=xcalloc((*pNchains), sizeof(monomernodes *));
	for(i=0;i<*pNchains;i++){
		fscanf(inp, "%i", &((*chainlength)[i]));
		((*monomerid)[i])=xcalloc(((*chainlength)[i]), sizeof(monomernodes));
		for(j=0;j<(*chainlength)[i];j++){
			fscanf(inp, "%i %i", &((*monomerid)[i][j].backbone), &((*monomerid)[i][j].sidechain));
		}
	}
	double_triple com, sep, newpos;
	int realnodes=0;
	for(i=0;i<(*pNnodes);i++){
		fscanf(inp, "%lf %lf %lf %lf %lf %lf %i %i %i %i %i", &((*coordarray)[i].r.x), &((*coordarray)[i].r.y), &((*coordarray)[i].r.z), &((*coordarray)[i].n.x), &((*coordarray)[i].n.y), &((*coordarray)[i].n.z), &((*coordarray)[i].nodetype), &((*coordarray)[i].chainid), &((*coordarray)[i].monomerid), &((*coordarray)[i].orientationtype), &((*coordarray)[i].leafid));
        if((*coordarray)[i].nodetype>=0){
            realnodes++;
            if(i==0){
                com=(*coordarray)[i].r;
            }
            else{
                sep=subtract_double_triple((*coordarray)[i].r, (*coordarray)[i-1].r);
                recenter_double_triple(&sep, *pbox_dimension);
                (*coordarray)[i].r=add_double_triple((*coordarray)[i-1].r, sep);
                com=add_double_triple(com, (*coordarray)[i].r);
            }
        }
	}
    fscanf(inp, "%lli %i %i", &(*pt), &(*pframe), &(*prunningtime));
	com=scalar_multiply_double_triple(com, (1./realnodes));
	
	for(i=0;i<(*pNnodes);i++){
        if((*coordarray)[i].nodetype>=0){			
            (*coordarray)[i].r=subtract_double_triple((*coordarray)[i].r, com);
			recenter_double_triple(&((*coordarray)[i].r), (*pbox_dimension));				//	putting coordinates near 0 before changing box_dimension
        }
	}
    
	(*pbox_dimension).x=(*pbox_dimension).y=newwidth;
	(*pvacuumthickness)=newvacuumthickness;
	(*pbox_dimension).z=newheight+(*pvacuumthickness);
	newpos.x=0.5*(*pbox_dimension).x;
	newpos.y=0.5*(*pbox_dimension).y;
	newpos.z=(*pbox_dimension).z-0.5*(*pvacuumthickness)-newdepth;	
	for(i=0;i<(*pNnodes);i++){
        if((*coordarray)[i].nodetype>=0){
			(*coordarray)[i].r=add_double_triple((*coordarray)[i].r, newpos);
        }
	}
    
}

void initialize_interface_mesh(solvation_parameters *psolvationparameters, double_triple box_dimension, double vacuumthickness, double_triple **pmeshpositions, triangledata **pmeshtriangledata, bonddata **pmeshbonddata, vertexdata **pmeshdata){
	int i, j;
	(*psolvationparameters).meshsize.x=(int) floor(box_dimension.x/(*psolvationparameters).meshbondlength);
	(*psolvationparameters).meshsize.y=2*((int) floor(0.5*box_dimension.y/(0.5*sqrt(3.)*(*psolvationparameters).meshbondlength)));      //  divisible by 2
    (*psolvationparameters).meshnumber=2*(*psolvationparameters).meshsize.x*(*psolvationparameters).meshsize.y;
    
	allocate_array(double_triple, (*pmeshpositions), (*psolvationparameters).meshnumber);
	allocate_array(vertexdata, (*pmeshdata), (*psolvationparameters).meshnumber);
	allocate_array(triangledata, (*pmeshtriangledata), 2*(*psolvationparameters).meshnumber);
	allocate_array(bonddata, (*pmeshbonddata), 3*(*psolvationparameters).meshnumber);
	double_triple spacing;
	spacing.x=box_dimension.x/(*psolvationparameters).meshsize.x;
	spacing.y=box_dimension.y/(*psolvationparameters).meshsize.y;
	for(i=0;i<(*psolvationparameters).meshsize.y;i++){
		for(j=0;j<(*psolvationparameters).meshsize.x;j++){
			(*pmeshpositions)[i*(*psolvationparameters).meshsize.x+j].y=i*spacing.y;
			(*pmeshpositions)[i*(*psolvationparameters).meshsize.x+j].x=((mod(i, 2))*0.5+j)*spacing.x;
			(*pmeshpositions)[i*(*psolvationparameters).meshsize.x+j].z=0.5*vacuumthickness;
            (*pmeshtriangledata)[2*(i*(*psolvationparameters).meshsize.x+j)].area=0.5*spacing.x*spacing.y;
            (*pmeshtriangledata)[2*(i*(*psolvationparameters).meshsize.x+j)+1].area=0.5*spacing.x*spacing.y;
            (*pmeshtriangledata)[2*(i*(*psolvationparameters).meshsize.x+j)].columnvolume=0.5*spacing.x*spacing.y*(0.5*vacuumthickness);
            (*pmeshtriangledata)[2*(i*(*psolvationparameters).meshsize.x+j)+1].columnvolume=0.5*spacing.x*spacing.y*(0.5*vacuumthickness);
            (*pmeshtriangledata)[2*(i*(*psolvationparameters).meshsize.x+j)].normal.x=0;
            (*pmeshtriangledata)[2*(i*(*psolvationparameters).meshsize.x+j)].normal.y=0;
            (*pmeshtriangledata)[2*(i*(*psolvationparameters).meshsize.x+j)].normal.z=1;
            (*pmeshtriangledata)[2*(i*(*psolvationparameters).meshsize.x+j)+1].normal.x=0;
            (*pmeshtriangledata)[2*(i*(*psolvationparameters).meshsize.x+j)+1].normal.y=0;
            (*pmeshtriangledata)[2*(i*(*psolvationparameters).meshsize.x+j)+1].normal.z=1;
            (*pmeshbonddata)[3*(i*(*psolvationparameters).meshsize.x+j)].length=spacing.x;
            (*pmeshbonddata)[3*(i*(*psolvationparameters).meshsize.x+j)+1].length=sqrt(pow(0.5*spacing.x, 2)+pow(spacing.y, 2));
            (*pmeshbonddata)[3*(i*(*psolvationparameters).meshsize.x+j)+2].length=sqrt(pow(0.5*spacing.x, 2)+pow(spacing.y, 2));
            (*pmeshbonddata)[3*(i*(*psolvationparameters).meshsize.x+j)].meancurvature=0;
            (*pmeshbonddata)[3*(i*(*psolvationparameters).meshsize.x+j)+1].meancurvature=0;
            (*pmeshbonddata)[3*(i*(*psolvationparameters).meshsize.x+j)+2].meancurvature=0;
			(*pmeshdata)[i*(*psolvationparameters).meshsize.x+j].meancurvature=0;
			(*pmeshdata)[i*(*psolvationparameters).meshsize.x+j].meancurvaturesquared=0;
			(*pmeshdata)[i*(*psolvationparameters).meshsize.x+j].area=3*spacing.x*spacing.y;
		}
	}
	for(i=0;i<(*psolvationparameters).meshsize.y;i++){
		for(j=0;j<(*psolvationparameters).meshsize.x;j++){
			(*pmeshpositions)[((*psolvationparameters).meshsize.y+i)*(*psolvationparameters).meshsize.x+j].y=i*spacing.y;
			(*pmeshpositions)[((*psolvationparameters).meshsize.y+i)*(*psolvationparameters).meshsize.x+j].x=((mod(i, 2))*0.5+j)*spacing.x;
			(*pmeshpositions)[((*psolvationparameters).meshsize.y+i)*(*psolvationparameters).meshsize.x+j].z=box_dimension.z-0.5*vacuumthickness;
            (*pmeshtriangledata)[2*(((*psolvationparameters).meshsize.y+i)*(*psolvationparameters).meshsize.x+j)].area=0.5*spacing.x*spacing.y;
            (*pmeshtriangledata)[2*(((*psolvationparameters).meshsize.y+i)*(*psolvationparameters).meshsize.x+j)+1].area=0.5*spacing.x*spacing.y;
            (*pmeshtriangledata)[2*(((*psolvationparameters).meshsize.y+i)*(*psolvationparameters).meshsize.x+j)].columnvolume=0.5*spacing.x*spacing.y*(box_dimension.z-0.5*vacuumthickness);
            (*pmeshtriangledata)[2*(((*psolvationparameters).meshsize.y+i)*(*psolvationparameters).meshsize.x+j)+1].columnvolume=0.5*spacing.x*spacing.y*(box_dimension.z-0.5*vacuumthickness);
            (*pmeshtriangledata)[2*(((*psolvationparameters).meshsize.y+i)*(*psolvationparameters).meshsize.x+j)].normal.x=0;
            (*pmeshtriangledata)[2*(((*psolvationparameters).meshsize.y+i)*(*psolvationparameters).meshsize.x+j)].normal.y=0;
            (*pmeshtriangledata)[2*(((*psolvationparameters).meshsize.y+i)*(*psolvationparameters).meshsize.x+j)].normal.z=1;
            (*pmeshtriangledata)[2*(((*psolvationparameters).meshsize.y+i)*(*psolvationparameters).meshsize.x+j)+1].normal.x=0;
            (*pmeshtriangledata)[2*(((*psolvationparameters).meshsize.y+i)*(*psolvationparameters).meshsize.x+j)+1].normal.y=0;
            (*pmeshtriangledata)[2*(((*psolvationparameters).meshsize.y+i)*(*psolvationparameters).meshsize.x+j)+1].normal.z=1;
            (*pmeshbonddata)[3*(((*psolvationparameters).meshsize.y+i)*(*psolvationparameters).meshsize.x+j)].length=spacing.x;
            (*pmeshbonddata)[3*(((*psolvationparameters).meshsize.y+i)*(*psolvationparameters).meshsize.x+j)+1].length=sqrt(pow(0.5*spacing.x, 2)+pow(spacing.y, 2));
            (*pmeshbonddata)[3*(((*psolvationparameters).meshsize.y+i)*(*psolvationparameters).meshsize.x+j)+2].length=sqrt(pow(0.5*spacing.x, 2)+pow(spacing.y, 2));
            (*pmeshbonddata)[3*(((*psolvationparameters).meshsize.y+i)*(*psolvationparameters).meshsize.x+j)].meancurvature=0;
            (*pmeshbonddata)[3*(((*psolvationparameters).meshsize.y+i)*(*psolvationparameters).meshsize.x+j)+1].meancurvature=0;
            (*pmeshbonddata)[3*(((*psolvationparameters).meshsize.y+i)*(*psolvationparameters).meshsize.x+j)+2].meancurvature=0;
			(*pmeshdata)[((*psolvationparameters).meshsize.y+i)*(*psolvationparameters).meshsize.x+j].meancurvature=0;
			(*pmeshdata)[((*psolvationparameters).meshsize.y+i)*(*psolvationparameters).meshsize.x+j].meancurvaturesquared=0;
			(*pmeshdata)[((*psolvationparameters).meshsize.y+i)*(*psolvationparameters).meshsize.x+j].area=3*spacing.x*spacing.y;
		}
	}
	(*psolvationparameters).totalvolume=(*psolvationparameters).targetvolume=box_dimension.x*box_dimension.y*(box_dimension.z-vacuumthickness);
    ((*psolvationparameters).totallength)[0]=((*psolvationparameters).totallength)[3]=(*psolvationparameters).meshsize.y*box_dimension.x;
    ((*psolvationparameters).totallength)[1]=((*psolvationparameters).totallength)[2]=((*psolvationparameters).totallength)[4]=((*psolvationparameters).totallength)[5]=(*psolvationparameters).meshsize.x*box_dimension.y*sqrt(spacing.y*spacing.y+0.25*spacing.x*spacing.x)/spacing.y;
	(*psolvationparameters).longestmeshbond=sqrt(pow(0.5*spacing.x, 2)+pow(spacing.y, 2));
	if(spacing.x>(*psolvationparameters).longestmeshbond) (*psolvationparameters).longestmeshbond=spacing.x;
}

void output_xyz_mesh(double_triple *meshpositions, char *xyzname, char *sourcename, int *pframe, double siteradius, double_triple box_dimension, solvation_parameters solvp){
    if(solvp.interface!=1) return;
	int newsourcefile=0, i, j, l, index;
    declare_array_nozero(double_triple, lastrow, solvp.meshsize.x);
    declare_array_nozero(double_triple, zeropoint, 2);
    double_triple thispos;
    for(i=0;i<2;i++) zeropoint[i].x=zeropoint[i].y=0;
    zeropoint[0].z=0;
    zeropoint[1].z=box_dimension.z;
    FILE *outp, *sourceoutp;
	if((*pframe)==0){
        newsourcefile=1;
		sourceoutp=fopen(sourcename, "w");
		fprintf(sourceoutp, "pbc set {%.6f %.6f %.6f} -first %i -last %i\n", box_dimension.x, box_dimension.y, box_dimension.z, *pframe, *pframe);
		fprintf(sourceoutp, "topo clearbonds\n");
		fprintf(sourceoutp, "set sel [atomselect top \" name I\"]\n$sel set radius %.6f\n", siteradius);
		fprintf(sourceoutp, "color Name I blue\n");
		outp=fopen(xyzname, "w");
	}
	else{        
		if ((sourceoutp = fopen(sourcename, "r")) == NULL) {
            newsourcefile=1;
			sourceoutp=fopen(sourcename, "w");
			fprintf(sourceoutp, "topo clearbonds\n");
			fprintf(sourceoutp, "set sel [atomselect top \" name I\"]\n$sel set radius %.6f\n", siteradius);
		} else {
			fclose(sourceoutp);
			sourceoutp=fopen(sourcename, "a");
		}
		fprintf(sourceoutp, "pbc set {%.6f %.6f %.6f} -first %i -last %i\n", box_dimension.x, box_dimension.y, box_dimension.z, *pframe, *pframe);
		if(newsourcefile==0) fclose(sourceoutp);
		outp=fopen(xyzname, "a");
	}
	fprintf(outp, "%i\n%i\t%i\t%.6f\t%.6f\t%.6f\n", solvp.meshnumber, solvp.meshsize.x, solvp.meshsize.y, box_dimension.x, box_dimension.y, box_dimension.z);
	for(l=0;l<2;l++){
        for(i=0;i<solvp.meshsize.y;i++){
            for(j=0;j<solvp.meshsize.x;j++){
                index=(l*solvp.meshsize.y+i)*solvp.meshsize.x+j;
                if(i==0){
                    if(j==0) thispos=nearest_image(meshpositions[index], zeropoint[l], box_dimension);
                    else thispos=nearest_image(meshpositions[index], lastrow[j-1], box_dimension);
                    lastrow[j]=thispos;
                }
                else{
                    thispos=nearest_image(meshpositions[index], lastrow[j], box_dimension);
                    lastrow[j]=thispos;
                }
                outputcoords_centered(outp, "I", thispos, box_dimension);				//	backbone N
                if(newsourcefile==1){
                    if(j>0){
                        fprintf(sourceoutp, "topo addbond %i %i\n", index-1, index);
                    }
                    if(i>0){
                        if(mod(i, 2)==0){
                            if(j>0){
                                fprintf(sourceoutp, "topo addbond %i %i\n", index-solvp.meshsize.x-1, index);
                            }
                            fprintf(sourceoutp, "topo addbond %i %i\n", index-solvp.meshsize.x, index);
                        }
                        else{
                            fprintf(sourceoutp, "topo addbond %i %i\n", index-solvp.meshsize.x, index);
                            if(j<solvp.meshsize.x-1){
                                fprintf(sourceoutp, "topo addbond %i %i\n", index-solvp.meshsize.x+1, index);
                            }
                        }
                    }
                }
            }
        }
    }
    if(newsourcefile==1){
        fclose(sourceoutp);
    }
    fclose(outp);
    (*pframe)++;
    free(lastrow);
    free(zeropoint);
}

void cellconstructasymmetric(double_triple *positions, int N, int reverseN, linkedlistasymmetric *plink, coord *coordarray){
	int i, j, k, n;
	for(i=0;i<(*plink).core.cellsperside.x;i++){
		for(j=0;j<(*plink).core.cellsperside.y;j++){
			for(k=0;k<(*plink).core.cellsperside.z;k++){
				(*plink).core.head[i][j][k]=-1;						//	empty
				(*plink).headreverse[i][j][k]=-1;						//	empty
			}
		}
	}
	for(n=0;n<N;n++) (*plink).reverselist[n]=-1;
	for(n=0;n<reverseN;n++) (*plink).reverselist_reverse[n]=-1;
	for(n=0;n<N;n++){
		i=(int) floor(positions[n].x/(*plink).cellwidth.x);
		j=(int) floor(positions[n].y/(*plink).cellwidth.y);
		k=(int) floor(positions[n].z/(*plink).cellwidth.z);
		i=mod(i, (*plink).core.cellsperside.x);
		j=mod(j, (*plink).core.cellsperside.y);
		k=mod(k, (*plink).core.cellsperside.z);
		(*plink).core.cell[n].x=i;
		(*plink).core.cell[n].y=j;
		(*plink).core.cell[n].z=k;
		(*plink).core.list[n]=(*plink).core.head[i][j][k];
        if((*plink).core.list[n]!=-1) (*plink).reverselist[(*plink).core.list[n]]=n;
		(*plink).core.head[i][j][k]=n;
	}
	for(n=0;n<reverseN;n++){
		i=(int) floor(coordarray[n].r.x/(*plink).cellwidth.x);
		j=(int) floor(coordarray[n].r.y/(*plink).cellwidth.y);
		k=(int) floor(coordarray[n].r.z/(*plink).cellwidth.z);
		i=mod(i, (*plink).core.cellsperside.x);
		j=mod(j, (*plink).core.cellsperside.y);
		k=mod(k, (*plink).core.cellsperside.z);
		(*plink).cellreverse[n].x=i;
		(*plink).cellreverse[n].y=j;
		(*plink).cellreverse[n].z=k;
		(*plink).listreverse[n]=(*plink).headreverse[i][j][k];
        if((*plink).listreverse[n]!=-1) (*plink).reverselist_reverse[(*plink).listreverse[n]]=n;
		(*plink).headreverse[i][j][k]=n;
	}
}

void allocate_linklist_asymmetric(linkedlistasymmetric *plink, int Nnodes, int reversenodes){
	int i, j, k;
	allocate_3d_tensor(int, (*plink).core.head, (*plink).core.cellsperside.x, (*plink).core.cellsperside.y, (*plink).core.cellsperside.z);
	allocate_array(int, (*plink).core.list, Nnodes);
	allocate_array(int, (*plink).reverselist, Nnodes);
	allocate_array(int, (*plink).reverselist_reverse, reversenodes);
	allocate_array(int_triple, (*plink).core.cell, Nnodes);
	
	allocate_array(int, (*plink).core.number_neighbors, reversenodes);
	allocate_matrix(int, (*plink).core.neighbor, reversenodes, (*plink).maxneighbors);
    
	allocate_3d_tensor(int, (*plink).headreverse, (*plink).core.cellsperside.x, (*plink).core.cellsperside.y, (*plink).core.cellsperside.z);
	allocate_array(int, (*plink).listreverse, reversenodes);
	allocate_array(int_triple, (*plink).cellreverse, reversenodes);
	allocate_array(int, (*plink).number_neighbors_reverse, Nnodes);
	allocate_matrix(int, (*plink).neighbor_reverse, Nnodes, (*plink).maxreverseneighbors);
}

void free_linklist_asymmetric(linkedlistasymmetric link, int Nnodes, int reversenodes){
    int i, j;
    free_3d_tensor(link.core.head, link.core.cellsperside.x, link.core.cellsperside.y);
    free(link.core.list);
    free(link.reverselist);
    free(link.reverselist_reverse);
    free(link.core.cell);
    free(link.core.number_neighbors);
    free_matrix(link.core.neighbor, reversenodes);
    free_3d_tensor(link.headreverse, link.core.cellsperside.x, link.core.cellsperside.y);
    free(link.listreverse);
    free(link.cellreverse);
    free(link.number_neighbors_reverse);
    free_matrix(link.neighbor_reverse, Nnodes);
}

void constructneighborlistasymmetric(linkedlistasymmetric *plink, int N, int Nreverse, coord *coordarray, double_triple box_dimension, double_triple *position){
    int i, j, k, self, neighbor;
	double_triple selfpos, selfposshifted;
	int_triple selfcell, cellshift, neighborcell;
	for(self=0;self<Nreverse;self++){
		(*plink).number_neighbors_reverse[self]=0;
	}
	for(self=0;self<N;self++){
		(*plink).core.number_neighbors[self]=0;
		selfpos=coordarray[self].r;
        selfcell=(*plink).cellreverse[self];
		for(cellshift.x=-1;cellshift.x<=1;cellshift.x++){
			if((selfcell.x==0)&&(cellshift.x==-1)){
				neighborcell.x=(*plink).core.cellsperside.x-1;
				selfposshifted.x=selfpos.x+box_dimension.x;
			}
			else if((selfcell.x==(*plink).core.cellsperside.x-1)&&(cellshift.x==1)){
				neighborcell.x=0;
				selfposshifted.x=selfpos.x-box_dimension.x;
			}
			else{
				neighborcell.x=selfcell.x+cellshift.x;
				selfposshifted.x=selfpos.x;
			}
			for(cellshift.y=-1;cellshift.y<=1;cellshift.y++){
				if((selfcell.y==0)&&(cellshift.y==-1)){
					neighborcell.y=(*plink).core.cellsperside.y-1;
					selfposshifted.y=selfpos.y+box_dimension.y;
				}
				else if((selfcell.y==(*plink).core.cellsperside.y-1)&&(cellshift.y==1)){
					neighborcell.y=0;
					selfposshifted.y=selfpos.y-box_dimension.y;
				}
				else{
					neighborcell.y=selfcell.y+cellshift.y;
					selfposshifted.y=selfpos.y;
				}
				for(cellshift.z=-1;cellshift.z<=1;cellshift.z++){
					if((selfcell.z==0)&&(cellshift.z==-1)){
						neighborcell.z=(*plink).core.cellsperside.z-1;
						selfposshifted.z=selfpos.z+box_dimension.z;
					}
					else if((selfcell.z==(*plink).core.cellsperside.z-1)&&(cellshift.z==1)){
						neighborcell.z=0;
						selfposshifted.z=selfpos.z-box_dimension.z;
					}
					else{
						neighborcell.z=selfcell.z+cellshift.z;
						selfposshifted.z=selfpos.z;
					}
					neighbor=(*plink).core.head[neighborcell.x][neighborcell.y][neighborcell.z];
					while(neighbor!=-1){
						if(norm(subtract_double_triple(position[neighbor], selfposshifted))<(*plink).mincellwidth){
							(*plink).core.neighbor[self][(*plink).core.number_neighbors[self]]=neighbor;
							(*plink).core.number_neighbors[self]++;
							if((*plink).core.number_neighbors[self]>(*plink).maxneighbors){
								printf("constructneighborlistasymmetric: number_neighbors[%i]>%i!\n", self, (*plink).maxneighbors);
								exit(1);
							}
							(*plink).neighbor_reverse[neighbor][(*plink).number_neighbors_reverse[neighbor]]=self;
							(*plink).number_neighbors_reverse[neighbor]++;
							if((*plink).number_neighbors_reverse[neighbor]>(*plink).maxneighbors){
								printf("constructneighborlistasymmetric: number_neighbors_reverse[%i]>%i!\n", neighbor, (*plink).maxneighbors);
								exit(1);
							}
						}
						neighbor=(*plink).core.list[neighbor];
					}
				}
			}
		}
	}
}

double calc_onebody_energy(coord coord_new, onebodyparam p, solvation_parameters solvp){
	double affinityheight, interfaceheight;
    affinityheight=coord_new.height-p.z0;
    interfaceheight=coord_new.height-p.zinterface;
	double newaffinity=-p.solvationenergy-roomTinkcals*log(1+(exp(-p.solvationenergy/roomTinkcals)-1)/(1+exp(-4*affinityheight/solvp.interfacethickness)));
	double newinterface=-p.uinterface*exp(-0.5*pow(interfaceheight/p.sigmainterface, 2));
	return (newaffinity+newinterface);
}

double height_relative_to_interface_recreatecell(coord mycoord, solvation_parameters solvp, int n, linkedlistasymmetric meshlink, double_triple box_dimension, double_triple *meshpositions, int **mapmeshtoneighbors){
    fmod_double_triple(&(mycoord.r), box_dimension);
	int_triple newcell=cellofpos(mycoord.r, meshlink.cellwidth), imagecell;
    int trialneighbor, a, b, c, minindex=-1, i, sign, closesttriangle, leftneighbor;
    double result=0, height, sep2, distance, distancealongedge, leftmag, rightleftmag;
    double mindistance2=pow(solvp.meshsurfacecutoff, 2);
    double_triple sep, relativepos, center, left, right, leftdirector, rightdirector, normal, updirector, postriangleframe, nearestpoint, rightsep, leftsep, lefttriangleframe, righttriangleframe, rightleftdirector, rightleftsep;
    
    for(a=-1;a<=1;a++){                                                                 //  as in height_relative_to_interface, but looping through cells not neighbors
		if((newcell.x==0)&&(a==-1)) imagecell.x=meshlink.core.cellsperside.x-1;
		else if((newcell.x==meshlink.core.cellsperside.x-1)&&(a==1)) imagecell.x=0;
		else imagecell.x=newcell.x+a;
		for(b=-1;b<=1;b++){
			if((newcell.y==0)&&(b==-1)) imagecell.y=meshlink.core.cellsperside.y-1;
			else if((newcell.y==meshlink.core.cellsperside.y-1)&&(b==1)) imagecell.y=0;
			else imagecell.y=newcell.y+b;
			for(c=-1;c<=1;c++){
				if((newcell.z==0)&&(c==-1)) imagecell.z=meshlink.core.cellsperside.z-1;
				else if((newcell.z==meshlink.core.cellsperside.z-1)&&(c==1)) imagecell.z=0;
				else imagecell.z=newcell.z+c;
				trialneighbor=meshlink.core.head[imagecell.x][imagecell.y][imagecell.z];
				while(trialneighbor>=0){
                    sep=subtract_double_triple(mycoord.r, meshpositions[trialneighbor]);
                    recenter_double_triple(&sep, box_dimension);
                    sep2=dot_product(sep, sep);
                    if(sep2<mindistance2){
                        mindistance2=sep2;
                        minindex=trialneighbor;
                        relativepos=sep;
                    }
					trialneighbor=meshlink.core.list[trialneighbor];
				}
			}
		}
	}
	if(minindex==-1){						//	assign height=-meshsurfacecutoff (ASSUME NEGATIVE)
        
        height=-solvp.meshsurfacecutoff;
		return height;
	}
	
	//	project onto each of six triangles; calculate minimum distance to each triangle; keep smallest (following D. A. Simon's Ph. D. dissertation)
	
	center=meshpositions[minindex];
	for(i=0;i<6;i++){
		right=meshpositions[mapmeshtoneighbors[minindex][i]];
		rightsep=subtract_double_triple(right, center);
		recenter_double_triple(&rightsep, box_dimension);
		rightdirector=rightsep;
		righttriangleframe.x=norm(rightdirector);
		rightdirector=scalar_multiply_double_triple(rightdirector, 1./righttriangleframe.x);
		
		if(i==5) left=meshpositions[mapmeshtoneighbors[minindex][0]];
		else left=meshpositions[mapmeshtoneighbors[minindex][i+1]];
		leftsep=subtract_double_triple(left, center);
		recenter_double_triple(&leftsep, box_dimension);
        leftmag=norm(leftsep);
        leftdirector=scalar_multiply_double_triple(leftsep, 1./leftmag);
		
		normal=cross_product(rightdirector, leftdirector);
		normalize(&normal);
		updirector=cross_product(normal, rightdirector);
		
		lefttriangleframe.x=dot_product(rightdirector, leftsep);
		lefttriangleframe.y=dot_product(updirector, leftsep);
		lefttriangleframe.z=0;
		righttriangleframe.y=0;
		righttriangleframe.z=0;
		
		//	calculate nearest point in triangle, in the frame of the triangle
		
		postriangleframe.x=dot_product(rightdirector, relativepos);
		postriangleframe.y=dot_product(updirector, relativepos);
		postriangleframe.z=dot_product(normal, relativepos);
		nearestpoint.z=0;
		if(postriangleframe.y<0){
			if(postriangleframe.x<0){
                distance=norm(postriangleframe);
                if(postriangleframe.z<0) distance*=-1;
			}
            else if(postriangleframe.x>righttriangleframe.x){
                distance=sqrt(pow(postriangleframe.x-righttriangleframe.x, 2)+pow(postriangleframe.y, 2)+pow(postriangleframe.z, 2));
                if(postriangleframe.z<0) distance*=-1;
            }
			else{
                distance=sqrt(pow(postriangleframe.y, 2)+pow(postriangleframe.z, 2));
                if(postriangleframe.z<0) distance*=-1;
            }
		}
		else if(postriangleframe.y>postriangleframe.x*lefttriangleframe.y/lefttriangleframe.x){		// to the left of center-left edge
            leftdirector=lefttriangleframe;     //   now in triangle frame
            normalize(&leftdirector);            
			distancealongedge=dot_product(postriangleframe, leftdirector);            
			if(distancealongedge<0){
                distance=norm(postriangleframe);
                if(postriangleframe.z<0) distance*=-1;
			}
			else if(distancealongedge>leftmag){
                distance=sqrt(pow(postriangleframe.x-lefttriangleframe.x, 2)+pow(lefttriangleframe.y-lefttriangleframe.y, 2)+pow(postriangleframe.z, 2));
                if(postriangleframe.z<0) distance*=-1;
			}
			else{
				nearestpoint=scalar_multiply_double_triple(leftdirector, distancealongedge);        //   do use it here
                distance=sqrt(pow(postriangleframe.x-nearestpoint.x, 2)+pow(lefttriangleframe.y-nearestpoint.y, 2)+pow(postriangleframe.z, 2));
                if(postriangleframe.z<0) distance*=-1;
			}
		}
		else if(postriangleframe.y>-(postriangleframe.x-righttriangleframe.x)*lefttriangleframe.y/(righttriangleframe.x-lefttriangleframe.x)){		//	to the outside of left-right edge
			rightleftsep=subtract_double_triple(lefttriangleframe, righttriangleframe);
            rightleftmag=norm(rightleftsep);
            rightleftdirector=scalar_multiply_double_triple(rightleftsep, 1./rightleftmag);
			distancealongedge=dot_product(subtract_double_triple(postriangleframe, righttriangleframe), rightleftdirector);
			if(distancealongedge<0){
                distance=sqrt(pow(postriangleframe.x-righttriangleframe.x, 2)+pow(postriangleframe.y, 2)+pow(postriangleframe.z, 2));
                if(postriangleframe.z<0) distance*=-1;
            }
			else if(distancealongedge>rightleftmag){
                distance=sqrt(pow(postriangleframe.x-lefttriangleframe.x, 2)+pow(postriangleframe.y-lefttriangleframe.y, 2)+pow(postriangleframe.z, 2));
                if(postriangleframe.z<0) distance*=-1;
            }
			else{
                nearestpoint=add_double_triple(righttriangleframe, scalar_multiply_double_triple(rightleftdirector, distancealongedge));
                distance=sqrt(pow(postriangleframe.x-nearestpoint.x, 2)+pow(lefttriangleframe.y-nearestpoint.y, 2)+pow(postriangleframe.z, 2));
                if(postriangleframe.z<0) distance*=-1;
            }
		}
		else{
            distance=postriangleframe.z;
		}
        if(i==0){
            height=distance;
        }
        else if(distance*distance<height*height){
            height=distance;
        }
	}    
    if(minindex<solvp.meshsize.x*solvp.meshsize.y){
        height*=-1;
    }
    if(height>solvp.meshsurfacecutoff) height=solvp.meshsurfacecutoff;
    else if(height<-solvp.meshsurfacecutoff) height=-solvp.meshsurfacecutoff;
    return height;
}

void calc_cellstruct_leafid_fixedinterface(coord *coordarray, bonded_params my_bonded_params, nonbonded_params my_nonbonded_params, int Nchains, int *chainlength, monomernodes **monomerid, double_triple box_dimension, energycomponents *pmy_energycomponents, linkedlistset linkset, double *avheight, double_triple *pavboxdimension){
	int i, j, left=-1, twoleft=-1, threeleft=-1, fourleft=-1, backboneid, sidechainid;
	double height;
	(*pavboxdimension).x+=box_dimension.x;
	(*pavboxdimension).y+=box_dimension.y;
	(*pavboxdimension).z+=box_dimension.z;
	for(i=0;i<Nchains;i++){
		if(chainlength[i]>2){
			for(j=0;j<chainlength[i];j++){
				
                backboneid=monomerid[i][j].backbone;
				if(my_nonbonded_params.solvationparams.interface==1) calc_onebody_energy_components(coordarray[backboneid], chooseonebody(coordarray[backboneid].nodetype, my_nonbonded_params), my_nonbonded_params.solvationparams, pmy_energycomponents);
				
				//	using old sidechainid:
				//	backbone only has short-ranged interactions (zero)
                
				sidechainid=monomerid[i][j].sidechain;
                if(coordarray[sidechainid].nodetype>=0){
                    if(my_nonbonded_params.solvationparams.interface==1){
                        height=coordarray[backboneid].r.z;
                        if(height>0.5*(my_nonbonded_params.solvationparams.interfaceheight1+my_nonbonded_params.solvationparams.interfaceheight2)) height=height-my_nonbonded_params.solvationparams.interfaceheight2;
                        else height=my_nonbonded_params.solvationparams.interfaceheight1-height;
                        avheight[coordarray[backboneid].nodetype]+=height;
                        height=coordarray[sidechainid].r.z;
                        if(height>0.5*(my_nonbonded_params.solvationparams.interfaceheight1+my_nonbonded_params.solvationparams.interfaceheight2)) height=height-my_nonbonded_params.solvationparams.interfaceheight2;
                        else height=my_nonbonded_params.solvationparams.interfaceheight1-height;
                        avheight[coordarray[sidechainid].nodetype]+=height;
                    }
                    if(my_nonbonded_params.solvationparams.interface==1) calc_onebody_energy_components(coordarray[sidechainid], chooseonebody(coordarray[sidechainid].nodetype, my_nonbonded_params), my_nonbonded_params.solvationparams, pmy_energycomponents);
                    
                    if(coordarray[sidechainid].nodetype==1){
						calc_nonbonded_energy_components_cellstruct_leafid(sidechainid, backboneid, -1, -1, linkset.phenylphenyl.core, coordarray, box_dimension, my_nonbonded_params, pmy_energycomponents);
                    }
					if(coordarray[sidechainid].nodetype>1){
						calc_nonbonded_energy_components_cellstruct_leafid(sidechainid, backboneid, -1, -1, linkset.chargedcharged.core, coordarray, box_dimension, my_nonbonded_params, pmy_energycomponents);
					}
                    if(coordarray[sidechainid].nodetype>=0) (*pmy_energycomponents).sidechain+=calc_sidechain_bonded_energy(my_bonded_params.sidechain[coordarray[sidechainid].nodetype], coordarray[backboneid], coordarray[sidechainid], box_dimension);
                }
				if(j==2){
                    (*pmy_energycomponents).backbone+=calc_backbone_bonded_energy(my_bonded_params, coordarray[left], coordarray[twoleft], coordarray[backboneid], box_dimension);
                }
				else if(j>=3){
                    (*pmy_energycomponents).backbone+=calc_backbone_bonded_energy_norim1i(my_bonded_params, coordarray[left], coordarray[twoleft], coordarray[backboneid], box_dimension);
                }
				fourleft=threeleft;
				threeleft=twoleft;
				twoleft=left;
				left=backboneid;
			}
            
			//	backbone only has short-ranged interactions (zero)
            
		}
		else if(chainlength[i]==2){
			for(j=0;j<chainlength[i];j++){
				backboneid=monomerid[i][j].backbone;
				if(my_nonbonded_params.solvationparams.interface==1) calc_onebody_energy_components(coordarray[backboneid], chooseonebody(coordarray[backboneid].nodetype, my_nonbonded_params), my_nonbonded_params.solvationparams, pmy_energycomponents);
				sidechainid=monomerid[i][j].sidechain;
                if(coordarray[sidechainid].nodetype>=0){
                    if(my_nonbonded_params.solvationparams.interface==1){
                        height=coordarray[backboneid].r.z;
                        if(height>0.5*(my_nonbonded_params.solvationparams.interfaceheight1+my_nonbonded_params.solvationparams.interfaceheight2)) height=height-my_nonbonded_params.solvationparams.interfaceheight2;
                        else height=my_nonbonded_params.solvationparams.interfaceheight1-height;
                        avheight[coordarray[backboneid].nodetype]+=height;
                        height=coordarray[sidechainid].r.z;
                        if(height>0.5*(my_nonbonded_params.solvationparams.interfaceheight1+my_nonbonded_params.solvationparams.interfaceheight2)) height=height-my_nonbonded_params.solvationparams.interfaceheight2;
                        else height=my_nonbonded_params.solvationparams.interfaceheight1-height;
                        avheight[coordarray[sidechainid].nodetype]+=height;
                    }
                    if(my_nonbonded_params.solvationparams.interface==1) calc_onebody_energy_components(coordarray[sidechainid], chooseonebody(coordarray[sidechainid].nodetype, my_nonbonded_params), my_nonbonded_params.solvationparams, pmy_energycomponents);
                    
                    if(coordarray[sidechainid].nodetype==1) calc_nonbonded_energy_components_cellstruct_leafid(sidechainid, backboneid, -1, -1, linkset.phenylphenyl.core, coordarray, box_dimension, my_nonbonded_params, pmy_energycomponents);
                    if(coordarray[sidechainid].nodetype>1) calc_nonbonded_energy_components_cellstruct_leafid(sidechainid, backboneid, -1, -1, linkset.chargedcharged.core, coordarray, box_dimension, my_nonbonded_params, pmy_energycomponents);
                    
                    if(coordarray[sidechainid].nodetype>=0) (*pmy_energycomponents).sidechain+=calc_sidechain_bonded_energy(my_bonded_params.sidechain[coordarray[sidechainid].nodetype], coordarray[backboneid], coordarray[sidechainid], box_dimension);
                }
				if(j==1){
					
					//	backbone only has short-ranged interactions (zero)
                    
					(*pmy_energycomponents).backbone+=calc_backbone_bonded_energy_onlyleft(my_bonded_params, coordarray[backboneid], coordarray[left], box_dimension);
				}
				left=backboneid;
			}
			
			//	backbone only has short-ranged interactions (zero)
            
		}
		else{
			for(j=0;j<chainlength[i];j++){
				backboneid=monomerid[i][j].backbone;
				if(my_nonbonded_params.solvationparams.interface==1) calc_onebody_energy_components(coordarray[backboneid], chooseonebody(coordarray[backboneid].nodetype, my_nonbonded_params), my_nonbonded_params.solvationparams, pmy_energycomponents);
				sidechainid=monomerid[i][j].sidechain;
                if(coordarray[sidechainid].nodetype>=0){
                    if(my_nonbonded_params.solvationparams.interface==1){
                        height=coordarray[backboneid].r.z;
                        if(height>0.5*(my_nonbonded_params.solvationparams.interfaceheight1+my_nonbonded_params.solvationparams.interfaceheight2)) height=height-my_nonbonded_params.solvationparams.interfaceheight2;
                        else height=my_nonbonded_params.solvationparams.interfaceheight1-height;
                        avheight[coordarray[backboneid].nodetype]+=height;
                        height=coordarray[sidechainid].r.z;
                        if(height>0.5*(my_nonbonded_params.solvationparams.interfaceheight1+my_nonbonded_params.solvationparams.interfaceheight2)) height=height-my_nonbonded_params.solvationparams.interfaceheight2;
                        else height=my_nonbonded_params.solvationparams.interfaceheight1-height;
                        avheight[coordarray[sidechainid].nodetype]+=height;
                    }
                    if(my_nonbonded_params.solvationparams.interface==1) calc_onebody_energy_components(coordarray[sidechainid], chooseonebody(coordarray[sidechainid].nodetype, my_nonbonded_params), my_nonbonded_params.solvationparams, pmy_energycomponents);
                    
                    if(coordarray[sidechainid].nodetype==1) calc_nonbonded_energy_components_cellstruct_leafid(sidechainid, backboneid, -1, -1, linkset.phenylphenyl.core, coordarray, box_dimension, my_nonbonded_params, pmy_energycomponents);
                    if(coordarray[sidechainid].nodetype>1) calc_nonbonded_energy_components_cellstruct_leafid(sidechainid, backboneid, -1, -1, linkset.chargedcharged.core, coordarray, box_dimension, my_nonbonded_params, pmy_energycomponents);
                    
                    if(coordarray[sidechainid].nodetype>=0) (*pmy_energycomponents).sidechain+=calc_sidechain_bonded_energy(my_bonded_params.sidechain[coordarray[sidechainid].nodetype], coordarray[backboneid], coordarray[sidechainid], box_dimension);
                }
			}
        }
	}
}

void output_timeseries_fixedinterface(char *totalfilename, energycomponents *pmy_energycomponents, int cycle, monomertypes monomercount, double factor, int Nnodetypes, double *avheight, double_triple *pavboxdimension, int runningtime){
	int i;
	FILE *outp;
	outp=fopen(totalfilename, "a");
	fprintf(outp, "%i\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f", cycle, (*pmy_energycomponents).backbone*factor/monomercount.monomers, (*pmy_energycomponents).sidechain*factor/monomercount.monomers, 0.5*(*pmy_energycomponents).nnsamepoly*factor/monomercount.monomers, 0.5*(*pmy_energycomponents).cnsamepoly*factor/monomercount.monomers, 0.5*(*pmy_energycomponents).cclikesamepoly*factor/monomercount.monomers, 0.5*(*pmy_energycomponents).ccunlikesamepoly*factor/monomercount.monomers, 0.5*(*pmy_energycomponents).nn*factor/monomercount.monomers, 0.5*(*pmy_energycomponents).cn*factor/monomercount.monomers, 0.5*(*pmy_energycomponents).cclike*factor/monomercount.monomers, 0.5*(*pmy_energycomponents).ccunlike*factor/monomercount.monomers);
	for(i=0;i<Nnodetypes;i++){
        fprintf(outp, "\t%.8f\t%.8f", (*pmy_energycomponents).solvation[i]*factor/monomercount.monomers, (*pmy_energycomponents).interface[i]*factor/monomercount.monomers);
		(*pmy_energycomponents).solvation[i]=0;
		(*pmy_energycomponents).interface[i]=0;
    }
	for(i=0;i<Nnodetypes;i++){
		fprintf(outp, "\t%.8f", avheight[i]*factor/monomercount.type[i]);
        avheight[i]=0;
	}
	fprintf(outp, "\t%.8f\t%.8f\t%.8f", (*pavboxdimension).x*factor, (*pavboxdimension).y*factor, (*pavboxdimension).z*factor);
	(*pavboxdimension).x=(*pavboxdimension).y=(*pavboxdimension).z=0;
	fprintf(outp, "\t%i\n", runningtime);
	fclose(outp);
	(*pmy_energycomponents).backbone=0;
	(*pmy_energycomponents).sidechain=0;
	(*pmy_energycomponents).nn=0;
	(*pmy_energycomponents).ccunlike=0;
	(*pmy_energycomponents).cclike=0;
	(*pmy_energycomponents).cn=0;
	(*pmy_energycomponents).nnsamepoly=0;
	(*pmy_energycomponents).ccunlikesamepoly=0;
	(*pmy_energycomponents).cclikesamepoly=0;
	(*pmy_energycomponents).cnsamepoly=0;
}

void calc_interface_avheight(int meshnumber, double_triple *meshpositions, double *avmeshheight){
    int i, j;
    for(i=0;i<2;i++){
        for(j=0;j<meshnumber/2;j++){
            avmeshheight[i]+=meshpositions[i*meshnumber/2+j].z;
        }
    }
}

void input_coords_hex(char *filename, int *pNnodes, int *pNchains, int **chainlength, monomernodes ***monomerid, coord **coordarray, double_triple *pbox_dimension, double *pvacuumthickness, int *pNnodetypes, monomertypes *pmonomercount, long long int *pt, int *pframe, int *prunningtime, double_triple **pmeshpositions, triangledata **pmeshtriangledata, bonddata **pmeshbonddata, vertexdata **pmeshdata, solvation_parameters *psolvp){
	int i, j;
    double val;
	FILE *inp;
	inp=fopen(filename, "r");
	fscanf(inp, "%lf %lf %lf", &((*pbox_dimension).x), &((*pbox_dimension).y), &((*pbox_dimension).z));
	fscanf(inp, "%i", &((*psolvp).interface));
	if(((*psolvp).interface)==1)fscanf(inp, " %lf", pvacuumthickness);
	fscanf(inp, "%i %i %i %i %i %i", &(*pNnodes), &(*pNchains), &((*pmonomercount).charged), &((*pmonomercount).nonpolar), &((*pmonomercount).monomers), &(*pNnodetypes));
	for(i=0;i<(*pNnodetypes);i++) fscanf(inp, " %i", &((*pmonomercount).type[i]));
	(*coordarray)=xcalloc((*pNnodes), sizeof(coord));
	(*chainlength)=xcalloc((*pNchains), sizeof(int));
	(*monomerid)=xcalloc((*pNchains), sizeof(monomernodes *));
	for(i=0;i<*pNchains;i++){
		fscanf(inp, "%i", &((*chainlength)[i]));
		((*monomerid)[i])=xcalloc(((*chainlength)[i]), sizeof(monomernodes));
		for(j=0;j<(*chainlength)[i];j++){
			fscanf(inp, "%i %i", &((*monomerid)[i][j].backbone), &((*monomerid)[i][j].sidechain));
		}
	}
	for(i=0;i<(*pNnodes);i++){
		fscanf(inp, "%lf %lf %lf %lf %lf %lf %i %i %i %i %i", &((*coordarray)[i].r.x), &((*coordarray)[i].r.y), &((*coordarray)[i].r.z), &((*coordarray)[i].n.x), &((*coordarray)[i].n.y), &((*coordarray)[i].n.z), &((*coordarray)[i].nodetype), &((*coordarray)[i].chainid), &((*coordarray)[i].monomerid), &((*coordarray)[i].orientationtype), &((*coordarray)[i].leafid));
	}
    fscanf(inp, "%lli %i %i", &(*pt), &(*pframe), &(*prunningtime));
	if(((*psolvp).interface)==1){
		fscanf(inp, "%i %i", &(((*psolvp).meshsize).x), &(((*psolvp).meshsize).y));
		((*psolvp).meshnumber)=2*((*psolvp).meshsize).x*((*psolvp).meshsize).y;
		allocate_array(double_triple, (*pmeshpositions), ((*psolvp).meshnumber));
		allocate_array(triangledata, (*pmeshtriangledata), (2*(*psolvp).meshnumber));
		allocate_array(bonddata, (*pmeshbonddata), (3*(*psolvp).meshnumber));
		allocate_array(vertexdata, (*pmeshdata), ((*psolvp).meshnumber));
		fscanf(inp, "%lf", &((*psolvp).meshbondlength));
		fscanf(inp, "%lf", &((*psolvp).maxmeshbond));
		fscanf(inp, "%lf", &((*psolvp).longestmeshbond));
		fscanf(inp, "%lf", &((*psolvp).totalvolume));
		fscanf(inp, "%lf", &((*psolvp).targetvolume));
		
		fscanf(inp, "%lf %lf", &((*psolvp).totallength[0]), &((*psolvp).totallength[1]));
        fscanf(inp, "%lf", &val);
        if(val>0){
            (*psolvp).totallength[2]=val;
            fscanf(inp, "%lf %lf %lf", &((*psolvp).totallength[3]), &((*psolvp).totallength[4]), &((*psolvp).totallength[5]));
            j=0;
        }
        else{
            i=0;
            j=0;
            (*pmeshpositions)[i*((*psolvp).meshsize).x+j].x=val;
            fscanf(inp, "%lf %lf", &((*pmeshpositions)[i*((*psolvp).meshsize).x+j].y), &((*pmeshpositions)[i*((*psolvp).meshsize).x+j].z));
            j=1;
        }
		for(i=0;i<((*psolvp).meshsize).y;i++){
			for(;j<((*psolvp).meshsize).x;j++){
				fscanf(inp, "%lf %lf %lf", &((*pmeshpositions)[i*((*psolvp).meshsize).x+j].x), &((*pmeshpositions)[i*((*psolvp).meshsize).x+j].y), &((*pmeshpositions)[i*((*psolvp).meshsize).x+j].z));
			}
            j=0;
		}
		for(i=0;i<((*psolvp).meshsize).y;i++){
			for(j=0;j<((*psolvp).meshsize).x;j++){
				fscanf(inp, "%lf %lf %lf", &((*pmeshpositions)[(((*psolvp).meshsize).y+i)*((*psolvp).meshsize).x+j].x), &((*pmeshpositions)[(((*psolvp).meshsize).y+i)*((*psolvp).meshsize).x+j].y), &((*pmeshpositions)[(((*psolvp).meshsize).y+i)*((*psolvp).meshsize).x+j].z));
			}
		}
        for(i=0;i<2*(*psolvp).meshnumber;i++){
            fscanf(inp, "%lf %lf %lf %lf %lf", &((*pmeshtriangledata)[i].area), &((*pmeshtriangledata)[i].columnvolume), &((*pmeshtriangledata)[i].normal.x), &((*pmeshtriangledata)[i].normal.y), &((*pmeshtriangledata)[i].normal.z));
        }
        for(i=0;i<3*(*psolvp).meshnumber;i++){
            fscanf(inp, "%lf %lf", &((*pmeshbonddata)[i].length), &((*pmeshbonddata)[i].meancurvature));
        }
		for(i=0;i<(*psolvp).meshnumber;i++){
            fscanf(inp, "%lf %lf %lf", &((*pmeshdata)[i].area), &((*pmeshdata)[i].meancurvature), &((*pmeshdata)[i].meancurvaturesquared));
        }
	}
}

void output_coords_hex(char *filename, int Nnodes, int Nchains, int *chainlength, monomernodes **monomerid, coord *coordarray, double_triple box_dimension, double vacuumthickness, int Nnodetypes, monomertypes monomercount, long long int t, int frame, int runningtime, double_triple *meshpositions, triangledata *meshtriangledata, bonddata *meshbonddata, vertexdata *meshdata, solvation_parameters solvationparameters){
	int i, j;
	FILE *outp;
	outp=fopen(filename, "w");
	fprintf(outp, "%llui %llui %llui\n", (unsigned long long int) box_dimension.x, (unsigned long long int) box_dimension.y, (unsigned long long int) box_dimension.z);
	fprintf(outp, "%i", solvationparameters.interface);
	if(solvationparameters.interface==1)fprintf(outp, " %llui\n", (unsigned long long int) vacuumthickness);
	else fprintf(outp, "\n");
	fprintf(outp, "%i %i %i %i %i %i", Nnodes, Nchains, monomercount.charged, monomercount.nonpolar, monomercount.monomers, Nnodetypes);
	for(i=0;i<Nnodetypes;i++) fprintf(outp, " %i", monomercount.type[i]);
	fprintf(outp, "\n");
	for(i=0;i<Nchains;i++){
		fprintf(outp, "%i", chainlength[i]);
		for(j=0;j<chainlength[i];j++){
			fprintf(outp, " %i %i", monomerid[i][j].backbone, monomerid[i][j].sidechain);
		}
		fprintf(outp, "\n");
	}
	for(i=0;i<Nnodes;i++){
		fprintf(outp, "%llui %llui %llui %llui %llui %llui %i %i %i %i %i\n", (unsigned long long int) coordarray[i].r.x, (unsigned long long int) (unsigned long long int) coordarray[i].r.y, (unsigned long long int) coordarray[i].r.z, (unsigned long long int) coordarray[i].n.x, (unsigned long long int) coordarray[i].n.y, (unsigned long long int) coordarray[i].n.z, coordarray[i].nodetype, coordarray[i].chainid, coordarray[i].monomerid, coordarray[i].orientationtype, coordarray[i].leafid);
	}
    fprintf(outp, "%lli %i %i\n", t, frame, runningtime);
	if(solvationparameters.interface==1){
		fprintf(outp, "%i %i\n", solvationparameters.meshsize.x, solvationparameters.meshsize.y);
		fprintf(outp, "%llui\n", (unsigned long long int) solvationparameters.meshbondlength);
		fprintf(outp, "%llui\n", (unsigned long long int) solvationparameters.maxmeshbond);
		fprintf(outp, "%llui\n", (unsigned long long int) solvationparameters.longestmeshbond);
		fprintf(outp, "%llui\n", (unsigned long long int) solvationparameters.totalvolume);
		fprintf(outp, "%llui\n", (unsigned long long int) solvationparameters.targetvolume);
		fprintf(outp, "%llui %llui %llui %llui %llui %llui\n", (unsigned long long int) solvationparameters.totallength[0], (unsigned long long int) solvationparameters.totallength[1], (unsigned long long int) solvationparameters.totallength[2], (unsigned long long int) solvationparameters.totallength[3], (unsigned long long int) solvationparameters.totallength[4], (unsigned long long int) solvationparameters.totallength[5]);
		for(i=0;i<solvationparameters.meshsize.y;i++){
			for(j=0;j<solvationparameters.meshsize.x;j++){
				fprintf(outp, "%llui %llui %llui\n", (unsigned long long int) meshpositions[i*solvationparameters.meshsize.x+j].x, (unsigned long long int) meshpositions[i*solvationparameters.meshsize.x+j].y, (unsigned long long int) meshpositions[i*solvationparameters.meshsize.x+j].z);
			}
		}
		for(i=0;i<solvationparameters.meshsize.y;i++){
			for(j=0;j<solvationparameters.meshsize.x;j++){
				fprintf(outp, "%llui %llui %llui\n", (unsigned long long int) meshpositions[(solvationparameters.meshsize.y+i)*solvationparameters.meshsize.x+j].x, (unsigned long long int) meshpositions[(solvationparameters.meshsize.y+i)*solvationparameters.meshsize.x+j].y, (unsigned long long int) meshpositions[(solvationparameters.meshsize.y+i)*solvationparameters.meshsize.x+j].z);
			}
		}
        for(i=0;i<2*solvationparameters.meshnumber;i++){
            fprintf(outp, "%llui %llui %llui %llui %llui\n", (unsigned long long int) meshtriangledata[i].area, (unsigned long long int) (unsigned long long int) meshtriangledata[i].columnvolume, (unsigned long long int) meshtriangledata[i].normal.x, (unsigned long long int) meshtriangledata[i].normal.y, (unsigned long long int) meshtriangledata[i].normal.z);
        }
        for(i=0;i<3*solvationparameters.meshnumber;i++){
            fprintf(outp, "%llui %llui\n", (unsigned long long int) meshbonddata[i].length, (unsigned long long int) meshbonddata[i].meancurvature);
        }
        for(i=0;i<solvationparameters.meshnumber;i++){
            fprintf(outp, "%llui %llui %llui\n",(unsigned long long int)  meshdata[i].area, (unsigned long long int) meshdata[i].meancurvature, (unsigned long long int) meshdata[i].meancurvaturesquared);
        }
	}	
	fclose(outp);
}

void output_coords_long(char *filename, int Nnodes, int Nchains, int *chainlength, monomernodes **monomerid, coord *coordarray, double_triple box_dimension, double vacuumthickness, int Nnodetypes, monomertypes monomercount, long long int t, int frame, int runningtime, double_triple *meshpositions, triangledata *meshtriangledata, bonddata *meshbonddata, vertexdata *meshdata, solvation_parameters solvationparameters){
	int i, j;
	FILE *outp;
	outp=fopen(filename, "w");
	fprintf(outp, "%.40f %.40f %.40f\n", box_dimension.x, box_dimension.y, box_dimension.z);
	fprintf(outp, "%i", solvationparameters.interface);
	if(solvationparameters.interface==1)fprintf(outp, " %.40f\n", vacuumthickness);
	else fprintf(outp, "\n");
	fprintf(outp, "%i %i %i %i %i %i", Nnodes, Nchains, monomercount.charged, monomercount.nonpolar, monomercount.monomers, Nnodetypes);
	for(i=0;i<Nnodetypes;i++) fprintf(outp, " %i", monomercount.type[i]);
	fprintf(outp, "\n");
	for(i=0;i<Nchains;i++){
		fprintf(outp, "%i", chainlength[i]);
		for(j=0;j<chainlength[i];j++){
			fprintf(outp, " %i %i", monomerid[i][j].backbone, monomerid[i][j].sidechain);
		}
		fprintf(outp, "\n");
	}
	for(i=0;i<Nnodes;i++){
		fprintf(outp, "%.40f %.40f %.40f %.40f %.40f %.40f %i %i %i %i %i\n", coordarray[i].r.x, coordarray[i].r.y, coordarray[i].r.z, coordarray[i].n.x, coordarray[i].n.y, coordarray[i].n.z, coordarray[i].nodetype, coordarray[i].chainid, coordarray[i].monomerid, coordarray[i].orientationtype, coordarray[i].leafid);
	}
    fprintf(outp, "%lli %i %i\n", t, frame, runningtime);
	if(solvationparameters.interface==1){
		fprintf(outp, "%i %i\n", solvationparameters.meshsize.x, solvationparameters.meshsize.y);
		fprintf(outp, "%.40f\n", solvationparameters.meshbondlength);
		fprintf(outp, "%.40f\n", solvationparameters.maxmeshbond);
		fprintf(outp, "%.40f\n", solvationparameters.longestmeshbond);
		fprintf(outp, "%.40f\n", solvationparameters.totalvolume);
		fprintf(outp, "%.40f\n", solvationparameters.targetvolume);
		fprintf(outp, "%.40f %.40f %.40f %.40f %.40f %.40f\n", solvationparameters.totallength[0], solvationparameters.totallength[1], solvationparameters.totallength[2], solvationparameters.totallength[3], solvationparameters.totallength[4], solvationparameters.totallength[5]);
		for(i=0;i<solvationparameters.meshsize.y;i++){
			for(j=0;j<solvationparameters.meshsize.x;j++){
				fprintf(outp, "%.40f %.40f %.40f\n", meshpositions[i*solvationparameters.meshsize.x+j].x, meshpositions[i*solvationparameters.meshsize.x+j].y, meshpositions[i*solvationparameters.meshsize.x+j].z);
			}
		}
		for(i=0;i<solvationparameters.meshsize.y;i++){
			for(j=0;j<solvationparameters.meshsize.x;j++){
				fprintf(outp, "%.40f %.40f %.40f\n", meshpositions[(solvationparameters.meshsize.y+i)*solvationparameters.meshsize.x+j].x, meshpositions[(solvationparameters.meshsize.y+i)*solvationparameters.meshsize.x+j].y, meshpositions[(solvationparameters.meshsize.y+i)*solvationparameters.meshsize.x+j].z);
			}
		}
        for(i=0;i<2*solvationparameters.meshnumber;i++){
            fprintf(outp, "%.40f %.40f %.40f %.40f %.40f\n", meshtriangledata[i].area, meshtriangledata[i].columnvolume, meshtriangledata[i].normal.x, meshtriangledata[i].normal.y, meshtriangledata[i].normal.z);
        }
        for(i=0;i<3*solvationparameters.meshnumber;i++){
            fprintf(outp, "%.40f %.40f\n", meshbonddata[i].length, meshbonddata[i].meancurvature);
        }
        for(i=0;i<solvationparameters.meshnumber;i++){
            fprintf(outp, "%.40f %.40f %.40f\n", meshdata[i].area, meshdata[i].meancurvature, meshdata[i].meancurvaturesquared);
        }
	}	
	fclose(outp);
}




