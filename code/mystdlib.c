#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "mystdlib.h"

void *xcalloc (size_t nelem,  size_t elsize){
    void *new_mem = (void*)calloc(nelem ? nelem : 1, elsize ? elsize : 1);
    if (new_mem == NULL) {
        fprintf(stderr,
                "xcalloc: request for %lu elements of size %lu failed.\n",
                (unsigned long)nelem, (unsigned long)elsize);
        exit(EXIT_FAILURE);
    }    
    return new_mem;
}

anytype loadparam(int argc, char flag[100], char *argv[], char defaultchar[100]){
	int j;
	char *temp;
	temp=defaultchar;
	for(j=1;j<argc;j++){
		if(argv[j][0]=='-'){
			if(strcmp(flag, argv[j])==0){
				temp=argv[j+1];
			}
		}
	}
	anytype tempstruct;
	tempstruct.s=temp;
	tempstruct.i=atoi(temp);
	tempstruct.ll=atoll(temp);
	tempstruct.f=atof(temp);
	return tempstruct;
}

void loadseries(int argc, char flag[100], char *argv[], char *chararray[100]){		// different from inelastic_hard_sphere
	int i, j, slot;
	for(i=0;i<argc;i++){
		if(argv[i][0]=='-'){
			if(strcmp(flag, argv[i])==0){
				slot=i+1;
				j=0;
				while(argv[slot][0]!='-'){
					chararray[j]=argv[slot];
					j++;
					slot++;
				}
			}
		}
	}
}

void my_exit(char *mychar){
	printf("%s\n", mychar);
	exit(1);
}

int rand_int(int max){
    //printf("in rand_int, max %i, rand_max %i\n", max, rand_max);
	//int max_allowed=RAND_MAX-RAND_MAX%max;
	int max_allowed=rand_max-rand_max%max;
	int temp=max_allowed;
	while(temp>=max_allowed){
		temp=rand()%max;
	}
	return temp;
}

int mygetline(char s[], int lim, FILE *inp){
	int c=0, i;
	for(i=0;i<lim-1&&(c=getc(inp))!=EOF&&c!='\n';i++){
		s[i]=c;
	}
	if(c=='\n'){
		s[i]=c;
		++i;
	}
	s[i]='\0';
	return i;
}

//  Random number generator from Matteo Italia
//  http://stackoverflow.com/questions/10198758/how-to-get-current-seed-from-c-rand
//  Use so I can retrieve the seed, useful for checkpointing, restarting, and debugging
//  But using glibc parameters from http://en.wikipedia.org/wiki/Linear_congruential_generator
//  Note that rand_max is set in mystdlib.h, though I don't seem to be using it

int myrand()
{
    //mynextseed = mynextseed * 1103515245 + 12345;
    //return (unsigned int)(mynextseed/65536) % 32768;
    mynextseed = mynextseed * 1103515245 + 12345;
    return (unsigned int)(mynextseed) % 2147483648;
}
    
void reseed(unsigned long seed)
{
    mynextseed = seed;
}
    
unsigned long getseed()
{
    return mynextseed;
}
