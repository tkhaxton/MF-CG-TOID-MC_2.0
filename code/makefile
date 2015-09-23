
COMPILER=gcc
CFLAGS=-lm -O3

#CFLAGS=-lm -O0 -g

all : monolayer peptoid post_analysis bilayer_from_monolayer_enumerate bilayer_from_monolayer 

peptoid : peptoid.c mystdlib.c mystdlib.h mymath.c mymath.h peptoid_functions.c peptoid_functions.h
	${COMPILER} ${CFLAGS} -o peptoid peptoid.c peptoid_functions.c mystdlib.c mymath.c

post_analysis : post_analysis.c mystdlib.c mystdlib.h mymath.c mymath.h post_analysis_functions.c post_analysis_functions.h peptoid_functions.c peptoid_functions.h
	${COMPILER} ${CFLAGS} -o post_analysis post_analysis.c post_analysis_functions.c mymath.c mystdlib.c peptoid_functions.c

monolayer : monolayer.c mystdlib.c mystdlib.h mymath.c mymath.h peptoid_functions.c peptoid_functions.h search_functions.c search_functions.h
	${COMPILER} ${CFLAGS} -o monolayer monolayer.c peptoid_functions.c mymath.c mystdlib.c search_functions.c

bilayer_from_monolayer_enumerate : bilayer_from_monolayer_enumerate.c mystdlib.c mystdlib.h mymath.c mymath.h peptoid_functions.c peptoid_functions.h search_functions.c search_functions.h
	${COMPILER} ${CFLAGS} -o bilayer_from_monolayer_enumerate bilayer_from_monolayer_enumerate.c peptoid_functions.c mymath.c mystdlib.c search_functions.c

bilayer_from_monolayer : bilayer_from_monolayer.c mystdlib.c mystdlib.h mymath.c mymath.h peptoid_functions.c peptoid_functions.h search_functions.c search_functions.h
	${COMPILER} ${CFLAGS} -o bilayer_from_monolayer bilayer_from_monolayer.c peptoid_functions.c mymath.c mystdlib.c search_functions.c

clean:
	rm monolayer
	rm bilayer_from_monolayer
	rm bilayer_from_monolayer_enumerate
	rm peptoid
	rm post_analysis
