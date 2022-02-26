#ifndef GLOBDEF_H
#define GLOBDEF_H

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#ifndef RAND_MAX
#define RAND_MAX    32767
#endif
#define RaNdOmIzE   srand((unsigned)time(NULL))
#define RaNdOm(num) ((int)(( (long)rand()*(num) )/(RAND_MAX+1)))
#define MaX(n1,n2) ( (n1)>(n2) ? (n1) : (n2) )
#define MiN(n1,n2) ( (n1)<(n2) ? (n1) : (n2) )

typedef double TOMBnx101[30][101];
typedef double TOMBnx21[30][21];
typedef TOMBnx101* TOMB_nx101;
typedef TOMBnx21* TOMB_nx21;

#endif


