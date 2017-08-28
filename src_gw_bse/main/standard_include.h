#ifndef _STANDARD_INCLUDE_
#define _STANDARD_INCLUDE_

#define CHARM_ON
#define PUP_ON
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <ctime>
#include <iostream>

#ifdef PUP_ON
#include <pup.h>
#endif
#ifdef CHARM_ON
#include "charm++.h"
#endif

#include "piny_constants.h"
#include "gwbse_constants.h"
#include "proto_friend_lib_entry.h"

#ifdef USE_LAPACK
#include "cblas.h"
#endif

#endif
