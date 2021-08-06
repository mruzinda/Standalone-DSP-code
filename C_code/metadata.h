#include <stdint.h>
#include <stdio.h>

#define NBEAMS 65
#define NANTS 64

struct beamCoord {
  float ra[NBEAMS];
  float dec[NBEAMS];
};

struct antCoord {
  float east[NANTS];
  float north[NANTS];
  float up[NANTS];
};