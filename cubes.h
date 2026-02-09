#ifndef CUBES_H
#define CUBES_H

#if defined CUBES_ARRAY && !defined CELLS
#error CUBES_ARRAY CELLS=?
#endif

#include <iostream>
#include "dim.h"

#include "cubes-vect.h"
#include "cubes-rel.h"

#ifdef CUBES_REL
using Cubes = CubesRel;
using CubesIterator = CubesRel::CIterator;
#else
using Cubes = CubesVect;
using CubesIterator = CubesVectIterator;
#endif

#include "cubes-border.h"

#endif // CUBES_H
