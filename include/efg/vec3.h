#ifndef EFG_VEC3_H
#define EFG_VEC3_H

/*
 * Minimal 3-D vector type.
 *
 * This file is part of the project skeleton. Vector operations will be added
 * only when a later step needs them for the equations documented in docs/.
 */
typedef struct EfgVec3 {
    double x;
    double y;
    double z;
} EfgVec3;

int efg_vec3_module_ready(void);

#endif
