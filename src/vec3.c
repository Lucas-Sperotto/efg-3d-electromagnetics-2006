#include "efg/vec3.h"

/*
 * Link-time smoke symbol for the vector module.
 *
 * The real vector helpers will stay small and formula-free unless a later
 * equation from docs/ explicitly requires more behavior.
 */
int efg_vec3_module_ready(void)
{
    return 1;
}
