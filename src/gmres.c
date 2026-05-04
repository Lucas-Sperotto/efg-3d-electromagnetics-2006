#include "efg/gmres.h"

/*
 * Link-time smoke symbol for GMRES.
 *
 * TODO: implement the iterative solver only after the linear system interface
 * exists; see docs/02_formulacao_matematica.md, section 2.1.
 */
int efg_gmres_module_ready(void)
{
    return 1;
}
