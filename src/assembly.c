#include "efg/assembly.h"

/*
 * Link-time smoke symbol for weak-form assembly.
 *
 * TODO: add matrix/vector assembly after storage types are defined; use
 * docs/02_formulacao_matematica.md, section 2.1, as the reference.
 */
int efg_assembly_module_ready(void)
{
    return 1;
}
