#include "efg/analytical.h"
#include "efg/assembly.h"
#include "efg/error_norm.h"
#include "efg/gauss.h"
#include "efg/gmres.h"
#include "efg/mls.h"
#include "efg/nodes.h"
#include "efg/vec3.h"
#include "efg/weight.h"

/*
 * Minimal executable for the repository skeleton.
 *
 * It intentionally does not solve the EFG problem. It only includes every
 * public header and links every placeholder module so the build can fail early
 * if the project structure is broken.
 */
int main(void)
{
    const EfgVec3 origin = {0.0, 0.0, 0.0};
    const int modules_ready =
        efg_analytical_module_ready() &&
        efg_assembly_module_ready() &&
        efg_error_norm_module_ready() &&
        efg_gauss_module_ready() &&
        efg_gmres_module_ready() &&
        efg_mls_module_ready() &&
        efg_nodes_module_ready() &&
        efg_vec3_module_ready() &&
        efg_weight_module_ready();

    return (modules_ready && origin.x == 0.0 && origin.y == 0.0 && origin.z == 0.0) ? 0 : 1;
}
