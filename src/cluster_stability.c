#include "VCLDefines.h"

/*.
	file contains functions that speed up computations of cluster stablity
*/

// cluster size - 2-nd version 
// (one argument changed - function ready to work with R functions via .Call interface)

SEXP clv_clustersSizeExt(const SEXP cluster_tab_sxp, const SEXP clust_num_sxp)
{
	return clv_clustersSize(cluster_tab_sxp, INTEGER(clust_num_sxp)[0]);
}