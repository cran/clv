
# "get clusters" wrappers 

agnes.wrap <- function(data.tree, clust.num=0)
{
	return( cutree(data.tree, clust.num) )
}

diana.wrap = agnes.wrap
mona.wrap = agnes.wrap # other wrapper !!!
hclust.wrap = agnes.wrap

kmeans.wrap <- function(data.cls, clust.num=0)
{
	return( data.cls$cluster )
}

pam.wrap <- function(data.cls, clust.num=0)
{
	return( data.cls )
}

clara.wrap <- function(data.cls, clust.num=0)
{
	return( data.cls$clustering )
}

# cluster methods wrappers

agnes.clust <- function(data, clust.num=0, method="empty", ... )
{
	return( agnes(data, method=method, ... ) )
}

diana.clust <- function(data, clust.num=0, method="empty", ... )
{
	return( diana(data, ... ) )
}

mona.clust <- function(data, clust.num=0, method="empty", ... )
{
	return( mona(data, ... ) )
}

hclust.clust <- function(data, clust.num=0, method="empty", ... )
{
	return( hclust( dist(data), method=method, ... ) )
}

kmeans.clust <- function(data, clust.num=0, method="empty", ... )
{
	return( kmeans(data, centers=clust.num, ... ) )
}

pam.clust <- function(data, clust.num=0, method="empty", ... )
{
	return( pam(data, clust.num, cluster.only = TRUE,  ... ) )
}

clara.clust <- function(data, clust.num=0, method="empty", ... )
{
	return( clara(data, clust.num, ... ) )
}
