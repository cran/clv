
confusion.matrix <- function(clust1,clust2)
{
	clust1 = cls.vect.validity(clust1, "clust1")
	clust2 = cls.vect.validity(clust2, "clust2")

	if( length(clust1) != length(clust2) )
		stop("Bad input data: both input vectors should have the same length.")

	clust_num1 = max(clust1)
	clust_num2 = max(clust2)
	
	
	result = .Call("confusionMatrix",
					clust1,
					clust2,
					as.integer( c(clust_num1, clust_num2) ),
					PACKAGE="clv"
				  )
	return(result)
}

similarity.index <- function(cnf.mx)
{
	cnf.mx = data.validity.int(cnf.mx, "cnf.mx")

	if( TRUE %in% (cnf.mx < 0) )
		stop("Bad input data: each 'cnf.mx' matrix element should be equal or greater than 0.")
	if( dim(cnf.mx)[1] > dim(cnf.mx)[2] ) cnf.mx = t(cnf.mx)
	
	opt.assignment = ( .Call("clv_optimalAssignment", cnf.mx, PACKAGE="clv") + 1 )

	app = 0	
	opt.asgn.len = length(opt.assignment)
	for( i in 1:opt.asgn.len ) app = app + cnf.mx[i,opt.assignment[i]]
	
	result = (app - 1)/(sum(cnf.mx) - 1)
	return(result)
}
