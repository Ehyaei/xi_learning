library(igraph)


# This function computes the xi correlation coefficient.

xicorln = function (xvec, yvec) {
    n <- length(xvec)
    PI <- rank(xvec, ties.method = "random")
    fr <- rank(yvec, ties.method = "max")/n
    gr <- rank((-yvec), ties.method = "max")/n
    ord <- order(PI)
    fr <- fr[ord]
    A1 <- sum(abs(fr[1:(n - 1)] - fr[2:n]))/(2 * n)
    CU <- mean(gr * (1 - gr))
    xi <- 1 - A1/CU
    return(xi)
}


# This is the main function that computes the tree from data. 
# The input is a matrix x whose rows are the data vectors. 
# The sample size n is the number of rows.
# The number of variables p is the number of columns
# The function outputs the tree g.

condeptree = function(x) {
	p = ncol(x)
	n = nrow(x)
 	d = matrix(0, nrow = p, ncol = p)
	for (i in 1:(p-1)) {
		for (j in (i+1):p) {
			d[i,j] = xicorln(x[,i], x[,j])
			d[j,i] = xicorln(x[,j], x[,i])
		}
	}
	m = d
	maxval = max(d)
	for (i in 1:(p-1)) {
		for (j in (i+1):p) {
			u = d[,i]
			v = d[,j]
			w1 = setdiff(which(u >= d[j,i]), c(i,j))
			w2 = setdiff(which(v >= d[i,j]), c(i,j))
			if (length(intersect(w1,w2)) > 0) {
				m[i,j] = 0
				m[j,i] = 0
			}
			if (m[i,j] != 0 & m[j,i] != 0) {
				u1 = maxval + 1 - m[i,j]
				u2 = maxval + 1 - m[j,i]
				u = max(u1,u2)
				m[i,j] = u
				m[j,i] = u
			} 
		}
	}
	m = graph_from_adjacency_matrix(m, mode = "undirected", weighted = TRUE)
	g = mst(m)
	return(g)
}




# This function generates data from a linear tree of size p.

linear_tree = function(p) {
	x = rep(0,p)
	x[1] = rnorm(1)
	for (i in 2:p) {
		x[i] = (x[i-1] + rnorm(1))/sqrt(2)
	}
	return(x)
}

# This function generates data from n linear trees of size p and returns the data and the reconstructed trees from the condeptree algorithm.

linear_tree_sample = function(n,p) {
	x = matrix(nrow = n, ncol = p)
	for (i in 1:n) {
		x[i,] = linear_tree(p)
	}
	tree = condeptree(x)
	return(list(x = x, tree = tree))
}


# This function tests how accurately linear trees are reconstructed by the condeptree algorithm.

linear_tree_test = function(n,p,iter) {
	s = 0
	for (i in 1:iter) {
		q = linear_tree_sample(n,p)
		tree = q$tree
		edges = get.edgelist(tree)
		u = edges[,1]
		v = edges[,2]
		a = 0
		for (j in 1:(p-1)) {
			w = which((u == j & v == j + 1) | (u == j + 1 & v == j))
			l = length(w)
			a = a + l
		}
		s = s + a/(p-1)
	}
	s = s/iter
	return(s)
}


# The remaining functions do all of the above for the other tree types.

binary_tree = function(p) {
	k = log(p + 1, base = 2)
	k = ceiling(k)
	p = 2^k - 1
	x = rep(0,p)
	x[1] = rnorm(1)
	x[2] = (x[1] + rnorm(1))/sqrt(2)
	x[3] = (x[1] + rnorm(1))/sqrt(2)	
	for (i in 2:(k-1)) {
		for (j in 1:(2^i)) {
			x[2^i - 1 + j] = (x[2^(i-1) - 1 + ceiling(j/2)] + rnorm(1))/sqrt(2)
		}
	}
	return(x)
}

binary_tree_sample = function(n,p) {
	k = log(p + 1, base = 2)
	k = ceiling(k)
	p = 2^k - 1
	x = matrix(nrow = n, ncol = p)
	for (i in 1:n) {
		x[i,] = binary_tree(p)
	}
	tree = condeptree(x)
	return(list(x = x, tree = tree))
}

binary_tree_test = function(n,p,iter) {
	k = log(p + 1, base = 2)
	k = ceiling(k)
	p = 2^k - 1
	s = 0
	true1 = rep(0, p-1)
	true2 = rep(0, p-1)
	true1[1] = 1
	true2[1] = 2
	true1[2] = 1
	true2[2] = 3
	for (i in 2:(k-1)) {
		for (j in 1:(2^i)) {
			true1[2^i - 2 + j] = 2^(i-1) - 1 + ceiling(j/2)
			true2[2^i - 2 + j] = 2^i - 1 + j			
		}
	}	
	for (i in 1:iter) {
		q = binary_tree_sample(n,p)
		tree = q$tree
		edges = get.edgelist(tree)
		u = edges[,1]
		v = edges[,2]
		a = 0
		for (j in 1:(p-1)) {
			w = which((u == true1[j] & v == true2[j]) | (u == true2[j] & v == true1[j]))
			l = length(w)
			a = a + l
		}
		s = s + a/(p-1)
	}
	s = s/iter
	return(s)
}


star_tree = function(p) {
	x = rep(0,p)
	x[1] = rnorm(1)
	for (i in 2:p) {
		x[i] = (x[1] + rnorm(1))/sqrt(2)
	}
	return(x)
}


star_tree_sample = function(n,p) {
	x = matrix(nrow = n, ncol = p)
	for (i in 1:n) {
		x[i,] = star_tree(p)
	}
	tree = condeptree(x)
	return(list(x = x, tree = tree))
}

star_tree_test = function(n,p,iter) {
	s = 0
	true1 = rep(0, p-1)
	true2 = rep(0, p-1)
	for (i in 1:(p-1)) {
		true1[i] = 1
		true2[i] = i+1
	}
	for (i in 1:iter) {
		q = star_tree_sample(n,p)
		tree = q$tree
		edges = get.edgelist(tree)
		u = edges[,1]
		v = edges[,2]
		a = 0
		for (j in 1:(p-1)) {
			w = which((u == true1[j] & v == true2[j]) | (u == true2[j] & v == true1[j]))
			l = length(w)
			a = a + l
		}
		s = s + a/(p-1)
	}
	s = s/iter
	return(s)
}



rev_binary_tree = function(p) {
	k = log(p + 1, base = 2)
	k = ceiling(k)
	p = 2^k - 1
	x = rep(0,p)
	x[(2^(k-1)):(2^k-1)] = rnorm(2^(k-1))
	for (i in rev(0:(k-2))) {
		for (j in 1:(2^i)) {
			x[2^i - 1 + j] = (x[2*(2^i - 1 + j)] + x[2*(2^i - 1 + j) + 1] + rnorm(1))/sqrt(3)
		}
	}
	return(x)
}

rev_binary_tree_sample = function(n,p) {
	k = log(p + 1, base = 2)
	k = ceiling(k)
	p = 2^k - 1
	x = matrix(nrow = n, ncol = p)
	for (i in 1:n) {
		x[i,] = rev_binary_tree(p)
	}
	tree = condeptree(x)
	return(list(x = x, tree = tree))
}

rev_binary_tree_test = function(n,p,iter) {
	k = log(p + 1, base = 2)
	k = ceiling(k)
	p = 2^k - 1
	s = 0
	true1 = rep(0, p-1)
	true2 = rep(0, p-1)
	true1[1] = 1
	true2[1] = 2
	true1[2] = 1
	true2[2] = 3
	for (i in 2:(k-1)) {
		for (j in 1:(2^i)) {
			true1[2^i - 2 + j] = 2^(i-1) - 1 + ceiling(j/2)
			true2[2^i - 2 + j] = 2^i - 1 + j			
		}
	}	
	for (i in 1:iter) {
		q = rev_binary_tree_sample(n,p)
		tree = q$tree
		edges = get.edgelist(tree)
		u = edges[,1]
		v = edges[,2]
		a = 0
		for (j in 1:(p-1)) {
			w = which((u == true1[j] & v == true2[j]) | (u == true2[j] & v == true1[j]))
			l = length(w)
			a = a + l
		}
		s = s + a/(p-1)
	}
	s = s/iter
	return(s)
}

